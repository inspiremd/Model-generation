from simtk.openmm import app
import simtk.openmm as mm
from simtk import unit

def AddAlchemyForces(system, ligand_ind):
    """
    Input an OpenMM 'system' object and the indices of the ligand
    No output. Function adds alchemical nonbonded forces to the system
    """
    forces = {force.__class__.__name__ : force for force in system.getForces()}
    nbforce = forces['NonbondedForce']
    
    ligand  = ligand_ind
    protein = set(range(system.getNumParticles())) - ligand
    
    alchemical_energy  = 'lambda*4*epsilon*x*(x-1.0); x = (sigma/reff_sterics)^6;'
    alchemical_energy += 'reff_sterics = sigma*(0.5*(1.0-lambda) + (r/sigma)^6)^(1/6);'
    alchemical_energy += 'sigma = 0.5*(sigma1+sigma2); epsilon = sqrt(epsilon1*epsilon2);'
    
    alchemical_force = mm.CustomNonbondedForce(alchemical_energy)
    alchemical_force.addGlobalParameter('lambda',1.0)
    alchemical_force.addPerParticleParameter('sigma')
    alchemical_force.addPerParticleParameter('epsilon')
    for atm in range(system.getNumParticles()):
        # Get the atom's default nonbonded parameters
        [charge, sigma, epsilon] = nbforce.getParticleParameters(atm)
        alchemical_force.addParticle([sigma,epsilon])
        if atm in ligand:
            # TODO: Not actually sure which '*0' options are necessary
            nbforce.setParticleParameters(atm, charge*0, sigma, epsilon)
    alchemical_force.addInteractionGroup(ligand,protein)
    system.addForce(alchemical_force)    
    
def SimulateAlchemy(path, niter, nsteps_per_iter, nlambda):
    """Calculates the binding free energy of a ligand names 'UNL' using alchemy.
    One step corresponds to two femtoseconds.
    """
    prmtop = app.AmberPrmtopFile(f'{path}/com.prmtop')
    inpcrd = app.AmberInpcrdFile(f'{path}/com.inpcrd')
    system = prmtop.createSystem(implicitSolvent=app.GBn2,
                                 nonbondedMethod=app.CutoffNonPeriodic,
                                 nonbondedCutoff=1.0*unit.nanometers,
                                 constraints=app.HBonds,
                                 rigidWater=True,
                                 ewaldErrorTolerance=0.0005)
    
    # Detect ligand indices
    ligand_ind = []
    for atm in prmtop.topology.atoms():
        # OpenEye make the ligand name 'UNL'
        if atm.residue.name == 'UNL':
            ligand_ind.append(atm.index)
    ligand_ind = set(ligand_ind)
    AddAlchemyForces(system, ligand_ind)
    
    integrator = mm.LangevinIntegrator(300*unit.kelvin,
                                       1.0/unit.picoseconds,
                                       2.0*unit.femtoseconds)
    integrator.setConstraintTolerance(0.00001)
    # TODO: The issues here are the same as the mmgbsa.py script
    # TODO: This should just recognize whatever the computer is capable of, not force CUDA.
    # TODO: I am not sure if mixed precision is necessary. Just need to be consistent
    platform = mm.Platform.getPlatformByName('CUDA')
    properties = {'CudaPrecision': 'mixed'}
    simulation = app.Simulation(prmtop.topology, system, integrator, platform)
    simulation.context.setPositions(inpcrd.positions)
    simulation.minimizeEnergy()

    ### Now simulate system
    import numpy as np
    from pymbar import MBAR, timeseries
    lambdas = np.linspace(1.0, 0.0, nlambda)
    # Save the potential energies for MBAR
    u_kln = np.zeros([nlambda, nlambda, niter])
    kT = unit.AVOGADRO_CONSTANT_NA * unit.BOLTZMANN_CONSTANT_kB * integrator.getTemperature()
    # TODO: This runs in series. Someone comfortable with MPI should help parallelize this.
    for k in range(nlambda):
        for i in range(niter):
            print('state %5d iteration %5d / %5d' % (k, i, niter))
            simulation.context.setParameter('lambda',lambdas[k])
            integrator.step(nsteps_per_iter)
            for l in range(nlambda):
                simulation.context.setParameter('lambda',lambdas[l])
                u_kln[k,l,i] = simulation.context.getState(getEnergy=True).getPotentialEnergy() / kT

    # Subsample to reduce variation
    N_k = np.zeros([nlambda], np.int32) # number of uncorrelated samples
    for k in range(nlambda):
        [t0, g, Neff_max] = timeseries.detectEquilibration(u_kln[k,k,:])
        # TODO: maybe should use 't0:' instead of ':' in third index
        indices = timeseries.subsampleCorrelatedData(u_kln[k,k,:], g=g)
        N_k[k] = len(indices)
        u_kln[k,:,0:N_k[k]] = u_kln[k,:,indices].T
    # Calculate the energy difference
    # TODO: I've never worked with pymbar beyond the timeseries function. I'm not sure how the error in DeltaF is calculated, and I don't know what Theta is right now.
    mbar = MBAR(u_kln,N_k)
    [DeltaF_ij, dDeltaF_ij, Theta_ij] = mbar.getFreeEnergyDifferences()
    return DeltaF_ij[0][-1], dDeltaF_ij[0][-1]


