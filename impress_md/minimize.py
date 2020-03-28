from simtk.openmm import app
import simtk.openmm as mm
from simtk import unit
import sys
sys.path.append('./impress_md')
import solventlessPdbReporter as nosol

def MinimizedEnergy(filepath):
    prmtop = app.AmberPrmtopFile(f'{filepath}.prmtop')
    inpcrd = app.AmberInpcrdFile(f'{filepath}.inpcrd')
    system = prmtop.createSystem(implicitSolvent=app.GBn2,
                                 nonbondedMethod=app.CutoffNonPeriodic,
                                 nonbondedCutoff=1.0*unit.nanometers,
                                 constraints=app.HBonds,
                                 rigidWater=True,
                                 ewaldErrorTolerance=0.0005)

    integrator = mm.LangevinIntegrator(300*unit.kelvin,
                                       1.0/unit.picoseconds,
                                       2.0*unit.femtoseconds)
    integrator.setConstraintTolerance(0.00001)
    # TODO: This should just recognize whatever the computer is capable of, not force CUDA.
    platform = mm.Platform.getPlatformByName('CUDA')
    # TODO: I am not sure if mixed precision is necessary. It dramatically changes the results.
    properties = {'CudaPrecision': 'mixed'}
    
    simulation = app.Simulation(prmtop.topology, system, integrator, platform)
    simulation.context.setPositions(inpcrd.positions)
    
    simulation.minimizeEnergy()
    energy = simulation.context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilojoule/unit.mole)
    return energy


# def MinimizedEnergyWithParam(filepath):
#     prmtop = app.AmberPrmtopFile(f'{filepath}.prmtop')
#     inpcrd = app.AmberInpcrdFile(f'{filepath}.inpcrd')
#     system = prmtop.createSystem(implicitSolvent=app.GBn2,
#                                  nonbondedMethod=app.CutoffNonPeriodic,
#                                  nonbondedCutoff=1.0 * unit.nanometers,
#                                  constraints=app.HBonds,
#                                  rigidWater=True,
#                                  ewaldErrorTolerance=0.0005)
#
#     integrator = mm.LangevinIntegrator(300 * unit.kelvin,
#                                        1.0 / unit.picoseconds,
#                                        2.0 * unit.femtoseconds)
#     integrator.setConstraintTolerance(0.00001)
#     # TODO: This should just recognize whatever the computer is capable of, not force CUDA.
#     platform = mm.Platform.getPlatformByName('CUDA')
#     # TODO: I am not sure if mixed precision is necessary. It dramatically changes the results.
#     properties = {'CudaPrecision': 'mixed'}
#
#     simulation = app.Simulation(prmtop.topology, system, integrator, platform)
#     simulation.context.setPositions(inpcrd.positions)
#
#     simulation.minimizeEnergy()
#     energy = simulation.context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilojoule / unit.mole)
#     return energy

def simulation(filepath, outpath, nsteps):
    prmtop = app.AmberPrmtopFile(f'{filepath}.prmtop')
    inpcrd = app.AmberInpcrdFile(f'{filepath}.inpcrd')
    forcefield = app.ForceField('amber14-all.xml', 'amber14/tip3p.xml')
    modeller = app.Modeller(prmtop.topology, inpcrd.positions)
    modeller.addSolvent(forcefield, padding=1.4*unit.nanometer)
    system = forcefield.createSystem(modeller.topology, nonbondedMethod=app.PME, nonbondedCutoff=1.0*unit.nanometer,
            constraints=app.HBonds)
    integrator = mm.LangevinIntegrator(300*unit.kelvin, 1.0/unit.picosecond, 0.002*unit.picosecond)
    platform = mm.Platform.getPlatformByName('CUDA')
    properties = {'Precision': 'double'}
    simulation = app.Simulation(modeller.topology, system, integrator, platform, properties)
    simulation.context.setPositions(modeller.positions)
    simulation.minimizeEnergy()
    if nsteps != 0:
        simulation.reporters.append(app.DCDReporter(f'{outpath}_traj.dcd', 25000)) # snapshot at every 50 ps 
        simulation.reporters.append(nosol.NewPDBReporter(f'{outpath}_system_nosol.pdb', 25000)) # snapshot at every 50 ps 
        simulation.reporters.append(app.StateDataReporter(f'{outpath}_sim.log', 25000, step=True,
        potentialEnergy=True, temperature=True)) # reporting at every 50 ps
        simulation.reporters.append(app.CheckpointReporter(f'{outpath}_traj.chk', 250000)) # checkpoint at every 0.5 ns
        simulation.step(nsteps)
        positions = simulation.context.getState(getPositions=True).getPositions()
        app.PDBFile.writeFile(simulation.topology, positions, open(f'{outpath}_output.pdb', 'w'))

# Return potential energy at the end of the simulation
#    potential = simulation.context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilojoule/unit.mole)
#    return potential

# Return mean potential energy during the simulation
    with open(f'{outpath}_sim.log','r') as log_file: 
        lines = log_file.readlines()
        mean = 0
    for i in range(1,len(lines)):
        potential = float(lines[i].split(',')[1])
        mean = mean + potential
    potential = mean/(len(lines)-1)
    return potential

def simulation_ESMACS(filepath, outpath, nsteps):
    prmtop = app.AmberPrmtopFile(f'{filepath}_sol.prmtop')
    inpcrd = app.AmberInpcrdFile(f'{filepath}_sol.inpcrd')
    system = prmtop.createSystem(nonbondedMethod=app.PME, nonbondedCutoff=1.0*unit.nanometer,
            constraints=app.HBonds)
    integrator = mm.LangevinIntegrator(300*unit.kelvin, 1.0/unit.picosecond, 0.002*unit.picosecond)
    platform = mm.Platform.getPlatformByName('CUDA')
    properties = {'Precision': 'double'}
    simulation = app.Simulation(prmtop.topology, system, integrator, platform, properties)
    simulation.context.setPositions(inpcrd.positions)
    if inpcrd.boxVectors is not None:
        simulation.context.setPeriodicBoxVectors(*inpcrd.boxVectors)
    simulation.minimizeEnergy()
    simulation.reporters.append(app.DCDReporter(f'{outpath}/traj.dcd', 50000)) # snapshot at every 100 ps 
    simulation.reporters.append(app.StateDataReporter(f'{outpath}/sim.log', 5000, step=True,
    potentialEnergy=True, temperature=True)) # reporting at every 10 ps
    simulation.reporters.append(app.CheckpointReporter(f'{outpath}/traj.chk', 250000)) # checkpoint at every 0.5 ns
    simulation.step(nsteps)
    positions = simulation.context.getState(getPositions=True).getPositions()
    app.PDBFile.writeFile(simulation.topology, positions, open(f'{outpath}/output.pdb', 'w'))

def MMGBSA(inpath, tmppath, comp):
    prmtop = app.AmberPrmtopFile(f'{inpath}/{comp}.prmtop')
    system = prmtop.createSystem(implicitSolvent=app.GBn2, nonbondedMethod=app.CutoffNonPeriodic, 
            nonbondedCutoff=1.0*unit.nanometer, constraints=app.HBonds)
    integrator = mm.LangevinIntegrator(300*unit.kelvin, 1.0/unit.picosecond, 0.002*unit.picosecond)
    platform = mm.Platform.getPlatformByName('CUDA')
    properties = {'Precision': 'double'}
    simulation = app.Simulation(prmtop.topology, system, integrator, platform, properties)
    
    import os, glob
    mean = 0
    cnt = 0
    for file in glob.glob(os.path.join(tmppath, comp + '_MMGBSA_snap*.pdb')):
        pdb = app.PDBFile(file)
        simulation.context.setPositions(pdb.positions)
        potential = simulation.context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilojoule/unit.mole)
        mean = mean + potential
        cnt += 1
    potential = mean/cnt
    return potential

	
