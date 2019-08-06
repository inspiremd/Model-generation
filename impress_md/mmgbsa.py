import numpy as np
from pymbar import timeseries
from math import sqrt

def simulate(inpcrd_filenames, prmtop_filenames, niterations=1000, implicit=True):
    """
    Run simple MM-PBSA like workflow to estimate free energy of binding from enthalpies of interaction.
    This method can accept either explicit or implicit solvent.
    .. todo :: Allow user to leave out 'receptor' if desired.

    Parameters
    ----------
    inpcrd_filename : dict of str
        inpcrd_filename[phase] is the Amber inpcrd filename for phase in 'complex', 'ligand', and (optionally) 'receptor'
    prmtop_filenames : dict of str
        prmtop_filename[phase] is the Amber prmtop filename for phase in 'complex', 'ligand', and (optionally) 'receptor'
    niterations : int, optional, default=100
        Maximum number of iterations to run.

    Returns
    -------
    enthalpies: float array
        A sample of GB enthalpies from the trajectory
    """
    from simtk import unit, openmm
    from simtk.openmm import app

    temperature = 300 * unit.kelvin
    collision_rate = 1.0 / unit.picoseconds
    hydrogen_mass = 3.5 * unit.amu
    timestep = 4.0 * unit.femtoseconds
    pressure = 1.0 * unit.atmospheres

    if implicit:
        method = app.CutoffNonPeriodic
    else:
        method = app.PME

    # Store thermal energy
    kT = unit.BOLTZMANN_CONSTANT_kB * temperature

    # TODO: Autodetect nonbondedMethod

    # Termination criteria
    nsteps_per_iteration = 250 # 1 picosecond
    # TODO: Use statistical error threshold to drive this to convergence
    phases = inpcrd_filenames.keys()

    enthalpies = dict()
    for phase in phases:
        # Load Amber system
        prmtop = app.AmberPrmtopFile(prmtop_filenames[phase])
        inpcrd = app.AmberInpcrdFile(inpcrd_filenames[phase])

        if not implicit:
            box = inpcrd.getBoxVectors()
            context.setPeriodicBoxVectors(a=box[0], b=box[1], c=box[2])

        cutoff = 2*unit.nanometers

        # Create system
        # Implicit solvent
        system = prmtop.createSystem(implicitSolvent=app.GBn2,
                                     nonbondedMethod=method,
                                     constraints=app.HBonds,
                                     nonbondedCutoff=cutoff,
                                     hydrogenMass=hydrogen_mass)

        # Create integrator and context
        integrator = openmm.LangevinIntegrator(temperature, collision_rate, timestep)
        context = openmm.Context(system, integrator)
        context.setPositions(inpcrd.getPositions())
        if not implicit:
            context.setPeriodicBoxVectors(inpcrd.getPeriodicBoxVectors())

        # Sample and store enthalpies
        enthalpies[phase] = np.zeros([niterations])

        openmm.LocalEnergyMinimizer.minimize(context)

        for iteration in range(niterations):
            # print("On iteration",iteration,"of",niterations)
            integrator.step(nsteps_per_iteration)
            state = context.getState(getEnergy=True)
            potential_energy = state.getPotentialEnergy()

            if implicit:
                enthalpies[phase][iteration] = potential_energy.value_in_unit(unit.kilojoules_per_mole)
            else:
                volume = state.getPeriodicBoxVolume()
                enthalpy = potential_energy + (pressure * volume).in_unit_system(unit.si_unit_system) / unit.mole
                enthalpies[phase][iteration] = enthalpy.value_in_unit(unit.kilojoules_per_mole)
                #enthalpies[phase][iteration] = (potential_energy + pressure*volume) / kT
    # [mean_dG, var_dG] = mmpbsa(enthalpies)
    return enthalpies

def subsample(enthalpies):
    """
    Subsamples the enthalpies using John Chodera's code.
    I want to compare results to a cutoff
        Returns a smaller array
    """
    for phase in enthalpies:
        [t0, g, Neff_max] = timeseries.detectEquilibration(enthalpies[phase])
        enthalpies[phase] = enthalpies[phase][t0:]
        indices = timeseries.subsampleCorrelatedData(enthalpies[phase], g=g)
        enthalpies[phase] = enthalpies[phase][indices]


def mmgbsa(enthalpies):
    """
    DeltaG, errDeltaG : float
        Estimated free energy of binding (in kT)
    TODO: it calculates standard error rather than variance
    """
    # Use automatic equilibration detection and pymbar.timeseries to subsample
    DeltaH = dict()
    varDeltaH = dict()
    errDeltaH = dict()
    for phase in enthalpies:
        DeltaH[phase] = enthalpies[phase].mean()
        varDeltaH[phase] = enthalpies[phase].std()**2
        errDeltaH[phase] = varDeltaH[phase]/len(enthalpies[phase])

    # Shouldn't have to use names?
    try:
        DeltaH['diff'] = 2*DeltaH['complex']
    except:
        DeltaH['diff'] = 2*DeltaH['com']
    varDeltaH['diff'] = 0
    errDeltaH['diff'] = 0
    for phase in enthalpies:
        DeltaH['diff'] -= DeltaH[phase]
        varDeltaH['diff'] += varDeltaH[phase]
        errDeltaH['diff'] += errDeltaH[phase]

    errDeltaH['diff'] = sqrt(errDeltaH['diff'])
    return DeltaH, errDeltaH

