import numpy as np
from pymbar import timeseries
from math import sqrt


def simulate(inpcrd_filenames, prmtop_filenames, niterations=1000, implicit=True):
    from simtk.openmm import app
    import simtk.openmm as mm
    from simtk import unit
    
    phases = inpcrd_filenames.keys()
    nsteps_per_iteration = 500 # 1 picosecond
    
    enthalpies = dict()
    for phase in phases:
        prmtop = app.AmberPrmtopFile(prmtop_filenames[phase])
        inpcrd = app.AmberInpcrdFile(inpcrd_filenames[phase])
        
        system = prmtop.createSystem(implicitSovlent=app.GBn2, 
                 nonbondedMethod=app.CutoffNonPeriodic,
                 nonbondedCutoff=2.0*unit.nanometers, 
                 constraints=app.HBonds)
        integrator = mm.LangevinIntegrator(300*unit.kelvin, 1.0/unit.picoseconds, 2.0*unit.femtoseconds)
        integrator.setConstraintTolerance(0.00001)
        
        platform = mm.Platform.getPlatformByName('CUDA')
        properties = {'CudaPrecision': 'mixed'}
        simulation = app.Simulation(prmtop.topology, system, integrator, platform, properties)
        simulation.context.setPositions(inpcrd.positions)
        
        # Minimize & equilibrate
        simulation.minimizeEnergy()
        simulation.context.setVelocitiesToTemperature(300*unit.kelvin)
        simulation.step(100)
        
        # Run simulation
        for iteration in range(niterations):
            simulation.step(nsteps_per_iteration)
            state = simulation.context.getState(getEnergy=True)
            potential_energy = state.getPotentialEnergy()
            enthalpies[phase][iteration] = potential_energy.value_in_units(unit.kilojoules_per_mole)
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

