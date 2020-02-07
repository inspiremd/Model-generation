from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *

class NewPDBReporter(object):
    """NewPDBReporter outputs a series of frames from a Simulation to a PDB file ignoring all the solvent molecules (water and counterions).

    To use it, create a PDBReporter, then add it to the Simulation's list of reporters.
    """

    def __init__(self, file, reportInterval, enforcePeriodicBox=None):
        """Create a PDBReporter.

        Parameters
        ----------
        file : string
            The file to write to
        reportInterval : int
            The interval (in time steps) at which to write frames
        enforcePeriodicBox: bool
            Specifies whether particle positions should be translated so the center of every molecule
            lies in the same periodic box.  If None (the default), it will automatically decide whether
            to translate molecules based on whether the system being simulated uses periodic boundary
            conditions.
        """
        self._reportInterval = reportInterval
        self._enforcePeriodicBox = enforcePeriodicBox
        self._out = open(file, 'w')
        self._topology = None
        self._nextModel = 0

    def describeNextReport(self, simulation):
        """Get information about the next report this object will generate.

        Parameters
        ----------
        simulation : Simulation
            The Simulation to generate a report for

        Returns
        -------
        tuple
            A six element tuple. The first element is the number of steps
            until the next report. The next four elements specify whether
            that report will require positions, velocities, forces, and
            energies respectively.  The final element specifies whether
            positions should be wrapped to lie in a single periodic box.
        """
        steps = self._reportInterval - simulation.currentStep%self._reportInterval
        return (steps, True, False, False, False, self._enforcePeriodicBox)

    def report(self, simulation, state):
        """Generate a report.

        Parameters
        ----------
        simulation : Simulation
            The Simulation to generate a report for
        state : State
            The current state of the simulation
        """
        if self._nextModel == 0:
            # create a new topology just for the protein .pdb
            # more specifically, we are filtering our here the water and counterion residues
            protein_top = Topology()

            chain = None
            residue = None
            selected_indices = []
            for atom in simulation.topology.atoms():
                # select all but water HOH and counter ions (only the 9 counterions allowed in OpenMM modeller included)
                if atom.residue.name in {'WAT', 'HOH', 'TP4', 'TP5', 'T4E', 'NA', 'K', 'LI', 'CS', 'RB', 'CL', 'BR', 'F', 'IOD'}:
                    continue

                # create a new chain each time for each new chain
                if chain is not atom.residue.chain:
                    chain = protein_top.addChain(atom.residue.chain)
                # create a new residue for each new residue
                if residue is None or residue.id != atom.residue.id:
                    residue = protein_top.addResidue(atom.residue.name, chain)

                # update the protein topology
                protein_top.addAtom(atom.name, atom.element, residue)
                # remember the indices of the selected atoms
                selected_indices.append(atom.index)

            # save our topology and indices
            self._topology = protein_top
            self._selected_indices = selected_indices
            PDBFile.writeHeader(self._topology, self._out)
            self._nextModel += 1

        # get the 3D coordinates
        selected_positions_np = state.getPositions(asNumpy=True)[self._selected_indices]
        # update the .pdb file
        PDBFile.writeModel(self._topology, selected_positions_np, self._out, self._nextModel)
        self._nextModel += 1
        if hasattr(self._out, 'flush') and callable(self._out.flush):
            self._out.flush()

    def __del__(self):
        if self._topology is not None:
            PDBFile.writeFooter(self._topology, self._out)
        self._out.close()

