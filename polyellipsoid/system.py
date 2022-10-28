from polyellipsoid import Ellipsoid
from polyellipsoid.utils import base_units

import hoomd
import mbuild as mb
from mbuild.formats.hoomd_forcefield import to_hoomdsnapshot
from mbuild.lib.recipes.polymer import Polymer
import numpy as np

units = base_units.base_units()


class System:
    """
    """
    def __init__(
            self,
            n_chains,
            chain_lengths,
            bead_length,
            bead_mass,
            density,
            bond_length=.01,
            seed=42,
    ):
        if not isinstance(n_chains, list):
            n_chains = [n_chains]
        if not isinstance(chain_lengths, list):
            chain_lengths = [chain_lengths]
        assert (
                len(n_chains) == len(chain_lengths)
        ), "n_chains and chain_lengths must be equal in length"

        self.n_chains = n_chains
        self.chain_lengths = chain_lengths
        self.density = density
        self.bead_mass = bead_mass
        self.bead_length = bead_length
        self.bond_length = bond_length
        self.n_beads = sum([i*j for i,j in zip(n_chains, chain_lengths)])
        self.system_mass = bead_mass * self.n_beads
        self.target_box = None
        self.mb_system = None
        
        self.chains = []
        for n, l in zip(n_chains, chain_lengths):
            for i in range(n):
                ellipsoid = Ellipsoid(
                        mass=self.bead_mass, length=self.bead_length
                )
                chain = Polymer()
                chain.add_monomer(
                        ellipsoid,
                        indices=[0, 1],
                        orientation=[[1,0,0], [-1,0,0]],
                        replace=False,
                        separation=self.bond_length
                )
                chain.build(n=l, add_hydrogens=False)
                chain.freud_generate_bonds(
                        name_a="B",
                        name_b="B",
                        dmin=self.bead_length/2, 
                        dmax=self.bead_length/2 + bond_length + 0.1 
                )
                self.chains.append(chain)

    def pack(self, box_expand_factor=5):
        """Uses mBuild's fill_box function to fill a cubic box
        with the ellipsoid chains. It may be necessary to expand
        the system volume during the packing step, and handling
        shrinking towards a target density during the simulation.

        Parameters
        ----------
        box_expand_factor : float, default=5
            The factor by which to expand the box edge lengths during
            the packing step. If PACKMOL fails, you may need to
            increase this parameter.

        """
        if self.target_box is None:
            self.set_target_box()
        pack_box = self.target_box * box_expand_factor
        self.mb_system = mb.packing.fill_box(
            compound=self.chains,
            n_compounds=[1 for i in self.chains],
            box=list(pack_box),
            overlap=0.5,
            edge=0.5,
            fix_orientation=True
        )
        self.mb_system.label_rigid_bodies(discrete_bodies="dimer")

    def stack(self, y, z, n, vector, x_axis_adjust=1.0):
        """Arranges chains in layers on an n x n lattice.

        """
        if sum(self.n_chains) != n*n*2:
            raise ValueError(
                    "Using this method creates a system of n x n "
                    "unit cells with each unit cell containing 2 molecules. "
                    "The number of molecules in the system should equal "
                    f"2*n*n. You have {sum(self.n_chains)} number of chains."
            )
        next_idx = 0
        self.mb_system = mb.Compound()
        for i in range(n):
            layer = mb.Compound()
            for j in range(n): # Add chains to the layer along the y dir
                try:
                    chain1 = self.chains[next_idx]
                    chain2 = self.chains[next_idx + 1]
                    translate_by = np.array(vector)*(0, y, z)
                    chain2.translate(translate_by)
                    cell = mb.Compound(subcompounds=[chain1, chain2])
                    cell.translate((0, y*j, 0))
                    layer.add(cell)
                    next_idx += 2
                except IndexError:
                    pass
            layer.translate((0, 0, z*1)) # shift layers along x dir
            self.mb_system.add(layer)

        bounding_box = np.array(self.mb_system.get_boundingbox().lengths)
        bounding_box *= 3 
        target_x = bounding_box[0] * x_axis_adjust
        self.mb_system.box = mb.box.Box(bounding_box)
        self.set_target_box(x_constraint=target_x)
        self.mb_system.translate_to(
                (
                    self.mb_system.box.Lx/2,
                    self.mb_system.box.Ly/2,
                    self.mb_system.box.Lz/2
                )
        )
        self.mb_system.label_rigid_bodies(discrete_bodies="dimer")

    def set_target_box(
            self,
            x_constraint=None,
            y_constraint=None,
            z_constraint=None
    ):
        """Set the target volume of the system during
        the initial shrink step.
        If no constraints are set, the target box is cubic.
        Setting constraints will hold those box vectors
        constant and adjust others to match the target density.

        Parameters
        -----------
        x_constraint : float, optional, defualt=None
            Fixes the box length along the x axis
        y_constraint : float, optional, default=None
            Fixes the box length along the y axis
        z_constraint : float, optional, default=None
            Fixes the box length along the z axis

        """
        if not any([x_constraint, y_constraint, z_constraint]):
            Lx = Ly = Lz = self._calculate_L()
        else:
            constraints = np.array([x_constraint, y_constraint, z_constraint])
            fixed_L = constraints[np.where(constraints!=None)]
            #Conv from nm to cm for _calculate_L
            fixed_L /= units["cm_to_nm"]
            L = self._calculate_L(fixed_L = fixed_L)
            constraints[np.where(constraints==None)] = L
            Lx, Ly, Lz = constraints

        self.target_box = np.array([Lx, Ly, Lz])

    def _calculate_L(self, fixed_L=None):
        """Calculates the required box length(s) given the
        mass of a sytem and the target density.

        Box edge length constraints can be set by set_target_box().
        If constraints are set, this will solve for the required
        lengths of the remaining non-constrained edges to match
        the target density.

        Parameters
        ----------
        fixed_L : np.array, optional, defualt=None
            Array of fixed box lengths to be accounted for
            when solving for L

        """
        M = self.system_mass * units["amu_to_g"]  # grams
        vol = (M / self.density) # cm^3
        if fixed_L is None:
            L = vol**(1/3)
        else:
            L = vol / np.prod(fixed_L)
            if len(fixed_L) == 1: # L is cm^2
                L = L**(1/2)
        L *= units["cm_to_nm"]  # convert cm to nm
        return L
