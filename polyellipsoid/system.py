from ellipsoid import Ellipsoid, Polymer

class System:
    def __init__(
            n_chains,
            chain_lengths,
            bead_mass,
            density,
            axis_length,
            bond_length,
            major_axis=[1,0,0],
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
        self.bead_mass = bead_mass
        self.density = density
        self.axis_length = axis_length
        self.major_axis = major_axis
        self.system_mass = bead_mass * n_chains * chain_lengths
        self.mb_system = None
        
        self.chains = []
        for l in chain_lengths:
            ellipsoid = Ellipsoid(
                    name="bead",
                    mass=self.bead_mass,
                    major_length=self.axis_length,
                    major_axis=self.major_axis,
            )
            chain = Polymer()
            chain.add_bead(
                    bead=ellipsoid, bond_axis="major", separation=bond_length
            )
            chain.build(n=l)
            self.chains.append(chain)
        assert len(self.chains) == len(self.n_chains)

    @mb_system.setter
    def pack(self, box_expand_factor=2):
        """Uses mBuild's fill_box function to fill a cubic box
        with the ellipsoid chains. It may be necessary to expand
        the system volume during the packing step, and handling
        shrinking towards a target density during the simulation.

        Parameters:
        -----------
        box_expand_factor : float, default=2
            The factor by which to expand the box edge lengths during
            the packing step. If PACKMOL fails, you may need to
            increase this parameter.

        """
        self.set_target_box()
        pack_box = self.target_box * box_expand_factor
        system = mb.packing.fill_box(
            compounds=self.chains,
            n_compounds=self.n_chains,
            box=list(pack_box),
            overlap=0.2,
            edge=0.9,
            fix_orientation=True
        )
        self.mb_system = system

    @mb_system.setter
    def stack(self, x, y, n, vector, z_axis_adjust=1.0):
        """Arranges chains in layers on an n x n lattice.

        """
        if len(self.chains) != n*n*2:
            raise ValueError(
                    "Using this method creates a system of n x n "
                    "unit cells with each unit cell containing 2 molecules. "
                    "The number of molecules in the system should equal "
                    "2*n*n."
            )
        next_idx = 0
        system = mb.Compound()
        for i in range(n):
            layer = mb.Compound()
            for j in range(n):
                try:
                    chain1 = self.chains[next_idx]
                    chain2 = self.chains[next_idx + 1]
                    translate_by = np.array(vector)*(x, y, 0)
                    chain2.translate_by(translate_by)
                    cell = mb.Compound(subcompounds=[chain1, chain2])
                    cell.translate((0, y*j, 0))
                    layer.add(cell)
                    next_idx += 2
                except IndexError:
                    pass
            layer.translate((x*i, 0, 0))
            system.add(layer)

        bounding_box = system.get_boundingbox().lengths
        target_z = bounding_box[-1] * z_axis_adjust
        self.mb_system = system
        self.set_target_box(z_constraint=target_z)

    @target_box.setter
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

        Parameters:
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
