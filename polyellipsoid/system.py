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


    def pack(self, box_expand_factor):
        system = mb.packing.fill_box(
            compounds=self.chains, n_compounds=self.n_chains
        )

    def stack(self, a, b, n, vector):
	    pass

    def _set_target_box(
            self,
            x_constraint=None,
            y_constraint=None,
            z_constraint=None
    ):
        """Set the target volume of the system during
        the initial shrink step.
        If no constraints are set, the target box is cubic.
        Setting constraints will hold certain box vectors
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
