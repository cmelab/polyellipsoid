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

