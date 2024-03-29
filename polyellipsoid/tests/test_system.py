from polyellipsoid import System
from base_test import BaseTest

import hoomd
import numpy as np
import pytest

class TestSystem(BaseTest):
    
    def test_n_beads(self, packed_system):
        assert packed_system.n_beads == 100 

    def test_mass(self, packed_system):
        assert packed_system.system_mass == packed_system.bead_mass * 100 

    def test_bad_n_chains(self):
        with pytest.raises(AssertionError):
            sys = System(
                    n_chains=[5, 10],
                    chain_lengths=5,
                    bead_mass=1000,
                    density=0.5,
                    bond_length=0.25,
                    bead_length=2
            )

        with pytest.raises(AssertionError):
            sys = System(
                    n_chains=5,
                    chain_lengths=[5, 10],
                    bead_mass=1000,
                    density=0.5,
                    bond_length=0.25,
                    bead_length=2
            )


    def test_make_chains(self, packed_system):
        sys = System(
                n_chains=[5, 5],
                chain_lengths=[4, 5],
                bead_mass=1000,
                density=0.5,
                bead_length=2,
                bond_length=0.2
        )
        assert len(sys.chains) == 10 
        assert sys.n_beads == 45


    def test_pack(self):
        sys = System(
                n_chains=5,
                chain_lengths=5,
                bead_mass=1000,
                density=0.5,
                bond_length=0.25,
                bead_length=2
        )
        sys.pack()
        assert isinstance(sys.target_box, np.ndarray)

    def test_stack(self):
        sys = System(
                n_chains=8,
                chain_lengths=5,
                bead_mass=1000,
                density=0.5,
                bond_length=0.25,
                bead_length=2
        )
        sys.stack(z=1, y=1, n=2, vector=[1,1,0])
        assert isinstance(sys.target_box, np.ndarray)

    def test_stack_wrong_num_chains(self):
        sys = System(
                n_chains=7,
                chain_lengths=5,
                bead_mass=1000,
                density=0.5,
                bond_length=0.25,
                bead_length=2
        )
        with pytest.raises(ValueError):
            sys.stack(z=1, y=1, n=2, vector=[1,1,0])
