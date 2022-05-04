import os

import pytest
from polyellipsoid import Ellipsoid, Polymer, System

class BaseTest:
    @pytest.fixture(autouse=True)
    def initdir(self, tmpdir):
        tmpdir.chdir()
    
    @pytest.fixture
    def bead(self):
        return Ellipsoid("A", 100, 2)

    @pytest.fixture
    def packed_system(self):
        sys = System(
                n_chains=5,
                chain_lengths=5,
                bead_mass=1000,
                density=0.5,
                bond_length=0.25,
                axis_length=2
        )
        sys.pack()
        return sys
