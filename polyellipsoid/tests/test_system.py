from polyellipsoid import System
from base_test import BaseTest

import numpy as np
import pytest

class TestSystem(BaseTest):
    
    def test_n_beads(self, packed_system):
        assert packed_system.n_beads == 25

    def test_rigid_body_count(self, packed_system):
        ids = packed_system.snapshot.particles.typeid
        num_rigid = np.count_nonzero(ids == 0)
        assert num_rigid == 25

        for i in range(0, num_rigid):
            assert packed_system.snapshot.particles.types[ids[i]] == "R"

