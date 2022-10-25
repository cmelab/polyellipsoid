from polyellipsoid import Simulation
from base_test import BaseTest

import numpy as np
import pytest

class TestSimulate(BaseTest):

    def test_init_sim(self, packed_system):
        sim = Simulation(
                system=packed_system,
                lperp=1.0,
                lpar=1.0,
                epsilon=1.0,
                tau=0.1,
                dt=0.001,
                r_cut=2.0,
                bond_k=1000,
                seed=42,
                gsd_write=1000,
                log_write=100
        )
        ids = sim.snapshot.particles.typeid
        num_rigid = np.count_nonzero(ids==0)
        assert num_rigid == 100

        for i in range(0, num_rigid):
            assert sim.snapshot.particles.types[ids[i]] == "R"
    
    def test_shrink(self, sim_init):
        sim_init.shrink(n_steps=2000, kT=1.0) 
        
    def test_quench(self, sim_init):
        sim_init.quench(n_steps=2000, kT=1.0) 
        
    def test_anneal(self, sim_init):
        sim_init.anneal(kT_init=1.0, kT_final=2.0, step_sequence=[200]*5)

    @pytest.mark.skip(reason="Getting particle out of box errors")
    def test_shrink_to_quench(self, sim_init):
        sim_init.shrink(n_steps=2000, kT=1.0)
        sim_init.quench(n_steps=2000, kT=2.0)

    @pytest.mark.skip(reason="Getting particle out of box errors")
    def test_shrink_to_anneal(self, sim_init):
        sim_init.shrink(n_steps=2000, kT=1.0)
        sim_init.anneal(kT_init=1.0, kT_final=2.0, step_sequence=[200]*5)
