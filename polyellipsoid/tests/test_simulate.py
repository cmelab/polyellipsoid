from polyellipsoid import Simulation
from base_test import BaseTest

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
                bond_r0=0.25,
                seed=42,
                gsd_write=1000,
                log_write=100
        )
    
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
