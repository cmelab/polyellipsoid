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
        sim_init.shrink(n_steps=1e2, kT=1.0) 
        
    def test_quench(self, sim_init):
        sim_init.quench(n_steps=1e3, kT=1.0) 
        
