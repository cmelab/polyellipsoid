from polyellipsoid import Simulation
from base_test import BaseTest

import hoomd
import numpy as np
import pytest

class TestSimulate(BaseTest):
    def test_change_dt(self, packed_system):
        sim = Simulation(
                system=packed_system,
                lperp=1.0,
                lpar=1.0,
                epsilon=1.0,
                dt=0.001,
                r_cut=2.0,
                bond_k=1000,
                seed=42,
                gsd_write=1000,
                log_write=100
        )
        assert sim.dt == 0.001
        sim.dt = 0.0001
        assert sim.dt == 0.0001

    def test_angles(self, packed_system):
        sim = Simulation(
                system=packed_system,
                lperp=1.0,
                lpar=1.0,
                epsilon=1.0,
                dt=0.001,
                r_cut=2.0,
                bond_k=1000,
                angle_k=200,
                angle_theta=2.0,
                seed=42,
                gsd_write=1000,
                log_write=100
        )
        assert sim.snapshot.angles.types[0] == "B-B-B"
        sim.run_NVT(n_steps=2000, kT=2.0, tau_kt=0.001)

    def test_init_sim(self, packed_system):
        sim = Simulation(
                system=packed_system,
                lperp=1.0,
                lpar=1.0,
                epsilon=1.0,
                dt=0.001,
                r_cut=2.0,
                bond_k=1000,
                seed=42,
                gsd_write=1000,
                log_write=100
        )
        ids = sim.snapshot.particles.typeid
        num_rigid = np.count_nonzero(ids==0)
        assert sim.dt == 0.001
        assert num_rigid == 100

        for i in range(0, num_rigid):
            assert sim.snapshot.particles.types[ids[i]] == "R"
    
    def test_shrink(self, sim_init):
        sim_init.run_shrink(n_steps=2000, kT=1.0, tau_kt=0.01) 
        
    def test_run_NVT(self, sim_init):
        sim_init.run_NVT(n_steps=2000, kT=1.0, tau_kt=0.01) 
        assert isinstance(sim_init.method, hoomd.md.methods.NVT)

    def test_run_NPT(self, sim_init):
        sim_init.run_NPT(
                n_steps=0, kT=1.0, tau_kt=0.01, pressure=0.001, tau_pressure=0.1
        ) 
        assert isinstance(sim_init.method, hoomd.md.methods.NPT)

    def test_run_NVE(self, sim_init):
        sim_init.run_NVE(n_steps=2000) 
        assert isinstance(sim_init.method, hoomd.md.methods.NVE)

    def test_update_method(self, sim_init):
        sim_init.run_NVT(n_steps=0, kT=1.0, tau_kt=0.01)
        assert isinstance(sim_init.method, hoomd.md.methods.NVT)
        sim_init.run_NPT(
                n_steps=0, kT=1.0, tau_kt=0.01, pressure=0.001, tau_pressure=0.1
        )
        assert isinstance(sim_init.method, hoomd.md.methods.NPT)
        sim_init.run_NVE(n_steps=0)
        assert isinstance(sim_init.method, hoomd.md.methods.NVE)
        sim_init.run_langevin(n_steps=0, kT=1.0, alpha=0.5)
        assert isinstance(sim_init.method, hoomd.md.methods.Langevin)

    def test_temperature_ramp(self, sim_init):
        kT_ramp = sim_init.temperature_ramp(
                n_steps=2000, kT_start=0.5, kT_final=1.0, period=100
        )
        sim_init.run_NVT(n_steps=2000, kT=kT_ramp, tau_kt=0.001)

    def test_shrink_to_quench(self, sim_init):
        sim_init.run_shrink(n_steps=2000, kT=1.0, tau_kt=0.001)
        sim_init.run_NVT(n_steps=2000, kT=2.0, tau_kt=0.001)
