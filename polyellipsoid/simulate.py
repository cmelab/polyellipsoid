from cmeutils.geometry import moit
from cmeutils.gsd_utils import create_rigid_snapshot, update_rigid_snapshot
from gmso.external.convert_mbuild import from_mbuild
from gmso.external.convert_parmed import to_parmed
from mbuild.formats.hoomd_forcefield import to_hoomdsnapshot

import hoomd
import numpy as np


class Simulation:
    """Simulation initialization class.

    Parameters
    ----------
    system : polyellipsoid.system.System, required
        System containing initial snapshot
    epsilon : float, required
        Energy value for Gay-Berne pair potential
    lperp : float, required
        Perpendicular length for Gay-Berne pair potential
    lpar : float, required
        Parallel length for Gay-Berne pair potential
    bond_k : float, required
        Spring constant for hoomd.md.bond.Harmonic force
    r_cut : float, required
        Cutoff radius for potentials (in simulation distance units)
    dt : float, optional, default 0.001
        Size of simulation timestep (in simulation time units)
    seed : int, optional, default 21
        Seed passed to integrator when randomizing velocities.
    gsd_write : int, default 1e4
        Period to write simulation snapshots to gsd file.
    log_write : int, default 1e3
        Period to write simulation data to the log file.

    Methods
    -------
    run_shrink : Runs a Hoomd simulation
        Run a shrink simulation to reduce simulation volume to match
        the target box set by system.target_box
    run_NPT : Runs a Hoomd simulation in the NPT ensemble
    run_NVT : Runs a Hoomd simulation in the NVT ensemble
    run_langevin : Runs a Hoomd using Langevin dynamics
    run_NVE : Runs a Hoomd simulation in the NVE ensemble
    temperature_ramp : Retruns a hoomd.variant.Ramp()
        Can be passed into the kT parameter for each of the run functions
    set_integrator_method : Sets an initial (or updates) integrator method
        This is called automatically within the run functions, but
        can also be used directly if needed.

    """
    def __init__(
            self,
            system,
            epsilon,
            lperp,
            lpar,
            bond_k,
            r_cut,
            angle_k=None,
            angle_theta=None,
            dt=0.0001,
            seed=21,
            gsd_write=1e4,
            log_write=1e3,
    ):
        self.system = system
        self._dt = dt
        # Snapshot with rigid center placeholders, no toplogy information
        init_snap = create_rigid_snapshot(system.mb_system)
        # Use GMSO to populate angle information before making snapshot
        gmso_system = from_mbuild(system.mb_system)
        gmso_system.identify_connections()
        parmed_system = to_parmed(gmso_system)
        # Atom types need to be set for angles to be correctly added
        for atom in parmed_system.atoms:
            atom.type = atom.name
        # Snapsot with complete toplogy information added
        self._snapshot, refs = to_hoomdsnapshot(
                parmed_system, hoomd_snapshot=init_snap
        )
        # Snapshot with info updated for rigid centers, used by Hoomd
        self.snapshot, self.rigid = update_rigid_snapshot(
                snapshot=self._snapshot, mb_compound=system.mb_system
        )
        self.gsd_write = gsd_write
        self.log_write = log_write
        self.log_quantities = [
            "kinetic_temperature",
            "potential_energy",
            "kinetic_energy",
            "volume",
            "pressure",
            "pressure_tensor"
        ]
        # Set up sim object
        self.device = hoomd.device.auto_select()
        self.sim = hoomd.Simulation(device=self.device, seed=seed)
        self.sim.create_state_from_snapshot(self.snapshot)
        self.forcefield = []
        # Set up forces, GB pair, harmonic bonds and angles:
        nl = hoomd.md.nlist.Cell(buffer=0.40)
        gb = hoomd.md.pair.aniso.GayBerne(nlist=nl, default_r_cut=r_cut)
        gb.params[('R', 'R')] = dict(epsilon=epsilon, lperp=lperp, lpar=lpar)
        zero_pairs = [('A','A'), ('B','B'), ('A','B'), ('A','R'), ('B','R')]
        for pair in zero_pairs:
            gb.params[pair] = dict(epsilon=0.0, lperp=0.0, lpar=0.0)
        self.forcefield.append(gb)

        harmonic_bond = hoomd.md.bond.Harmonic()
        harmonic_bond.params["A-A"] = dict(
                k=bond_k, r0=self.system.bond_length * 10
        )
        harmonic_bond.params["B-B"] = dict(
                k=0, r0=(self.system.bead_length * 10) / 2
        )
        self.forcefield.append(harmonic_bond)

        if all([angle_k, angle_theta]):
            harmonic_angle = hoomd.md.angle.Harmonic()
            harmonic_angle.params["B-B-B"] = dict(k=angle_k, t0=angle_theta)
            self.forcefield.append(harmonic_angle)

        self.all = hoomd.filter.Rigid(("center", "free"))
        # Set up gsd and log writers
        self.integrator = None
        gsd_writer, table_file = self._hoomd_writers()
        self.sim.operations.writers.append(gsd_writer)
        self.sim.operations.writers.append(table_file)

    @property
    def dt(self):
        return self._dt

    @property 
    def method(self):
        return self.sim.operations.integrator.methods[0]

    @dt.setter
    def dt(self, value):
        self._dt = value 
        if self.integrator:
            self.sim.operations.integrator.dt = self._dt

    def set_integrator_method(
            self,
            integrator_method,
            method_kwargs,
    ):
        # No integrator and method has been created yet
        if not self.integrator:
            self.integrator = hoomd.md.Integrator(
                    dt=self.dt, integrate_rotational_dof=True
            )
            self.integrator.rigid = self.rigid
            self.integrator.forces = self.forcefield
            self.sim.operations.add(self.integrator)
            new_method = integrator_method(**method_kwargs) 
            self.sim.operations.integrator.methods = [new_method]
        # Update the existing integrator with a new method
        else:
            self._update_integrator_method(integrator_method, method_kwargs)

    def _update_integrator_method(self, integrator_method, method_kwargs):
        self.integrator.methods.remove(self.method)
        new_method = integrator_method(**method_kwargs)
        self.integrator.methods.append(new_method)

    def run_shrink(
            self,
            kT,
            tau_kt,
            n_steps,
            shrink_period=10,
            thermalize_particles=True
    ):
        """Run a shrink simulation to reach a target volume."""
        # Set up box resizer
        box_resize_trigger = hoomd.trigger.Periodic(shrink_period)
        box_ramp = hoomd.variant.Ramp(
                A=0, B=1, t_start=0, t_ramp=int(n_steps)
        )
        # Convert from nm (mbuild units, used in System()) to Angstrom
        self.target_box = hoomd.Box(
                Lx=self.system.target_box[0] * 10,
                Ly=self.system.target_box[1] * 10,
                Lz=self.system.target_box[2] * 10,
        )
        box_resize=hoomd.update.BoxResize(
                box1=self.sim.state.box,
                box2=self.target_box,
                variant=box_ramp,
                trigger=box_resize_trigger,
                filter=self.all
        )
        self.sim.operations.updaters.append(box_resize)
        # Use NVT integrator during shrinking
        self.set_integrator_method(
                integrator_method=hoomd.md.methods.NVT,
                method_kwargs={"tau": tau_kt, "filter": self.all, "kT": kT},
        )
        if thermalize_particles:
            if isinstance(kT, hoomd.variant.Ramp):
                self.sim.state.thermalize_particle_momenta(
                        filter=self.all, kT=kT.range[0]
                )
            else:
                self.sim.state.thermalize_particle_momenta(
                        filter=self.all, kT=kT
                )
        self.sim.run(n_steps) 

    def run_langevin(
            self,
            n_steps,
            kT,
            alpha,
            tally_reservoir_energy=False,
            default_gamma=1.0,
            default_gamma_r=(1.0, 1.0, 1.0),
            thermalize_particles=True
    ):
        self.set_integrator_method(
                integrator_method=hoomd.md.methods.Langevin,
                method_kwargs={
                        "filter": self.all,
                        "kT": kT,
                        "alpha": alpha,
                        "tally_reservoir_energy": tally_reservoir_energy,
                        "default_gamma": default_gamma,
                        "default_gamma_r": default_gamma_r,
                    }
        )
        if thermalize_particles:
            if isinstance(kT, hoomd.variant.Ramp):
                self.sim.state.thermalize_particle_momenta(
                        filter=self.all, kT=kT.range[0]
                )
            else:
                self.sim.state.thermalize_particle_momenta(
                        filter=self.all, kT=kT
                )
        self.sim.run(n_steps)

    def run_NPT(
            self,
            n_steps,
            kT,
            pressure,
            tau_kt,
            tau_pressure,
            couple="xyz",
            box_dof=[True, True, True, False, False, False],
            rescale_all=False,
            gamma=0.0,
            thermalize_particles=True
    ):
        self.set_integrator_method(
                integrator_method=hoomd.md.methods.NPT,
                method_kwargs={
                    "kT": kT,
                    "S": pressure,
                    "tau": tau_kt,
                    "tauS": tau_pressure,
                    "couple": couple,
                    "box_dof": box_dof,
                    "rescale_all": rescale_all,
                    "gamma": gamma,
                    "filter": self.all,
                    "kT": kT
                }
        )
        if thermalize_particles:
            if isinstance(kT, hoomd.variant.Ramp):
                self.sim.state.thermalize_particle_momenta(
                        filter=self.all, kT=kT.range[0]
                )
            else:
                self.sim.state.thermalize_particle_momenta(
                        filter=self.all, kT=kT
                )
            self.sim.run(0)
            self.method.thermalize_thermostat_and_barostat_dof()
        self.sim.run(n_steps)
    
    def run_NVT(self, n_steps, kT, tau_kt, thermalize_particles=True):
        self.set_integrator_method(
                integrator_method=hoomd.md.methods.NVT,
                method_kwargs={"tau": tau_kt, "filter": self.all, "kT": kT},
        )
        if thermalize_particles:
            if isinstance(kT, hoomd.variant.Ramp):
                self.sim.state.thermalize_particle_momenta(
                        filter=self.all, kT=kT.range[0]
                )
            else:
                self.sim.state.thermalize_particle_momenta(
                        filter=self.all, kT=kT
                )
            self.sim.run(0)
            self.method.thermalize_thermostat_dof()
        self.sim.run(n_steps)

    def run_NVE(self, n_steps):
        self.set_integrator_method(
                integrator_method=hoomd.md.methods.NVE,
                method_kwargs={"filter": self.all},
        )
        self.sim.run(n_steps)

    def temperature_ramp(
            self,
            n_steps,
            kT_start,
            kT_final,
            period,
    ):  
        return hoomd.variant.Ramp(
                A=kT_start,
                B=kT_final,
                t_start=self.sim.timestep,
                t_ramp=int(n_steps)
        )

    def _hoomd_writers(self):
        """Creates gsd and log writers"""
        writemode = "w"
        gsd_writer = hoomd.write.GSD(
                filename="sim_traj.gsd",
                trigger=hoomd.trigger.Periodic(int(self.gsd_write)),
                mode=f"{writemode}b",
                dynamic=["momentum"]
        )

        logger = hoomd.logging.Logger(categories=["scalar", "string"])
        logger.add(self.sim, quantities=["timestep", "tps"])
        thermo_props = hoomd.md.compute.ThermodynamicQuantities(filter=self.all)
        self.sim.operations.computes.append(thermo_props)
        logger.add(thermo_props, quantities=self.log_quantities)
        for f in self.forcefield:
            logger.add(f, quantities=["energy"])

        table_file = hoomd.write.Table(
            output=open("sim_traj.txt", mode=f"{writemode}", newline="\n"),
            trigger=hoomd.trigger.Periodic(period=int(self.log_write)),
            logger=logger,
            max_header_len=None,
        )
        return gsd_writer, table_file
