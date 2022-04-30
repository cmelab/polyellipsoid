import hoomd

class Simulation:
    def __init__(
            system,
            epsilon,
            lperp,
            lpar,
            tau,
            dt,
            r_cut,
            bond_k,
            bond_r0,
            integrator,
            seed,
            gsd_write,
            log_write,
    ):
        self.system = system
        self.snapshot = system.snapshot
        self.epsilon = epsilon
        self.lperp = lperp
        self.lpar = lpar
        self.tau = tau
        self.dt = dt
        self.r_cut = r_cut
        self.bond_k = bond_k
        self.bond_r0 = bond_r0
        self.seed = seed
        self.gsd_write = gsd_write
        self.log_write = log_write
        self.ran_shrink = False
        
        # Set up rigid object for Hoomd
        inds = [system.n_beads, system.n_beads + 1]
        rigid = hoomd.md.constrain.Rigid()
        r_pos = self.snapshot.particles.position[0]
        c_pos = self.snapshot.particles.position[inds]
        c_pos -= r_pos
        c_pos = [tuple(i) for i in c_pos]
        c_types = [
                self.snapshot.particles.types[i] for
                i in self.snapshot.particles.typeid[inds]
        ]
        c_orient = [tuple(i) for i in self.snapshot.particles.orientaiton[inds]]
        c_charge = [i for i in self.snapshot.particles.charge[inds]]
        c_diam = [i for i in self.snapshot.particles.diameter[inds]]
        mass = np.array([self.system.bead_mass]*2)
        moits = moit(c_pos, mass)
        self.snapshot.particles.moment_inertia[0:system.n_beads] = moits

        rigid.body["R"] = {
                "constituent_types": c_types,
                "positions": c_pos,
                "orientations": c_orient,
                "charges": c_charge,
                "diameters": c_diam
        }
        
        # Set up forces, GB pair and harmonic bond:
        nl = hoomd.nlist.Cell(buffer=0.40)
        gb = hoomd.md.pair.aniso.GayBerne(nlist=nl, default_r_cut=r_cut)
        gb.params[('R', 'R')] = dict(epsilon=epsilon, lperp=lperp, lpar=lpar)
        zero_pairs = [
                ('CT','CT'), ('CH','CT'), ('CH','CH'), ('CT','R'), ('CH','R')
        ]
        for pair in zero_pairs:
            gb.params[pair] = dict(epsilon=0.0, lperp=0.0, lpar=0.0)
        
        # Set up harmonic bond force
        harmonic = hoomd.md.bond.Harmonic()
        harmonic.params["CT-CH"] = dict(k=bond_k, r0=bond_r0)
        harmonic.params["CH-CT"] = dict(k=bond_k, r0=bond_r0)

        # Set up remaining Hoomd objects
        self.centers = hoomd.filter.Rigid()
        #TODO: Do we need the _all filter?
        _all = hoomd.filter.All()
        self.integrator = hoomd.md.Integrator(
                dt=dt, integrate_rotational_dof=True
        )
        self.integrator.rigid=rigid
        self.integrator.forces = [gb, harmonic]
        # Set up sim object
        self.sim = hoomd.Simulation(
                device=hoomd.device.auto_select(), seed=seed
        )
        self.sim.create_state_from_snapshot(self.snapshot)

    def shrink(self, kT, tau, n_steps, shrink_period):
        # Set up box resizer
        box_resize_trigger = hoomd.trigger.Periodic(shrink_period)
        ramp = hoomd.variant.Ramp(
                A=0, B=1, t_start=0, t_ramp=int(n_steps)
        )
        target_box = hoomd.Box(
                Lx=self.system.target_box[0],
                Ly=self.system.target_box[1],
                Lz=self.system.target_box[2],
        )
        box_resize=hoomd.update.BoxResize(
                box1=sim.state.box,
                box2=target_box,
                variant=ramp,
                trigger=box_resize_trigger
        )
        self.sim.operations.updates.append(box_resize)
        # Use NVT integrator during shrinking
        integrator_method = hoomd.md.methods.NVT(
                filter=self.centers, kT=kT, tau=tau
        )
        self.integrator.methods = [integrator_method]
        self.sim.operations.add(self.integrator)
        self.sim.state.thermalize_particle_momenta(filter=self.centers, kT=kT)
        self.sim.run(n_steps + 1) 
        self.ran_shrink = True

    def quench(self, kT, tau, n_steps):
        if self.ran_shrink: # Shrink step ran, update temperature & tau
            self.integrator.methods[0].kT = kT
            self.integrator.methods[0].tau = tau
        else: # Shrink not ran, add integrator method
            integrator_method = hoomd.md.methods.NVT(
                    filter=self.centers, kT=kT, tau=tau
            )
            self.integrator.methods = [integrator_method]
            self.sim.operations.add(self.integrator)

        self.sim.state.thermalize_particle_momenta(filter=self.centres, kT=kT)
        sim.run(n_steps)

    def anneal(self, kT_init, kT_final, step_sequence, schedule):
        pass

