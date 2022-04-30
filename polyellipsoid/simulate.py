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
        c_orient = [tuple(i) for i in snapshot.particles.orientaiton[inds]]
        c_charge = [i for i in snapshot.particles.charge[inds]]
        c_diam = [i for i in snapshot.particles.diameter[inds]]
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
        centers = hoomd.filter.Rigid()
        _all = hoomd.filter.All()
        integrator = hoomd.md.Integrator(dt=dt, integrate_rotational_dof=True)
        integrator.rigid=rigid
         



