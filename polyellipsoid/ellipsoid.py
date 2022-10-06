from mbuild import Compound


class Ellipsoid(Compound):
    def __init__(
            self, mass, length, name="dimer"
    ):
        super(Ellipsoid, self).__init__(name=name)
        self.length = float(major_length)
        n_particles = 2
        # Create the constituent particles
        self.head = Compound(
                pos=[self.length/2, 0, 0],
                name="CH",
                mass=mass/2
        )
        self.tail = Compound(
                pos=[-self.length/2, 0, 0],
                name="CT",
                mass=mass/2
        )
        self.add([self.head, self.tail])
