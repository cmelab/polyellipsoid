from mbuild import Compound


class Ellipsoid(Compound):
    def __init__(self, mass, length, name="dimer"):
        super(Ellipsoid, self).__init__(name=name)
        self.length = float(length)
        # Create the constituent particles
        self.head = Compound(
                pos=[self.length/2, 0, 0],
                name="A",
                mass=mass/4
        )
        self.tail = Compound(
                pos=[-self.length/2, 0, 0],
                name="A",
                mass=mass/4
        )
        self.head_mid = Compound(
                pos=self.head.xyz[0] / 2,
                name="B",
                mass=mass/4
        )
        self.tail_mid = Compound(
                pos=self.tail.xyz[0] / 2,
                name="B",
                mass=mass/4
        )
        self.add([self.head, self.tail, self.head_mid, self.tail_mid])
