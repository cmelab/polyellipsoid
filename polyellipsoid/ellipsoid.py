from mbuild import Compound


class Ellipsoid(Compound):
    def __init__(self, mass, length):
        super(Ellipsoid, self).__init__(name="ellipsoid")
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
                pos=[self.length/2, 0, 0],
                name="B",
                mass=mass/4
        )
        self.tail_mid = Compound(
                pos=[-self.length/2, 0, 0],
                name="B",
                mass=mass/4
        )
        self.add([self.head, self.tail, self.head_mid, self.tail_mid])
