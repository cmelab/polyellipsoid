from mbuild import Compound


class Ellipsoid(Compound):
    def __init__(
            self, mass, length, name="dimer"
    ):
        super(Ellipsoid, self).__init__(name=name)
        self.length = float(length)
        # Create the constituent particles
        self.head = Compound(
                pos=[0, 0, self.length/2],
                name="CH",
                mass=mass/2
        )
        self.tail = Compound(
                pos=[0, 0, -self.length/2],
                name="CT",
                mass=mass/2
        )
        self.mid_head = Compound(
                pos=[0,0,self.length/4]
                name="CHM"
                mass=0
        )
        self.mid_tail = Compound(
                pos=[0,0,-self.length/4]
                name="CTM"
                mass=0
        )
        self.add([self.head, self.tail])
        self.add_bond([self.head, self.mid_head])
        self.add_bond([self.tail, self.mid_tail])
