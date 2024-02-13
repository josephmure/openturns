import openturns as ot

class Stupid(ot.DistributionImplementation):

    def __init__(self):
        super().__init__()
    
    def getRealization(self):
        return ot.Point([0.123])
    
s = Stupid()
print(s.getSample(3))

dist = ot.Distribution(s)
print(dist.getSample(3))

rv = ot.UsualRandomVector(s)
print(rv.getRealization())