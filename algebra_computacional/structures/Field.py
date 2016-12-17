from algebra_computacional.structures.EuclideanDomain import Euclid

class Field(Euclid):
    def __init__(self):
        super(Field, self).__init__()
    def __div__(self, op2):
        raise NotImplementedError
    def __invert__(self, op2):
        raise NotImplementedError
    