from structures.EuclideanDomain import Euclid
import factories.IntegerRing

class Integer(Euclid):
    """docstring for Int"""
    def __init__(self, value):
        super(Integer, self).__init__()
        self.number = long(value)
    def is_one(self):
        return self.number == 1
    def is_zero(self):
        return self.number == 0
    def __eq__(self, op2):
        return self.number == op2.number
    def __add__(self, op2):
        return Integer(self.number + op2.number)
    def __sub__(self, op2):
        return Integer(self.number - op2.number)
    def __neg__(self):
        return Integer(-self.number)
    def __mul__(self, op2):
        return Integer(self.number * op2.number)
    def __str__(self):
        return str(self.number)
    def __repr__(self):
        return repr(self.number)[:-1]
    def __div__(self, op2):
        return Integer(self.number // op2.number)
    @staticmethod
    def builder():
        return factories.IntegerRing.IntegerRing
