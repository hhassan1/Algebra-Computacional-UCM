from algebra_computacional.structures.EuclideanDomain import Euclid
from algebra_computacional.structures.Field import Field

class Fraction(Field):
    """docstring for Fraction"""
    def __init__(self, integral_domain, element):
        super(Fraction, self).__init__()
        self.integral = integral_domain
        if issubclass(integral_domain, Euclid):
            _,_,_,a,b = element[0].gcd(element[1])
            self.element = (a,b)
        else:
            self.element = element

    def is_zero(self):
        return self.element[0].is_zero()
    def is_one(self):
        return self.element[0].is_one() and self.element[1].is_one()
    def __eq__(self, op2):
        return self.element[0]*op2.element[1] == self.element[1]*op2.element[0]
    def __add__(self, op2):
        return Fraction(self.integral, (self.element[0]*op2.element[1] + self.element[1]*op2.element[0], self.element[1]*op2.element[1]))
    def __sub__(self, op2):
        return Fraction(self.integral, (self.element[0]*op2.element[1] - self.element[1]*op2.element[0], self.element[1]*op2.element[1]))
    def __neg__(self):
        return Fraction(self.integral, (-self.element[0],self.element[1]))
    def __mul__(self, op2):
        return Fraction(self.integral, (self.element[0]*op2.element[0],self.element[1]*op2.element[1]))
    def __div__(self, op2):
        return Fraction(self.integral, (self.element[0]*op2.element[1],self.element[1]*op2.element[0]))
    def __str__(self):
        if self.element[1] != self.integral.one():
            return str(self.element[0]) + '/' + str(self.element[1])
        return str(self.element[0])
