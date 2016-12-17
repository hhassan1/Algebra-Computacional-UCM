from EuclideanDomain import Euclid
from Field import Field
class QuotientField(object):
    def __init__(self, integral_domain):
        super(QuotientField, self).__init__()
        self.integral = integral_domain
    def __call__(self, expression):
        class Fraction(Field):
            """docstring for Fraction"""
            def __init__(self, integral_domain, expression):
                super(Fraction, self).__init__()
                self.integral = integral_domain
                if isinstance(expression, basestring):
                    element = fracparse(integral_domain, expression)
                elif isinstance(expression, tuple):
                    element = expression
                else:
                    raise TypeError('Expression should be string or tuple.')
                if isinstance(integral_domain, Euclid):
                    _,_,_,a,b = self.integral.gcd(element[0],element[1])
                    self.element = (a,b)
                else:
                    self.element = element
            def zero(self):
                return (self.integral.zero(),self.integral.one())
            def one(self):
                return (self.integral.one(),self.integral.one())
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
        return Fraction(self.integral,expression)
