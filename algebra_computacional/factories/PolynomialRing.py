from algebra_computacional.parsers.polynomial.polynomial import parse as polynomial_parse
from algebra_computacional.rings import Polynomial
import itertools
class PolynomialRing(object):
    def __init__(self, field):
        super(PolynomialRing, self).__init__()
        self.field = field
    def __call__(self, expression):
        if isinstance(expression, basestring):
            if expression:
                coefficients = polynomial_parse(self.field, expression)
            else:
                coefficients = {0:self.field.zero()}
        elif isinstance(expression, dict):
            coefficients = expression
        else:
            raise TypeError('Expression should be string or dictionary.')
        return Polynomial(self.field, coefficients)
    def zero(self):
        return Polynomial(self.field, {0:0L})
    def one(self):
        return Polynomial(self.field, {0:1L})
