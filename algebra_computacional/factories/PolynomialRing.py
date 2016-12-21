from parsers.polynomial.polynomial import parse as polynomial_parse
from rings.Polynomial import Polynomial
import itertools
class PolynomialRing(object):
    def __init__(self, field):
        super(PolynomialRing, self).__init__()
        self.field = field
    def __call__(self, expression):
        if isinstance(expression, basestring):
            if expression:
                coefficients = polynomial_parse(self.field.builder(), expression)
            else:
                coefficients = {0L:self.field.zero()}
        elif isinstance(expression, dict):
            coefficients = expression
        else:
            raise TypeError('Expression should be string or dictionary.')
        return Polynomial(self.field, coefficients, self)
    def zero(self):
        return Polynomial(self.field, {0L:self.field.builder().zero()}, self)
    def one(self):
        return Polynomial(self.field, {0L:self.field.builder().one()}, self)
    def monomial(self, deg, scalar):
        return Polynomial(self.field, {deg:scalar}, self)
