from algebra_computacional.rings.Fraction import Fraction
from algebra_computacional.parsers.fraction.fraction import FractionParser

class QuotientField(object):
    def __init__(self, integral_domain):
        super(QuotientField, self).__init__()
        self.integral = integral_domain
        self.fraction_parse = FractionParser(self, integral_domain)
    def __call__(self, expression):
        if isinstance(expression, basestring):
            element = self.fraction_parse(expression)
        elif isinstance(expression, tuple):
            if len(expression) == 2:
                element = expression
            elif len(expression) == 1:
                element = (expression[0], self.integral.one())
            else:
                raise IndexError('Expression have at least one element.')
        else:
            raise TypeError('Expression should be string or tuple.')
        return Fraction(self.integral, element)
    def zero(self):
        return (self.integral.zero(),self.integral.one())
    def one(self):
        return (self.integral.one(),self.integral.one())
