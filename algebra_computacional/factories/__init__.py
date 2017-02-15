#from .FractionFactory import *
from .factories import *

Z_X = IntegerPolynomialFactory()

def GF_X(p, expression=None, variable='x'):
    return GaloisPolynomialFactory(p, expression, variable)
def GF(p, expression=None):
    if expression is not None:
        return GaloisQuotientFactory(p, expression)
    else:
        return ModularIntegerFactory(p)