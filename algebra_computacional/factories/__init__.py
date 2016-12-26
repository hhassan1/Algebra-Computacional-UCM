#from .FractionFactory import *
from .factories import *

Z_X = PolynomialFactory(IntegerFactory())

def GF(p, expression=None):
    if expression is not None:
        return GaloisQuotientFactory(p, expression)
    else:
        return ModularIntegerFactory(p)