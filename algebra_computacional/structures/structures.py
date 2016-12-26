__all__ = ['Integral', 'EuclideanDomain', 'Field', 'GaloisField']

from functools import reduce
from algebra_computacional.utilities.haskell import scanl, scanr, zipWith


class Integral(object):
    """docstring for Integral"""
    def __init__(self):
        super(Integral, self).__init__()
        self.factory = None
    def builder(self):
        return self.factory
    def zero(self):
        return self.factory.zero()
    def one(self):
        return self.factory.one()
    def is_zero(self):
        return self == self.factory.zero()
    def is_one(self):
        return self == self.factory.one()
    def __eq__(self, op2):
        return not self != op2
    def __ne__(self, op2):
        return not self == op2
    def __add__(self, op2):
        raise NotImplementedError
    def __sub__(self, op2):
        raise NotImplementedError
    def __neg__(self):
        raise NotImplementedError
    def __mul__(self, op2):
        raise NotImplementedError
    def __xor__(self, exp):
        if exp == 0:
            return self.one()
        elif exp == 1:
            return self
        ret = self ^ (exp/2)
        ret *= ret
        if (exp % 2) == 1:
            ret *= self
        return ret
    def __div__(self, op2):
        raise NotImplementedError
    def __mod__(self, op2):
        raise NotImplementedError

class EuclideanDomain(Integral):
    """docstring for Euclid"""
    def __init__(self):
        super(EuclideanDomain, self).__init__()
    def eea(self, b):
        R = [self, b]
        S = [self.one(), self.zero()]
        T = [self.zero(), self.one()]
        i = 1
        while not R[i].is_zero():
            Q = R[i-1] / R[i]
            R.append(R[i-1] - Q*R[i])
            S.append(S[i-1] - Q*S[i])
            T.append(T[i-1] - Q*T[i])
            i = i + 1
        g = R[i-1]
        u, v = S[i-1], T[i-1]
        if (-g).is_one():
            return -g, -u, -v
        return g, u, v
    @classmethod
    def chinese_remainder(cls, equations):
        a_list, p_list = zip(*equations)
        forwards_products = scanl((lambda x,y: x*y), p_list[0], p_list[1:-1])
        backwards_products = list(reversed(scanr((lambda x,y: x*y), p_list[-1], p_list[1:-1])))
        N = [backwards_products[0]] + zipWith((lambda x,y: x*y), forwards_products[:-1], backwards_products[1:]) + [forwards_products[-1]]
        M = zipWith((lambda a,b: cls.eea(a,b)[1]), N, p_list)
        sums = zipWith((lambda x,y: x*y), a_list, zipWith((lambda x,y: x*y), M, N))
        return (reduce((lambda x,y: x+y), sums[1:], sums[0])) % (forwards_products[-1]*p_list[-1])

class Field(EuclideanDomain):
    def __init__(self):
        super(Field, self).__init__()

class GaloisField(Field):
    def __init__(self):
        super(GaloisField, self).__init__()