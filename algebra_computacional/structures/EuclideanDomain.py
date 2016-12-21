from IntegralDomain import Integral
from utilities.haskell import scanl, scanr, zipWith
from functools import reduce

class Euclid(Integral):
    """docstring for Euclid"""
    def __init__(self):
        Integral.__init__(self)
    def eea(self, b):
        Q = [self.zero(), self.zero()]
        R = [self, b]
        S = [self.one(), self.zero()]
        T = [self.zero(), self.one()]
        i = 1
        while not R[i].is_zero():
            print R[i]
            Q.append(R[i-1] / R[i])
            R.append(R[i-1] - Q[i+1]*R[i])
            S.append(S[i-1] - Q[i+1]*S[i])
            T.append(T[i-1] - Q[i+1]*T[i])
            print i
            i = i + 1
        g = R[i-1]
        u, v = S[i-1], T[i-1]
        return g, u, v
    @classmethod
    def chinese_remainder(cls, equations):
        a_list, p_list = zip(*equations)
        forwards_products = scanl(cls.__mul__, 1, p_list[:-1])
        backwards_products = scanr(cls.__mul__, 1, p_list[1:])
        N = zipWith(cls.__mul__, forwards_products, reversed(backwards_products))
        M = zipWith((lambda a,b: cls.eea(a,b)[1]), N, p_list)
        sums = zipWith(cls.__mul__, a_list, zipWith(cls.__mul__, M, N))
        return cls.__mod__(reduce(cls.__add__, sums, 0),forwards_products[-1]*p_list[-1])
    def __div__(self, op2):
        raise NotImplementedError
    def __mod__(self, op2):
        raise NotImplementedError