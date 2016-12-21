from structures.IntegralDomain import Integral
from structures.EuclideanDomain import Euclid
import itertools #sustituir por utilities haskell
class Polynomial(Euclid):
    """docstring for Polynomial"""
    def __init__(self, field, coefficients, builder):
        super(Polynomial, self).__init__()
        self.field = field
        self.coefficients = {k: v for k, v in coefficients.iteritems() if not v.is_zero() or k == 0}
        self.factory = builder
    def __add__(self, rop):
        return Polynomial(self.field,{ k: self.coefficients.get(k, self.field.builder().zero()) + rop.coefficients.get(k, self.field.builder().zero()) for k in set(self.coefficients) | set(rop.coefficients) }, self.factory)
    def __sub__(self, rop):
        return Polynomial(self.field,{ k: self.coefficients.get(k, self.field.builder().zero()) - rop.coefficients.get(k, self.field.builder().zero()) for k in set(self.coefficients) | set(rop.coefficients) }, self.factory)
    def __mul__(self, rop):
        a = self.factory.zero()
        for k in self.coefficients.iterkeys():
            a = a + self.factory({ j+k: self.coefficients[k]*rop.coefficients[j] for j in rop.coefficients.iterkeys() })
        return a
    def __neg__(self):
        return Polynomial(self.field,{ k: -v for k,v in self.coefficients.iteritems()},self.factory)
    def __eq__(self, rop):
        return not [0 for k in (set(self.coefficients) | set(rop.coefficients)) if not self.coefficients.get(k, self.field.zero()) == rop.coefficients.get(k, self.field.zero())]
    def divmod(self, rop):
        if rop.is_zero():
            raise ZeroDivisionError
        q, r = self.factory.zero(), self
        while not r.is_zero() and r.degree() >= rop.degree():
            print q, r
            t = self.factory.monomial(r.degree() - rop.degree(),r.leading_coefficient() / rop.leading_coefficient())
            q = q+t
            r = r - (t*rop)
            raw_input('')
        return q, r
    def __div__(self, rop):
        return self.divmod(rop)[0]
    def __mod__(self, rop):
        return self.divmod(rop)[1]
    def is_zero(self):
        return self.coefficients == {0L:self.field.builder().zero()}
    def is_one(self):
        return self.coefficients == {0L:self.field.builder().one()}
    def degree(self):
        return max(self.coefficients, key=int)
    def leading_coefficient(self):
        return self.coefficients[self.degree()]
    def builder(self):
        return self.factory
    def __repr__(self):
        return ' + '.join(reversed([ ( '(' + repr(item[1]) + ')' if not item[1].is_one() else '') + ('x' if item[0] > 0 else '') + ('^' + repr(item[0])[:-1]  if item[0] > 1 else '') for item in self.coefficients.iteritems()]))