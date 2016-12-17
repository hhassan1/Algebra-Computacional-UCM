from algebra_computacional.structures.IntegralDomain import Integral
import itertools #sustituir por utilities haskell
class Polynomial(Integral):
    """docstring for Polynomial"""
    def __init__(self, field, coefficients):
        super(Polynomial, self).__init__()
        self.field = field
        self.coefficients = coefficients
    def __add__(self, rop):
        return Polynomial(self.field,{ k: self.coefs.get(k, self.field.zero()) + rop.coefs.get(k, self.field.zero()) for k in set(self.coefs) | set(rop.coefs) })
    def __sub__(self, rop):
        return Polynomial(self.field,{ k: self.coefs.get(k, self.field.zero()) + rop.coefs.get(k, self.field.zero()) for k in set(self.coefs) | set(rop.coefs) })
    def __mul__(self, rop):
        return Polynomial(self.field,{ k*l: self.coefs.get(k, self.field.one()) * rop.coefs.get(l, self.field.one()) for k,l in itertools.product(set(self.coefs),set(rop.coefs)) })
    def __neg__(self):
        return Polynomial(self.field,{ k: -v for k,v in self.coefs.iteritems()})
    def __eq__(self, rop):
        return not [0 for k in (set(self.coefs) | set(rop.coefs)) if not self.coefs.get(k, self.field.zero()) == rop.coefs.get(k, self.field.zero())]
    def __div__(self, rop):
        if rop.equalZero():
            raise ZeroDivisionError
        q, r = self.zero(), self
        while not r.equalZero() and r.degree() >= rop.degree:
            t = r.lead()/rop.lead()
            q, r = q+t, r - (t*rop)
        return q
    def is_zero(self):
        return self.coefs == {0:0L}
    def is_one(self):
        return self.coefs == {0:1L}
    def degree(self):
        return max(self.coefs, key=int)
    def lead(self):
        return self.coefs.get(self.degree)