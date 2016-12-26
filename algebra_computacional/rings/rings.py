__all__ = ['Integer', 'PolynomialOverIntegral', 'PolynomialOverField',
           'PolynomialOverGalois', 'QuotientElement', 'ModularInteger',
           'GaussianInteger']

from algebra_computacional.structures import (EuclideanDomain,
                                              Field,
                                              GaloisField,
                                              Integral
                                             )
from algebra_computacional.utilities.haskell import zipWithAll, concat, eea
import itertools #sustituir por utilities haskell

class Integer(EuclideanDomain):
    """docstring for Int"""
    def __init__(self, value, factory):
        super(Integer, self).__init__()
        self.number = value
        self.factory = factory
    def is_one(self):
        return self.number == 1
    def is_zero(self):
        return self.number == 0
    def __mod__(self, op2):
        return self.factory(self.number % op2.number)
    def __eq__(self, op2):
        return self.number == op2.number
    def __add__(self, op2):
        return self.factory(self.number + op2.number)
    def __sub__(self, op2):
        return self.factory(self.number - op2.number)
    def __neg__(self):
        return self.factory(-self.number)
    def __mul__(self, op2):
        return self.factory(self.number * op2.number)
    def __div__(self, op2):
        return self.factory(self.number / op2.number)
    def __str__(self):
        return str(self.number)
    def __repr__(self):
        return str(self.number)

class QuotientElement(Field):
    """docstring for Int"""
    def __init__(self, value, factory):
        super(QuotientElement, self).__init__()
        self.value = value
        self.factory = factory
    def is_one(self):
        return (self.value - self.factory.one().value).is_zero()
    def is_zero(self):
        return (self.value - self.factory.zero().value).is_zero()
    def __mod__(self, op2):
        return self.factory(self.value % op2.value)
    def __eq__(self, op2):
        return (self.value - op2.value).is_zero()
    def __add__(self, op2):
        return self.factory(self.value + op2.value)
    def __sub__(self, op2):
        return self.factory(self.value - op2.value)
    def __neg__(self):
        return self.factory(-self.value)
    def __mul__(self, op2):
        return self.factory(self.value * op2.value)
    def inverse(self):
        g, u, v = self.value.eea(self.factory.quotient_factor)
        if not g.is_one():
            u = -u
        return self.factory(u)
    def __divmod__(self, op2):
        return self.factory(divmod(self.value,op2.value))
    def __div__(self, op2):
        return self.factory(self.value / op2.value)
    def __str__(self):
        return str(self.value)
    def __repr__(self):
        return repr(self.value)

class GaussianInteger(QuotientElement):
    """docstring for Int"""
    def __init__(self, value, factory):
        super(GaussianInteger, self).__init__(value, factory)
    def __div__(self, op2):
        if op2.is_zero():
            raise ZeroDivisionError
        conj = op2.conjugate()
        a = self * conj
        b = op2 * conj
        nat_a = a.value.coefficients.get(0,self.factory.inner_factory.inner_factory.zero()).number
        img_a = a.value.coefficients.get(1,self.factory.inner_factory.inner_factory.zero()).number
        nat_b = b.value.coefficients[0L].number
        nat_res = int(round(nat_a*(1.0/nat_b)))
        img_res = int(round(img_a*(1.0/nat_b)))
        return self.factory( '(' + str(nat_res) + ') + (' + str(img_res) + ')i' )
    def conjugate(self):
        copy = self.value + self.factory.inner_factory.zero()
        if 1L in self.value.coefficients:
            copy.coefficients[1L] = -self.value.coefficients[1L]
        return self.factory(copy)

class ModularInteger(QuotientElement, GaloisField):
    """docstring for Int"""
    def __init__(self, value, factory):
        super(ModularInteger, self).__init__(value, factory)
    def inverse(self):
        g, u, v = eea(self.value.number,self.factory.quotient_factor.number)
        if g == -1:
            u = -u
        return self.factory(self.factory.inner_factory(u))
    def __div__(self, op2):
        return self.factory((self * op2.inverse()).value)

class PolynomialOverIntegral(Integral):
    """docstring for Polynomial"""
    def __init__(self, coefficients, builder):
        super(PolynomialOverIntegral, self).__init__()
        self.coefficients = coefficients
        self.factory = builder
    def __add__(self, rop):
        return self.factory({k: self.coefficients.get(k, self.factory.inner_factory.zero()) + rop.coefficients.get(k, self.factory.inner_factory.zero()) for k in set(self.coefficients) | set(rop.coefficients)})
    def __sub__(self, rop):
        return self.factory({k: self.coefficients.get(k, self.factory.inner_factory.zero()) - rop.coefficients.get(k, self.factory.inner_factory.zero()) for k in set(self.coefficients) | set(rop.coefficients)})
    def __mul__(self, rop):
        a = self.factory.zero()
        for k in self.coefficients.iterkeys():
            a = a + self.factory({ j+k: self.coefficients[k]*rop.coefficients[j] for j in rop.coefficients.iterkeys() })
        return a
    def __neg__(self):
        return self.factory({ k: -v for k,v in self.coefficients.iteritems()})
    def __eq__(self, rop):
        return (self - rop).is_zero()
    def is_zero(self):
        return self.coefficients == {0L:self.factory.inner_factory.zero()}
    def is_one(self):
        return self.coefficients == {0L:self.factory.inner_factory.one()}
    def degree(self):
        return max(self.coefficients, key=int)
    def leading_coefficient(self):
        return self.coefficients[self.degree()]
    def builder(self):
        return self.factory
    def derivative(self):
        return self.factory({ k-1: self.factory.inner_factory(str(k)) * v for k,v in self.coefficients.iteritems() if k != 0})
    def __str__(self):
        return ' + '.join([(str(item[1]) if not item[1].is_one() or item[0] == 0 else '') +\
                     (self.factory.variable if item[0] > 0 else '') +\
                      ('^' + str(item[0])  if item[0] > 1 else '') \
                      for item in self.coefficients.iteritems()])
    def __call__(self, expression):
        subs = self.factory.inner_factory(expression)
        result = self.factory.inner_factory.zero()
        for k in xrange(self.degree(),-1,-1):
            scalar = self.coefficients.get(k,self.factory.inner_factory.zero())
            result = result*subs + scalar
    def __divmod__(self, rop):
        if rop.is_zero():
            raise ZeroDivisionError
        q = self.factory.zero()
        if (self.leading_coefficient() / rop.leading_coefficient())*rop.leading_coefficient() != self.leading_coefficient():
            r = self.factory.monomial(0, self.leading_coefficient()) * self
        else:
            r = self
        while not r.is_zero() and r.degree() >= rop.degree():
            t = self.factory.monomial(r.degree() - rop.degree(), r.leading_coefficient() / rop.leading_coefficient())
            q = q+t
            r = r - (t*rop)
        return q, r
    def __div__(self, rop):
        return divmod(self, rop)[0]
    def __mod__(self, rop):
        return divmod(self, rop)[1]



class PolynomialOverField(PolynomialOverIntegral, EuclideanDomain):
    """docstring for Polynomial"""
    def __init__(self, coefficients, builder):
        super(PolynomialOverField, self).__init__(coefficients, builder)
    def __divmod__(self, rop):
        if rop.is_zero():
            raise ZeroDivisionError
        q = self.factory.zero()
        r = self
        while not r.is_zero() and r.degree() >= rop.degree():
            t = self.factory.monomial(r.degree() - rop.degree(), r.leading_coefficient() / rop.leading_coefficient())
            q = q+t
            r = r - (t*rop)
        return q, r
    def __div__(self, rop):
        return divmod(self, rop)[0]
    def __mod__(self, rop):
        return divmod(self, rop)[1]


def list_power(l, exponent, one):
    return concat([[one() for _ in range(1,exponent)] + [x] for x in l])

class PolynomialOverGalois(PolynomialOverField, Field):
    def __init__(self, coefficients, builder):
        super(PolynomialOverGalois, self).__init__(coefficients, builder)
    def inverse(self):
        return self.factory.one() / self
    def __divmod__(self, rop):
        if rop.is_zero():
            raise ZeroDivisionError
        q = self.factory.zero()
        r = self
        while not r.is_zero() and r.degree() >= rop.degree():
            t = self.factory.monomial(r.degree() - rop.degree(), r.leading_coefficient() * rop.leading_coefficient().inverse())
            q = q+t
            r = r - (t*rop)
        return q, r
    def gcd(self, b):
        return self.eea(b)[0]

    def pth_root(self):
        return self.factory({ k/self.factory.p: v for k,v in self.coefficients.iteritems()})
    
    def frob_basis(self):
        base_size = self.degree()
        image = [[None]]*self.degree()
        for canonic in [self.factory.monomial(i) for i in range(base_size)]:
            image.append((canonic^(self.factory.q) - canonic)%self)
        result = gaussian_elimination(image)
        return [cls.one(), result]

    def berlekamp(self):
        basis = self.frob_basis()
        factors = [self]
        irreducibles = []
        while len(factors) + len(irreducibles) < len(basis):
            g = basis.pop()
            h = g.berlekamp_splitting(prebase)
            if h != g:
                factors.append(h)
                factors.append(g/h)
            else:
                irreducibles.append(g)
        return factors + irreducibles

    def squarefree(self):
        derivative = self.derivative()
        if derivative.is_zero():
            return list_power(self.pth_root().squarefree(),self.factory.p, self.factory.one)
        L = []
        g = [self.gcd(derivative)]
        w = [self / g[0]]
        while not w[-1].is_one():
            w.append(g[-1].gcd(w[-1]))
            g.append(g[-1]/w[-1])
            L.append(w[-2]/w[-1])
        if g[-1].is_one():
            return L
        return zipWithAll((lambda x,y: x * y), (lambda x: x), L, list_power(g[-1].pth_root().squarefree(), self.factory.p, self.factory.one))

    def distinctdegree(self):
        d = 0
        Result = []
        g = [self]
        x = self.factory.monomial(1)
        h = [x]
        while d <= (g[d].degree() / 2 - 1):
            d = d + 1
            h.append((h[d-1]^(self.factory.q)) % g[d-1])
            aux = g[d-1].gcd(h[d] - x)
            if not aux.is_one():
                Result.append((aux,d))
            g.append(g[d-1]/aux)
        if not g[d].is_one():
            Result.append((g[d],g[d].degree()))
        return Result

    def berlekamp_splitting(self, prebase):
        s = len(prebase)
        if s is 1:
            return self
        i = 2
        elements = self.factory.element_iterator()
        while i <= s:
            while elements.has_next():
                element = elements.next()
                g = self.gcd(prebase[i] - element)
                if not (g.is_one() or g == self):
                    return g
            i = i + 1
            elements.restart()
