from algebra_computacional.rings import (Integer,
                                         PolynomialOverIntegral,
                                         PolynomialOverField,
                                         PolynomialOverGalois,
                                         QuotientElement,
                                         ModularInteger,
                                         GaussianInteger,
                                         Fractional,
                                         FractionalPolynomial,
                                         PolynomialOverExtensionFractional
                                        )
import algebra_computacional.rings
from algebra_computacional.structures import (GaloisField,
                                              Field,
                                             )
from algebra_computacional.parsers import int_parse, poly_parse
from collections import OrderedDict
import random
import algebra_computacional.utilities.haskell as haskell


class IntegerFactory(object):
    """docstring for IntegerRing"""
    def __init__(self):
        super(IntegerFactory, self).__init__()
        self.product = Integer
        self.scalar_regex = r'\d+'
    def __call__(self, expression):
        if isinstance(expression, basestring):
            return int_parse(self, expression)
        elif isinstance(expression, (int, long)):
            return Integer(long(expression), self)
    def zero(self):
        return Integer(0L, self)
    def one(self):
        return Integer(1L, self)
    def scalar_translator(self, term):
        r'\d+'
        term.value = Integer(long(term.value), self)
        return term
    def random(self):
        return Integer(random.randint(-10000000L,10000000L),self)
class FractionalFactory(object):
    """docstring for IntegerRing"""
    def __init__(self):
        super(FractionalFactory, self).__init__()
        self.product = Fractional
        self.scalar_regex = r'\d+'
    def __call__(self, expression):
        if isinstance(expression, basestring):
            val = int_parse(self, expression)
            g = haskell.gcd(val.a,val.b)
            if g != 0:
                val.a /= g
                val.b /= g
            else:
                val.a = 0L
                val.b = 1L
            return val
        elif isinstance(expression, (int, long)):
            return Fractional(long(expression),1L, self)
        elif isinstance(expression, tuple):
            g = haskell.gcd(expression[0],expression[1])
            if g != 0:
                return Fractional(expression[0]/g,expression[1]/g, self)
            else:
                return Fractional(0L,1L, self)
    def zero(self):
        return Fractional(0L,1L, self)
    def one(self):
        return Fractional(1L,1L, self)
    def scalar_translator(self, term):
        r'\d+'
        term.value = Fractional(long(term.value),1L, self)
        return term
    def random(self):
        return self(random.randint(-10000000L,10000000L),random.randint(1,10000000L))


class PolynomialFactory(object):

    def __init__(self, factory, variable=r'x'):
        super(PolynomialFactory, self).__init__()
        self.inner_factory = factory
        self.variable = variable
        self.monic = True
        if factory.scalar_regex != r'\d+':
            self.scalar_regex = factory.scalar_regex + '|' + variable
        else:
            self.scalar_regex = variable
        if issubclass(factory.product, GaloisField):
            self.product = PolynomialOverGalois
        elif issubclass(factory.product, Field):
            self.product = PolynomialOverField
        else:
            self.product = PolynomialOverIntegral
    def __call__(self, expression):
        if isinstance(expression, basestring):
            if expression:
                polynomial = poly_parse(self, expression, self.variable)
            else:
                polynomial = self.zero()
            coefficients = polynomial.coefficients
            coefficients = OrderedDict(sorted({k: v for k, v in coefficients.iteritems() if not v.is_zero() or k == 0}.items(), reverse=True))
            if 0L in coefficients and coefficients[0L].is_zero() and len(coefficients) > 1:
                del coefficients[0L]
            if not coefficients:
                return self.zero()
            polynomial.coefficients = coefficients
            return polynomial
        elif isinstance(expression, self.inner_factory.product):
            return self.product({0L:expression}, self)
        elif isinstance(expression, dict):
            coefficients = OrderedDict(sorted({k: v for k, v in expression.iteritems() if not v.is_zero() or k == 0}.items(), reverse=True))
            if 0L in coefficients and coefficients[0L].is_zero() and len(coefficients) > 1:
                del coefficients[0L]
            if not coefficients:
                return self.zero()
            return self.product(coefficients, self)
        else:
            raise TypeError('Expression should be string or dictionary.')
    def zero(self):
        return self.product({0L:self.inner_factory.zero()}, self)
    def one(self):
        return self.product({0L:self.inner_factory.one()}, self)
    def monomial(self, degree, scalar=None):
        if not scalar:
            scalar = self.inner_factory.one()
            return self.product({degree:scalar}, self)
        if isinstance(scalar,self.inner_factory.product):
            return self.product({degree:scalar}, self)
        return self.product({degree:self.inner_factory(scalar)}, self)
    def random(self, lim=128):
        deg = random.randint(0,lim)
        return self.product({k:self.inner_factory.random() for k in range(deg+1)}, self)

class IntegerPolynomialFactory(PolynomialFactory):
    def __init__(self):
        super(IntegerPolynomialFactory, self).__init__(IntegerFactory())
        self.product = algebra_computacional.rings.IntegerPolynomial




class QuotientFactory(object):
    def __init__(self, factory, expression):
        super(QuotientFactory, self).__init__()
        self.inner_factory = factory
        self.scalar_regex = self.inner_factory.scalar_regex
        self.product = QuotientElement
        if isinstance(expression,factory.product):
            self.quotient_factor = expression
        else:
            self.quotient_factor = factory(expression)
    def __call__(self, expression):
        if isinstance(expression, (basestring,int,float)):
            ret = self.inner_factory(expression)
        else:
            ret =  expression
        return self.product(ret % self.quotient_factor, self)
    def zero(self):
        return self.product(self.inner_factory.zero(), self)
    def one(self):
        return self.product(self.inner_factory.one(), self)
    def monomial(self, degree, scalar=None):
        return self.product(self.inner_factory.monomial(degree, scalar) % self.quotient_factor, self)
    def random(self):
        return self.product(self.inner_factory.random() % self.quotient_factor, self)

class ModularIntegerFactory(QuotientFactory):
    def __init__(self, p):
        super(ModularIntegerFactory, self).__init__(IntegerFactory(), p)
        self.p = p
        self.k = 1
        self.q = p
        self.product = ModularInteger
    def random(self):
        num = random.randint(0,self.p)
        return self(num)
    def element_iterator(self):
        return GaloisIterator(self)
    def monomial(self, degree, scalar=None):
        if scalar:
            return self(scalar)
        return self('1')

class GaloisQuotientFactory(QuotientFactory):
    def __init__(self, p, expression):
        import re
        variable = re.search('([a-zA-Z])', expression).group(1)
        if not variable:
            variable = 'x'
        super(GaloisQuotientFactory, self).__init__(PolynomialFactory(ModularIntegerFactory(p), variable), expression)
        self.p = p
        self.k = self.quotient_factor.degree()
        self.q = p ** self.k
    def random(self):
        return self.product(self.inner_factory.random(self.k-1), self)
    def element_iterator(self):
        return GaloisIterator(self)

class GaloisPolynomialFactory(PolynomialFactory):
    def __init__(self, p, expression=None, variable='x',lifted = False):
        if expression:
            factory = GaloisQuotientFactory(p,expression)
            self.k = factory.k
            self.lifted = False
        else:
            factory = ModularIntegerFactory(p)
            self.k = 1
            self.lifted = lifted
        super(GaloisPolynomialFactory, self).__init__(factory,variable)
        self.quotient_factor = self.inner_factory.quotient_factor
        self.p = p
        self.q = p ** self.k
        self.product = PolynomialOverGalois
    def lift(self, n):
        if n <= 1:
            return self
        s = n*self.k
        if self.k > 1:
            f = GaloisPolynomialFactory(self.p)
        else:
            f = self
        item = f.monomial(s)
        iterators = [f.inner_factory.element_iterator() for i in range(s)]
        while len(item.factors2) > 1:
            i = 0
            sum_bool = False
            while sum_bool and i < s:
                item.coefficient[i] = iterators[i].next()
                item = item + f.monomial(i)
                if not iterators[i].has_next():
                    iterators[i].reset()
                else:
                    sum_bool = False
                i += 1
        return GaloisPolynomialFactory(self.p, str(item))


                

class GaussianIntegerFactory(QuotientFactory):
    def __init__(self):
        super(GaussianIntegerFactory, self).__init__(PolynomialFactory(IntegerFactory(),'i'),'i^2 + 1')
        self.product = GaussianInteger

class GaloisIterator(object):
    """docstring for GaloisIterator"""
    def __init__(self, factory):
        super(GaloisIterator, self).__init__()
        self.p = factory.p
        self.k = factory.k
        self.counter = 0
        self.limit = factory.q
        self.item = factory.zero()
        self.factory = factory
        self.counters = [0 for i in range(self.k)]
    def next(self):
        item = self.item
        sum_bool = True
        i = 0
        while sum_bool and i < self.k:
            self.counters[i] += 1
            self.item = self.item + self.factory.monomial(i)
            if self.counters[i] == self.p:
                self.counters[i] = 0
            else:
                sum_bool = False
            i += 1
        self.counter+=1
        return item
    def has_next(self):
        return self.counter < self.limit
    def reset(self):
        self.counter = 0
        self.counters = [0 for i in range(self.k)]
        self.item = self.factory.zero()

class FractionalPolynomialFactory(PolynomialFactory):
    def __init__(self, variable='x'):
        super(FractionalPolynomialFactory, self).__init__(FractionalFactory(),variable)
        self.product = FractionalPolynomial
        #self.monic = False

class FractionalQuotientFactory(QuotientFactory):
    def __init__(self, expression, factory=None):
        import re
        variable = re.search('([a-zA-Z])', expression).group(1)
        if not variable:
            variable = 'x'
        self.variable = variable

        if not factory:
            super(FractionalQuotientFactory, self).__init__(FractionalPolynomialFactory(variable), expression)
            self.delegate = None
            self.degree = self.quotient_factor.degree()
            self.dimension_chart = (self.degree,)
        else:
            super(FractionalQuotientFactory, self).__init__(PolynomialFactory(factory,variable), expression)
            self.delegate = factory
            self.degree = self.quotient_factor.degree()*factory.degree
            self.dimension_chart = (self.quotient_factor.degree(),) + factory.dimension_chart
    def separate_basis(self, poly):
        QX = FractionalPolynomialFactory()
        table = dict()
        for k in range(self.dimension_chart[0]):
            if not self.delegate:
                table[k] = QX(poly.value.coefficient(k))
            else:
                for j,u in self.delegate.separate_basis(poly.value.coefficient(k)).iteritems():
                    table[self.delegate.degree*k + j] = u
        return table
    def minimal_polynomial(self, element):
        QX = FractionalPolynomialFactory()
        basis = self.basis()
        L = [element*y for y in basis]
        j = 0
        s = len(L)
        M = dict()
        for i in L:
            T = self.separate_basis(i)
            T[j] = QX('-x')
            for k in range(s):
                M[j,k] = T[k]
            j+=1
        facs = determinant(M,s,s*s).factors()
        QabX = PolynomialFactory(self)
        for x in facs:
            result = QabX(str(x))
            if result(str(element)).is_zero():
                return x


    def basis(self):
        if isinstance(self.inner_factory,FractionalPolynomialFactory):
            for i in range(0,self.quotient_factor.degree()):
                yield self(self.variable+'^'+str(i))
        else:
            for i in range(0,self.quotient_factor.degree()):
                for j in self.inner_factory.inner_factory.basis():
                    yield self('(' + self.variable+'^'+str(i) + ')*('+ str(j) +')')

def determinant(m,n,s):
    if s == 1:
        return m[0,0]
    result = m[0,0].factory.zero()
    for i in range(n):
        minor = {(j-1,(k if k < i else k-1)):v for (j,k),v in m.iteritems() if k != i}
        if i % 2 == 0:
            result += m[0,i]*determinant(minor,n-1,s-2*n+1)
        else:
            result -= m[0,i]*determinant(minor,n-1,s-2*n+1)
    return result
    



class ExtensionFractionalPolynomialFactory(PolynomialFactory):
    def __init__(self, expression=None, variable='x',factory=None):
        if factory:
            super(ExtensionFractionalPolynomialFactory, self).__init__(factory,variable)
        else:
            if expression:
                factory = FractionalQuotientFactory(expression)
                self.quotient_factor = factory.quotient_factor
            else:
                factory = FractionalPolynomialFactory(variable)
            super(ExtensionFractionalPolynomialFactory, self).__init__(factory,variable)
            self.product = PolynomialOverExtensionFractional
    def quotient(self, q, var):
        return FractionalQuotientFactory(str(q),var,self)