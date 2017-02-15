from algebra_computacional.rings import (Integer,
                                         PolynomialOverIntegral,
                                         PolynomialOverField,
                                         PolynomialOverGalois,
                                         QuotientElement,
                                         ModularInteger,
                                         GaussianInteger,
                                         
                                        )
import algebra_computacional.rings
from algebra_computacional.structures import (GaloisField,
                                              Field,
                                             )
from algebra_computacional.parsers import int_parse, poly_parse
from collections import OrderedDict
import random


class IntegerFactory(object):
    """docstring for IntegerRing"""
    def __init__(self):
        super(IntegerFactory, self).__init__()
        self.product = Integer
        self.scalar_regex = r'\d+'
    def __call__(self, expression):
        if isinstance(expression, basestring):
            return Integer(int_parse(expression), self)
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
        return random.randint(-10000000L,10000000L)


class PolynomialFactory(object):

    def __init__(self, factory, variable=r'x'):
        super(PolynomialFactory, self).__init__()
        self.inner_factory = factory
        self.variable = variable
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
    def __init__(self, p, expression=None, variable='x'):
        if expression:
            factory = GaloisQuotientFactory(p,expression)
            self.k = factory.k
        else:
            factory = ModularIntegerFactory(p)
            self.k = 1
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