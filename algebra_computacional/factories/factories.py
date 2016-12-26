from algebra_computacional.rings import (Integer,
                                         PolynomialOverIntegral,
                                         PolynomialOverField,
                                         PolynomialOverGalois,
                                         QuotientElement,
                                         ModularInteger,
                                         GaussianInteger,
                                        )
from algebra_computacional.structures import (GaloisField,
                                              Field,
                                             )
from algebra_computacional.parsers import int_parse, poly_parse
from collections import OrderedDict


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

class QuotientFactory(object):
    def __init__(self, factory, expression):
        super(QuotientFactory, self).__init__()
        self.inner_factory = factory
        self.scalar_regex = self.inner_factory.scalar_regex
        self.product = QuotientElement
        self.quotient_factor = factory(expression)
    def __call__(self, expression):
        if isinstance(expression, basestring):
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

class ModularIntegerFactory(QuotientFactory):
    def __init__(self, p):
        super(ModularIntegerFactory, self).__init__(IntegerFactory(), p)
        self.p = p
        self.product = ModularInteger
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

class GaussianIntegerFactory(QuotientFactory):
    def __init__(self):
        super(GaussianIntegerFactory, self).__init__(PolynomialFactory(IntegerFactory(),'i'),'i^2 + 1')
        self.product = GaussianInteger