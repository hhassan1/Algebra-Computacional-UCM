__all__ = ['Integer', 'PolynomialOverIntegral', 'PolynomialOverField',
           'PolynomialOverGalois', 'QuotientElement', 'ModularInteger',
           'GaussianInteger', 'IntegerPolynomial','Fractional','PolynomialOverExtensionFractional','FractionalPolynomial']

from algebra_computacional.structures import (EuclideanDomain,
                                              Field,
                                              GaloisField,
                                              Integral
                                             )
from algebra_computacional.utilities.haskell import zipWithAll, concat, eea, foldl, gcd as normal_gcd, lcm as normal_lcm
import algebra_computacional.factories
import math
import fractions
import random
import itertools

class Integer(EuclideanDomain):
    """docstring for Int"""
    def __init__(self, value, factory):
        super(Integer, self).__init__()
        self.number = value
        self.factory = factory
    def is_perfect_power(self):
        b = 2
        c = math.log(self.number,2)
        while b <= c:
            a = self.number**(1./b)
            b+=1
            if math.floor(a) == a:
                return True
        return False
    def phi(self):
        amount = 0

        for k in range(1, self.number + 1):
            if fractions.gcd(self.number, k) == 1:
                amount += 1

        return amount
    def is_prime(self):
        if self.is_perfect_power():
            return False
        maxk=math.floor(math.log(self.number)**2)
        nextR=True
        r=2
        while nextR:
            nextR=False
            k = 1
            while (not nextR) and k <= maxk:
                aux = (self.number**k % r)
                nextR= aux == 1 or aux == 0
                k+=1
            r+=1
        r-=1
        for a in xrange(2,min(r,self.number-1)+1):
            if self.number % a == 0:
                return False
        if self.number <= r:
            return True
        GF_n = algebra_computacional.factories.GaloisPolynomialFactory(self.number)
        X = GF_n.monomial(1)
        X_n = GF_n.monomial(self.number)
        mod = GF_n.monomial(r) - GF_n.one()
        for a in xrange(1,  int(self.phi()*math.log(self.number,2))+1):
            A = GF_n.monomial(0,a)
            M = (((X + A)^self.number) - X_n - A) % mod
            if not M.is_zero():
                return False
        return True
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
class Fractional(Field):
    """docstring for Fractional"""
    def __init__(self, a, b, factory):
        super(Fractional, self).__init__()
        self.a = a
        self.b = b
        self.factory = factory
    def is_one(self):
        return self.a == 1 and self.b ==1
    def is_zero(self):
        return self.a == 0
    def __mod__(self, op2):
        return self.factory.zero()
    def __eq__(self, op2):
        return self.a*op2.b == op2.a*self.b
    def __add__(self, op2):
        return self.factory((self.a*op2.b + op2.a*self.b,self.b*op2.b))
    def __sub__(self, op2):
        return self.factory((self.a*op2.b - op2.a*self.b,self.b*op2.b))
    def __neg__(self):
        return self.factory((-self.a,self.b))
    def __mul__(self, op2):
        return self.factory((self.a*op2.a,self.b*op2.b))
    def __div__(self, op2):
        return self.factory((self.a*op2.b,self.b*op2.a))
    def __str__(self):
        if self.b == 1:
            return str(self.a)
        return str(self.a) + '/' + str(self.b)
    def __repr__(self):
        if self.b == 1:
            return str(self.a)
        return str(self.a) + '/' + str(self.b)




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
        return self * op2.inverse()
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
        if value.number > (factory.p)/2:
            value.number -= factory.p
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
        return self.factory({k: self.coefficient(k) + rop.coefficient(k) for k in set(self.coefficients) | set(rop.coefficients)})
    def __sub__(self, rop):
        return self.factory({k: self.coefficient(k) - rop.coefficient(k) for k in set(self.coefficients) | set(rop.coefficients)})
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
    def coefficient(self, k):
        return self.coefficients.get(k,self.factory.inner_factory.zero())
    def flip(self, n):
        return self.factory({ n-k: v for k,v in self.coefficients.iteritems() })
    def pp(self):
        if(self.is_zero()):
            return self
        gcd = foldl(self.factory.inner_factory.zero().gcd.__func__, self.coefficients.values()[0], self.coefficients.itervalues())
        return self.factory({k: v/gcd for k, v in self.coefficients.iteritems()})
    def pea(self, b):
        R = [self, b]
        i = 1
        while not R[i].is_zero():
            A = ((self.factory(str(R[i].leading_coefficient()))^(R[i-1].degree() - R[i].degree() + 1))*R[i-1] % R[i])
            R.append(A.pp())
            i = i + 1
        g = R[i-1]
        if (-g.leading_coefficient()).is_one():
            return -g
        return g
    def derivative(self):
        return self.factory({ k-1: self.factory.inner_factory(str(k)) * v for k,v in self.coefficients.iteritems() if k != 0})
    def __str__(self):
        aux = ' + '.join([(str(item[1]) if not item[1].is_one() or item[0] == 0 else '') +\
                     (self.factory.variable if item[0] > 0 else '') +\
                      ('^' + str(item[0])  if item[0] > 1 else '') \
                      for item in self.coefficients.iteritems()])
        aux = aux.replace('-1'+ self.factory.variable, '-'+self.factory.variable).replace('+ -','- ')
        if self.degree() > 0:
            aux = '(' + aux + ')'
        return aux
    def __call__(self, expression):
        subs = self.factory.inner_factory(expression)
        result = self.factory.inner_factory.zero()
        for k in xrange(self.degree(),-1,-1):
            scalar = self.coefficient(k)
            result = result*subs + scalar
        return result
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


class IntegerPolynomial(PolynomialOverIntegral):
    """docstring for IntegerPolynomial"""
    def __init__(self, coefficients, builder):
        super(IntegerPolynomial, self).__init__(coefficients, builder)
    def squarefree(self):
        Q_X = algebra_computacional.factories.FractionalPolynomialFactory()
        f = Q_X(str(self))
        f_prime = f.derivative()
        g = [f.gcd(f_prime)]
        h = [f/g[-1]]
        m = []
        while not h[-1].is_one():
            h.append(g[-1].gcd(h[-1]))
            g.append(g[-1]/h[-1])
            m.append(h[-2]/h[-1])
        if len(m) == 1:
            return [self]
        return [self.factory(str(x)) for x in m]

    def factors(self, rec=True, tuplify=False):
        if rec:
            first_factors = self.squarefree()
            i = 1
            result = []
            for x in first_factors:
                if tuplify:
                    for y in x.factors(rec=False):
                        result.append((y,i))
                else:
                    result = result + x.factors(rec=False)*i
                i+=1
            return result
        if self.degree() == 1:
            return [self]
        elif self.degree() == 0:
            return []
        for p in gen_primes():
            GFp = algebra_computacional.factories.GaloisPolynomialFactory(p)
            f_bar = GFp(str(self))
            f_bar_prime = f_bar.derivative()
            aux = f_bar.gcd(f_bar_prime)
            if not (self.leading_coefficient()%self.factory.inner_factory(p)).is_zero() and not f_bar_prime.is_zero() and f_bar.gcd(f_bar_prime).is_one():
               prime = p
               break
        f_bar_factors = f_bar.factors2()
        N = int(math.log((2**self.degree() + 1)*self.norm(),prime))
        len_factors = len(f_bar_factors)
        next_f_bar_factors = []
        f = f_bar
        lifted_f = lift(f_bar)[0]
        i = 1
        while i < N:         
            while len_factors > 0:
                g = f_bar_factors.pop(0)
                h = f/g
                gcd, s, t = g.eea(h)
                if gcd.is_one():
                    aux_f = h
                    len_factors -=1
                    aux, g, h, s, t = lift(f, g, h, s, t)
                    g, h, s, t = aux.hensel(g, h, s, t)
                    f = aux_f
                    next_f_bar_factors.insert(0,g)
                    if len_factors == 1:
                        next_f_bar_factors.insert(0,h)
                        break
                else:
                    f_bar_factors.append(g)
            f_bar_factors = next_f_bar_factors[:]
            len_factors = len(f_bar_factors)
            next_f_bar_factors = []
            f = lifted_f
            lifted_f = lift(f)[0]
            i*=2
        return self.true_factors(f_bar_factors,prime,N)

    def norm(self):
        res = 0
        for v in self.coefficients.itervalues():
            res += (v.number*v.number)
        return math.sqrt(res)
    def true_factors(self, factors, p, N):
        h = self
        L = factors[0].factory.monomial(0,str(self.leading_coefficient()))
        L_low = self.factory.monomial(0,str(self.leading_coefficient()))
        Result = []
        d = 1
        len_i = len(factors)
        I = set(range(len_i))
        while 2*d <= len_i:
            J = set(itertools.combinations(I, d))
            while J and 2*d <= len_i:
                S = J.pop()
                g_bar = L
                for g in [factors[i] for i in S]:
                    g_bar = g_bar*g
                g = self.factory({k: self.factory.inner_factory(str(v)) for k,v in g_bar.coefficients.iteritems()})
                if ((L_low*h) % g).is_zero():
                    gpp = g.pp()
                    Result.append(gpp)
                    h = h/gpp
                    I = I.difference(S)
                    J = J.difference([T for T in powerset(I) if set(T).intersection(S)])
                d += 1
        if h.degree() > 0:
            Result.append(h)
        return Result
def lift(*poly):
    new_GF = algebra_computacional.factories.GaloisPolynomialFactory(poly[-1].factory.q**2,lifted=True)
    return tuple([new_GF(str(p)) for p in poly])

def powerset(iterable):
    s = list(iterable)
    return itertools.chain.from_iterable(itertools.combinations(s, r) for r in range(len(s)+1))



class PolynomialOverField(PolynomialOverIntegral, EuclideanDomain):
    """docstring for Polynomial"""
    def __init__(self, coefficients, builder):
        super(PolynomialOverField, self).__init__(coefficients, builder)
    def __divmod__(self, rop):
        if rop.is_zero():
            raise ZeroDivisionError
        l1 = self.leading_coefficient()
        l2 = rop.leading_coefficient()
        if not (l1.is_one() and l2.is_one()):
            l1 = l1/l2
        l = self.factory.monomial(0,l1)
        q = self.factory.zero()
        r = self.monic()
        rop = rop.monic()

        while not r.is_zero() and r.degree() >= rop.degree():
            t = self.factory.monomial(r.degree() - rop.degree(), r.leading_coefficient() / rop.leading_coefficient())
            q = q+t
            r = r - (t*rop)
        return q*l, r*l
    def __div__(self, rop):
        return divmod(self, rop)[0]
    def __mod__(self, rop):
        return divmod(self, rop)[1]
    def monic(self):
        ld = self.leading_coefficient()
        return self.factory({ k: v/ld for k,v in self.coefficients.iteritems() })
        
    def gcd(self, b):
        aux = self.eea(b)[0]
        if self.factory.monic:
            return aux.monic()
        return aux
    def eea(self, b):
        g, s, t = EuclideanDomain.eea(self, b)
        if g.degree() == 0:
            s = s/g
            t = t/g
            g = self.factory.one();
        return g, s, t
class FractionalPolynomial(PolynomialOverField):
    def __init__(self, coefficients, builder):
        super(PolynomialOverField, self).__init__(coefficients, builder)

    def factors(self):
        zx_f = self.toZX()
        low_facs = zx_f.factors(tuplify=True)
        return [self.factory(str(x)) for x,i in low_facs]

    def toZX(self):
        ZX = algebra_computacional.factories.IntegerPolynomialFactory()
        if(self.is_zero()):
            return ZX.zero()
        denominators = [v.b for v in self.coefficients.itervalues()]
        lcm = normal_lcm(denominators)
        return ZX(str(self*self.factory(str(lcm))))



def list_power(l, exponent, one):
    return concat([[one() for _ in range(1,exponent)] + [x] for x in l])

class PolynomialOverGalois(PolynomialOverField, Field):
    def __init__(self, coefficients, builder):
        if isinstance(builder,algebra_computacional.factories.GaloisPolynomialFactory) and builder.lifted:
            coefficients = {k:v for k,v in coefficients.iteritems() if int(str(v)) == 0 or int(str(v)) == 1 or int(str(v)) == -1 or builder.p % int(str(v)) != 0}
            if not coefficients:
                coefficients = {0L:builder.inner_factory.zero()}
        super(PolynomialOverGalois, self).__init__(coefficients, builder)
    def inverse(self):
        return self.factory.one() / self
    def __divmod__(self, rop):
        if rop.is_zero():
            raise ZeroDivisionError
        l1 = self.leading_coefficient()
        l2 = rop.leading_coefficient()
        if not (l1.is_one() and l2.is_one()):
            l1 = l1/l2
        l = self.factory.monomial(0,l1)
        q = self.factory.zero()
        r = self.monic()
        rop = rop.monic()
        while not r.is_zero() and r.degree() >= rop.degree():
            t = self.factory.monomial(r.degree() - rop.degree(), r.leading_coefficient() * rop.leading_coefficient().inverse())
            q = q+t
            r = r - (t*rop)
        return q*l, r*l
    def monic(self):
        ld = self.leading_coefficient().inverse()
        return self.factory({ k: v*ld for k,v in self.coefficients.iteritems() })
    def pth_root(self):
        return self.factory({ k/self.factory.p: v for k,v in self.coefficients.iteritems()})
    
    def frob_basis(self):
        base_size = self.degree()
        image = []
        for i in range(base_size):
            canonic = self.factory.monomial(i)
            aux = (((canonic^(self.factory.q)) - canonic)%self).flip(base_size-1) * self.factory.monomial(base_size)
            image.append(aux + self.factory.monomial(base_size - 1 - i))
        return gaussian_elimination_polynomials(image,base_size-1)
    def factors1(self,tuplify=False):
        SQUAREFREE = self.squarefree()
        factors = []
        i = 1
        for sf in SQUAREFREE:
            if not sf.is_one():
                if tuplify:
                    for x in sf.bcz():
                        factors.append((x,i))
                else:
                    factors = factors + sf.bcz()*i
            i+=1
        return factors
    def factors2(self,tuplify=False):
        factors = []
        i = 1
        SQUAREFREE = self.squarefree()
        for sf in SQUAREFREE:
            if not sf.is_one():
                DISTINCTDEGREE = sf.distinctdegree()
                for dd in DISTINCTDEGREE:
                    if tuplify:
                        for x in dd[0].equaldegree(dd[1]):
                            factors.append((x,i))
                    else:
                        factors = factors + dd[0].equaldegree(dd[1])*i
            i+=1
        return factors
    def factors3(self,tuplify=False):
        factors = []
        SQUAREFREE = self.squarefree()
        i = 1
        for sf in SQUAREFREE:
            if not sf.is_one():
                if tuplify:
                    for x in sf.berlekamp():
                        factors.append((x,i))
                else:
                    factors = factors + sf.berlekamp()*i
            i+=1
        return factors
            

    def berlekamp(self):
        basis = self.frob_basis()
        factors = [self]
        irreducibles = []
        counter = 1
        len_basis = len(basis)
        while counter < len_basis:
            g = factors.pop()
            h = g.berlekamp_splitting(basis)
            if h != g:
                factors.append(h)
                factors.append(g/h)
                counter+=1
            else:
                irreducibles.append(g)
        return factors + irreducibles
    def bcz(self):
        basis = self.frob_basis()
        result = [self]
        counter = 1
        len_basis = len(basis)
        while counter < len_basis:
            idx = random.randrange(0,counter)
            g = result[idx]
            while(g.degree()<=1):
                idx = random.randrange(counter)
                g = result[idx]
            if idx < counter - 1:
                result[idx] = result.pop()
            else:
                result.pop()

            h = self.factory.zero()
            while h.is_zero():
                for x in basis:
                    h = h + self.factory.monomial(0,self.factory.inner_factory.random())*x
            aux = (h^((self.factory.q-1)/2)) - self.factory.one()
            w = g.gcd(aux)
            if not (w.is_one() or w == g):
                result.append(w)
                result.append(g/w)
                counter+=1
            else:
                result.append(g)
        return result

    def squarefree(self):
        derivative = self.derivative()
        if derivative.is_zero():
            return list_power(self.pth_root().squarefree(),self.factory.p, self.factory.one)
        L = []
        g = [self.gcd(derivative)]
        w = [self / g[0]]
        while not (w[-1].is_one()):
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
        while d <= ((g[d].degree() / 2) - 1):
            d = d + 1
            h.append((h[d-1]^(self.factory.q)) % g[d-1])
            aux = g[d-1].gcd(h[d] - x)
            if not aux.is_one():
                Result.append((aux,d))
            g.append(g[d-1]/aux)
        if not g[d].is_one():
            Result.append((g[d],g[d].degree()))
        return Result
    def equaldegree(self, k):
        def m_k(alpha, k):
            w = self.factory.k
            ret = alpha
            next_alpha =  alpha
            for i in range(1,w*k):
                next_alpha = next_alpha*next_alpha
                ret = ret + next_alpha
            return ret % self
        len_h = 1
        len_h_prime = 0
        r = self.degree()/k
        H = [self]
        H_prime = []
        while len_h < r:
            for h in H:
                a = self.factory.random() % h
                d = h.gcd(m_k(a,k))
                if d.is_one() or d == h:
                    H_prime.append(h)
                    len_h_prime+=1
                else:
                    H_prime.append(d)
                    H_prime.append(h/d)
                    len_h_prime+=2
            H = H_prime
            len_h = len_h_prime
            H_prime = []
            len_h_prime = 0
        return H
    def hensel(self, g, h, s, t):
        delta = self - g*h
        g_star = g*(self.factory.one() + ((s*delta)/h)) + t*delta
        h_star = h + ((s*delta) % h)
        sigma = s*g_star + t*h_star - self.factory.one()
        s_star = s - ((s*sigma)%h_star)
        t_star = (self.factory.one() - sigma)*t - g_star*((s*sigma)/h_star)
        return g_star, h_star, s_star, t_star
    def berlekamp_splitting(self, prebase):
        s = len(prebase)
        if s is 1:
            return self
        i = 1
        elements = self.factory.inner_factory.element_iterator()
        while i < s:
            while elements.has_next():
                element = self.factory.monomial(0,elements.next())
                g = self.gcd(prebase[i] - element)
                if not (g.is_one() or g == self):
                    return g
            i = i + 1
            elements.restart()
        return self

def gaussian_elimination_polynomials(image, max_degree):
    s = len(image)
    for i in range(s-1):
        if image[i].degree() < max_degree:
            found_bool = False
            for j in range(i+1,s):
                if image[j].degree() == max_degree:
                    image[i], image[j] = image[j], image[i]
                    found_bool = True
                    break
            if not found_bool:
                max_degree-=1
                continue
        image[i] = image[i].monic()
        for j in range(s):
            if j < i:
                if not image[j].coefficient(max_degree).is_zero():
                    image[j] = image[j] - image[j].factory.monomial(0,image[j].coefficient(max_degree)/image[i].leading_coefficient())*image[i]
            elif j > i:
                if image[j].degree() == image[i].degree():
                    image[j] = image[j] - image[j].factory.monomial(0,image[j].leading_coefficient()/image[i].leading_coefficient())*image[i]
        max_degree-=1
    return [x.flip(s-1) for x in image if x.degree() < s]

# Sieve of Eratosthenes
# Code by David Eppstein, UC Irvine, 28 Feb 2002
# http://code.activestate.com/recipes/117119/
def gen_primes():
    """ Generate an infinite sequence of prime numbers.
    """
    # Maps composites to primes witnessing their compositeness.
    # This is memory efficient, as the sieve is not "run forward"
    # indefinitely, but only as long as required by the current
    # number being tested.
    #
    D = {}
    
    # The running integer that's checked for primeness
    q = 2
    
    while True:
        if q not in D:
            # q is a new prime.
            # Yield it and mark its first multiple that isn't
            # already marked in previous iterations
            # 
            yield q
            D[q * q] = [q]
        else:
            # q is composite. D[q] is the list of primes that
            # divide it. Since we've reached q, we no longer
            # need it in the map, but we'll mark the next 
            # multiples of its witnesses to prepare for larger
            # numbers
            # 
            for p in D[q]:
                D.setdefault(p + q, []).append(p)
            del D[q]
        
        q += 1

class PolynomialOverExtensionFractional(PolynomialOverField, Field):
    def __init__(self, coefficients, builder):
        if isinstance(builder,algebra_computacional.factories.GaloisPolynomialFactory) and builder.lifted:
            coefficients = {k:v for k,v in coefficients.iteritems() if int(str(v)) == 0 or int(str(v)) == 1 or int(str(v)) == -1 or builder.p % int(str(v)) != 0}
        super(PolynomialOverExtensionFractional, self).__init__(coefficients, builder)
    def inverse(self):
        return self.factory.one() / self