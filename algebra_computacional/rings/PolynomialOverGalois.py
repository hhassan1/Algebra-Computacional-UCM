from algebra_computacional.structures.EuclideanDomain import Euclid
import PolynomialRing

class PolynomialOverGalois(PolynomialRing):
    def __init__(self, prime, exponent):
        super(PolynomialOverGalois, self).__init__()
        self.prime = prime
        self.exponent = exponent
    @classmethod
    def frob_basis(cls, polynomial, factors):
        image = []
        for canonic in [cls.monomial(i) for i in range(factors+1)]
            image.append((canonical^(prime^exponent) - canonical)%polynomial)
        result = gaussian_elimination(image)
        return [cls.one(), result]

    @classmethod
    def berlekamp(cls, polynomial):
        basis = frob_basis(polynomial)
        factors = [polynomial]
        irreducibles = []
        while len(factors) + len(irreducibles) < len(basis):
            g = basis.pop()
            h = berlekamp_splitting(g, prebase)
            if h != g:
                factors.append(h)
                factors.append(g/h)
            else:
                irreducibles.append(g)
        return factors + irreducibles
    @classmethod
    def squarefree(cls, polynomial):
        derivative = polynomial.derivative()
        if derivative.is_zero():
            return (squarefree(f^(1/p)))^p
        else:
            L = []
            g = [gcd(polynomial,derivative)]
            w = [f / g[0]]
            while not w[-1].is_one():
                w.append(gcd(g[-1],w[-1]))
                g.append(g[-1],w[-1])
                L.append(w[-2]/w[-1])
            if g[-1].is_one():
                return L
            return zipWith(cls.__mul__, L, squarefree(g[i]^(1/p))^p)
    @classmethod
    def distinctdegree(cls, polynomial):
        d = 0
        Result = []
        g = [polynomial]
        h = [x]
        while d<=(g[d].deg() / 2 - 1):
            d = d + 1
            h.append(h[d-1]^q % g[d-1])
            aux = gcd(g[d-1], h[d] - x)
            if not aux.is_one():
                Result.append((aux,d))
            g.append(g[d-1]/aux)
        if not g[d].is_one():
            Result.append((g[d],g[d].deg))
        return Result
    @classmethod
    def berlekamp_splitting(cls, polynomial, prebase):
        s = len(prebase)
        if s is 1:
            return polynomial
        i = 2
        while i <= s:
            for element in cls.elements:
                g = gcd(polynomial,prebase[i] - element)
                if not (g.is_one() or g == polynomial):
                    return g
            i = i + 1
    def __div__(self, op2):
        raise NotImplementedError
    def __invert__(self, op2):
        raise NotImplementedError
    