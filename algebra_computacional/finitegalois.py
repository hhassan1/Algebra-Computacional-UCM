import numpy as np
from copy import deepcopy

class FiniteGalois(Field):
    """docstring for FiniteGalois"""
    def __init__(self, p, k, a= None):
        super(FiniteGalois, self).__init__()
        self.p = p
        self.k = k
        self.a = a
    def gcd(self, a, b):
        M = np.empty((4,a.degree+ b.degree + 1))
        Q[:2] = [self.zero,self.zero]
        R[:2] = [a,b]
        S[:2] = [self.one, self.zero]
        T[:2] = [self.zero, self.one]
        i = 1
        while R[i] != self.zero:
            Q.append(R[i-1] // R[i])
            R.append(R[i-1] - Q[i+1]*R[i])
            S.append(S[i-1] - Q[i+1]*S[i])
            T.append(T[i-1] - Q[i+1]*T[i])
            i = i + 1
        return R[i-1], S[i-1], T[i-1], -(i % 2)*S[i], (1-(i % 2))*T[i]
    def poly(self, expression):
        coefficienparse_string(expression)
        
    class Polynomial(object):
        """docstring for Polynomial"""
        def __init__(self, coefs, offset = 0):
            super(FiniteGalois.Polynomial, self).__init__()
            self.content = coefs[0]
            for x in xrange(1,len(coefs)):
                self.content = self.gcd(self.content,abs(coefs[x]))
            self.offset = offset
            self.principal = [coefs[x]/self.content for x in range(len(coefs))]
            self.degree = len(coefs) - 1
        def __mul__(a,b):
            if b.degree > a.degree:
                a, b = b, a
            prin_product = []
            for x in range(a.degree+b.degree+2):
                prin_product.append(0)
                for i in xrange(max(x-a.degree,0),min(x, b.degree)+1),:
                    prin_product[x] += a.principal[i]*b.principal[x-i]
            aux = FiniteGalois.Polynomial(prin_product,a.offset+b.offset)
            aux.content *= a.content + b.content
            return aux
        def __add__(a,b):
            if b.degree > a.degree:
                a, b = b, a
            aux = []
            
            if a.offset > b.offset:
                min_offset, max_offset = b, a    
            else:
                min_offset, max_offset = a, b
            for i in xrange(min_offset.offset, max_offset.offset):
                aux.append(min_offset.content*min_offset.principal[i-min_offset])
            for i in range(b.degree+1):
                aux.append(min_offset.content*min_offset.principal[i+min_offset]+max_offset.content*max_offset.principal[i])
            for i in range(b.degree+1, a.degree+1):
                aux.append(b.content*b.principal[i])
            return FiniteGalois.Polynomial(aux,min_offset.offset)
        def __neg__(a):
            aux = FiniteGalois.Polynomial(map(lambda x: -x, a.principal), a.offset)
            aux.content = a.content
            return aux
        def __sub__(a,b):
            return a + (-b)
        def __floordiv__(a,b):
            if b == 0:
                raise ZeroDivisionError()
            q, r = 0, a
            while r != 0 and (r.degree+r.offset) >= (q.degree+q.offset):
                t = r.content*r.principal[r.degree] - d.content*d.principal[d.degree]
                q, r = q.add_constant(t), r - (d.scaled(t))
            return q, r
        def add_constant(self, t):
            aux = FiniteGalois.Polynomial(deepcopy(self.principal),0)
            for i in range(self.offset):
                aux.principal.insert(0,0)
            aux.principal[0] = t
            return aux
        def scale(self, t):
            aux = FiniteGalois.Polynomial(deepcopy(self.principal),self.offset)
            aux.content *= t
            return aux

def euclid(a, b):
    if a < b:
        (a, b) = (b, a)
    if b == 0:
        return a
    return euclid(b, a % b)
def extended_euclid(a,b):
    M = np.empty((4,np.ceil(np.log2(a)+np.log2(b))+1))
    M[:,0] = [0,a,1,0]
    M[:,1] = [0,b,0,1]
    i = 1
    while(M[1,i] != 0):
        M[0,i+1] = M[1,i-1] // M[1,i]
        M[1:,i+1] = M[1:,i-1] - M[0,i+1]*M[1:,i]
        i = i + 1
    return M[1,i-1], M[2,i-1], M[3,i-1], M[:,:i+1].T

def main():
    X = FiniteGalois(2,3)
    a = X.Polynomial([1,2])

    

if __name__ == '__main__':
    main()
    pass
