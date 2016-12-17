from Field import Field

class Galois(Field):
    def __init__(self, prime, exponent):
        super(Galois, self).__init__()
        self.prime = prime
        self.exponent = exponent
