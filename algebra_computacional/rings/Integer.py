from algebra_computacional.structures.EuclideanDomain import Euclid

class Integer(Euclid):
    """docstring for Int"""
    def __init__(self, value):
        super(Integer, self).__init__()
        self.number = long(value)
    @staticmethod
    def zero():
        return Integer(0L)
    @staticmethod
    def one():
        return Integer(1L)
    def __eq__(self, op2):
        return self.number == op2.number
    def __add__(self, op2):
        return Integer(self.number + op2.number)
    def __sub__(self, op2):
        return Integer(self.number - op2.number)
    def __neg__(self):
        return Integer(-self.number)
    def __mul__(self, op2):
        return Integer(self.number * op2.number)
    def __str__(self):
        return str(self.number)
    def __repr__(self):
        return repr(self.number)
    def __div__(self, op2):
        return Integer(self.number // op2.number)

    @staticmethod
    def scalar_translator(term):
        r'\d+'
        term.value = Integer(long(term.value))
        return term
