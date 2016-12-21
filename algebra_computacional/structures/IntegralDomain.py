class Integral(object):
    """docstring for Integral"""
    def __init__(self):
        super(Integral, self).__init__()
    def builder(self):
        raise NotImplementedError
    def zero(self):
        return self.builder().zero()
    def one(self):
        return self.builder().one()
    def is_zero(self):
        raise NotImplementedError
    def is_one(self):
        raise NotImplementedError
    def __eq__(self, op2):
        raise NotImplementedError
    def __ne__(self, op2):
        return not self == op2
    def __add__(self, op2):
        raise NotImplementedError
    def __sub__(self, op2):
        raise NotImplementedError
    def __neg__(self):
        raise NotImplementedError
    def __mul__(self, op2):
        raise NotImplementedError
    def __pow__(self, op2):
        raise NotImplementedError
