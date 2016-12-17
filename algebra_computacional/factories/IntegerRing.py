from algebra_computacional.rings.Integer import Integer
from algebra_computacional.parsers.integer.integer import parse as int_parse
class IntegerRing(object):
    """docstring for IntegerRing"""
    @staticmethod
    def __call__(self, expression):
        super(IntegerRing, self).__init__()
        if isinstance(expression, basestring):
            self.number = int_parse(expression)
        elif isinstance(expression, (int, long)):
            self.number = long(expression)
    @staticmethod
    def zero():
        return Integer(0L)
    @staticmethod
    def one():
        return Integer(1L)

    