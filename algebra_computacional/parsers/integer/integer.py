__all__ = ['int_parse']
import ply.lex as lex
import ply.yacc as yacc
class GenericIntegerLexer(object):
    """docstring for InnerParser"""
    
    def t_newline(self,t):
        r'\n+'
        t.lexer.lineno += t.value.count("\n")
    def t_error(self,t):
        print "Illegal character '%s'" % t.value[0]
        t.lexer.skip(1)
    def __init__(self):
        super(GenericIntegerLexer, self).__init__()
        self.tokens = (
            'NUMBER',
            'PLUS',
            'MINUS',
            'TIMES',
            'DIVIDE',
            'POWER',
            'LPAREN',
            'RPAREN',
            'SCALAR',
            )
        self.t_NUMBER = r'\d+'
        self.t_PLUS   = r'\+'
        self.t_MINUS  = r'-'
        self.t_TIMES  = r'\*'
        self.t_DIVIDE  = r'/'
        self.t_POWER  = r'\^'
        self.t_LPAREN = r'\('
        self.t_RPAREN = r'\)'
        self.t_ignore = " \t"
        self.lexer = lex.lex(object=self)
    def input(self, expression):
        self.lexer.input(expression)
    def token(self):
        return self.lexer.token()
class IntegerParser(object):

    def __init__(self, factory, lexer):
        super(IntegerParser, self).__init__()
        self.factory = factory
        self.parser = yacc.yacc(start='statement', module=self, debug=False)
        self.lexer = lexer
        self.result = None
    def __call__(self, expression):
        self.result = None
        self.parser.parse(expression, lexer=self.lexer)

        return self.result
    tokens = (
            'NUMBER',
            'PLUS',
            'MINUS',
            'TIMES',
            'DIVIDE',
            'POWER',
            'LPAREN',
            'RPAREN',
            'SCALAR',
            )
    precedence = (
                    ('left','PLUS','MINUS'),
                    ('left','TIMES','DIVIDE'),
                    ('right','POWER'),
                    ('right','UMINUS'),
                 )
    def p_statement_empty(self, t):
        'statement : '
        self.result = self.factory.zero()
    def p_statement_expr(self, t):
        'statement : scalar'
        self.result = t[1]

    def p_scalar_arith(self, t):
        '''scalar : scalar PLUS   scalar
                  | scalar MINUS  scalar
                  | scalar TIMES  scalar
                  | scalar DIVIDE scalar ''' 
        if t[2] == '+': t[0] = t[1] + t[3]
        elif t[2] == '-': t[0] = t[1] - t[3]
        elif t[2] == '*': t[0] = t[1] * t[3]
        elif t[2] == '/': t[0] = t[1] / t[3]

    def p_number_single(self, t):
        'number : NUMBER'
        t[0] = long(t[1])

    def p_scalar_uminus(self, t):
        'scalar : MINUS scalar %prec UMINUS'
        t[0] = -t[2]
    def p_scalar_parenthesis(self, t):
        'scalar : LPAREN scalar RPAREN'
        t[0] = t[2]
    def p_scalar_single(self, t):
        '''scalar : SCALAR
                  | NUMBER '''
        t[0] = self.factory(long(t[1]))
    def p_scalar_power(self, t):
        'scalar : scalar POWER number'
        t[0] = t[1] ^ t[3]

    def p_error(self, t):
        print("Syntax error at '%s'" % t.value)

__parsers = dict()

def ParserGenerator(factory):
    if factory not in __parsers:
        class SpecificIntegerLexer(GenericIntegerLexer):
            def __init__(self):
                self.t_SCALAR = factory.scalar_regex
                super(SpecificIntegerLexer, self).__init__()
        __parsers[factory] = IntegerParser(factory, SpecificIntegerLexer())
    return __parsers[factory]
def int_parse(factory, exp):
    return ParserGenerator(factory)(exp)