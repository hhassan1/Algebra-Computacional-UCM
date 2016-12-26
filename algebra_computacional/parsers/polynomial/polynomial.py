__all__ = ['poly_parse']
import ply.lex as lex
import ply.yacc as yacc

class GenericPolynomialLexer(object):
    """docstring for InnerParser"""
    
    def t_newline(self,t):
        r'\n+'
        t.lexer.lineno += t.value.count("\n")
    def t_error(self,t):
        print "Illegal character '%s'" % t.value[0]
        t.lexer.skip(1)
    def __init__(self):
        super(GenericPolynomialLexer, self).__init__()
        self.tokens = (
            'VAR',
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
        #self.t_VAR    = None
        self.t_PLUS   = r'\+'
        self.t_MINUS  = r'-'
        self.t_TIMES  = r'\*'
        self.t_POWER  = r'\^'
        self.t_LPAREN = r'\('
        self.t_RPAREN = r'\)'
        self.t_ignore = " \t"
        #self.t_SCALAR = None
        self.lexer = lex.lex(object=self)
    def input(self, expression):
        self.lexer.input(expression)
    def token(self):
        return self.lexer.token()

class PolynomialParser(object):

    def __init__(self, poly_factory, lexer):
        super(PolynomialParser, self).__init__()
        self.poly_factory = poly_factory
        self.parser = yacc.yacc(start='statement', module=self, debug=False)
        self.lexer = lexer
        self.result = None
    def __call__(self, expression):
        self.result = None
        self.parser.parse(expression, lexer=self.lexer)

        return self.result
    tokens = (
            'VAR',
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
        self.result = self.poly_factory.zero()
    def p_statement_expr(self, t):
        'statement : polynomial'
        self.result = t[1]

    def p_monomial_scalar(self, t):
        'monomial : scalar'
        t[0] = self.poly_factory.monomial(0, t[1])
    def p_monomial_single(self, t):
        'monomial : scalar VAR'
        t[0] = self.poly_factory.monomial(1, t[1])
    def p_monomial_exp(self, t):
        'monomial : scalar VAR POWER number'
        t[0] = self.poly_factory.monomial(t[4], t[1])
    def p_monomial_single_one(self, t):
        'monomial : VAR'
        t[0] = self.poly_factory.monomial(1)
    def p_monomial_exp_one(self, t):
        'monomial : VAR POWER number'
        t[0] = self.poly_factory.monomial(t[3])

    def p_polynomial_single(self, t):
        'polynomial : monomial'
        t[0] = t[1]
    def p_polynomial_arith(self, t):
        '''polynomial : polynomial PLUS   polynomial
                      | polynomial MINUS  polynomial
                      | polynomial TIMES  polynomial
                      | polynomial DIVIDE polynomial '''   
        if t[2] == '+': t[0] = t[1] + t[3]
        elif t[2] == '-': t[0] = t[1] - t[3]
        elif t[2] == '*': t[0] = t[1] * t[3]
        elif t[2] == '/': t[0] = t[1] / t[3]
    def p_polynomial_silentproduct(self, t):
        'polynomial : LPAREN polynomial RPAREN monomial'
        t[0] = t[2] * t[4]
    def p_polynomial_power(self, t):
        'polynomial : LPAREN polynomial RPAREN POWER number'
        t[0] = t[2] ^ t[5]
    def p_polynomial_parenthesis(self, t):
        'polynomial : LPAREN polynomial RPAREN'
        t[0] = t[2]
    def p_polynomial_uminus(self, t):
        'polynomial : MINUS polynomial %prec UMINUS'
        t[0] = -t[2]

    def p_number_single(self, t):
        'number : NUMBER'
        t[0] = long(t[1])

    def p_scalar_single(self, t):
        '''scalar : SCALAR
                  | NUMBER '''
        t[0] = self.poly_factory.inner_factory(t[1])
    def p_scalar_power(self, t):
        'scalar : scalar POWER number'
        t[0] = t[1] ^ t[3]

    def p_error(self, t):
        print("Syntax error at '%s'" % t.value)

__parsers = dict()

def ParserGenerator(factory, variable):
    if factory not in __parsers:
        class SpecificPolynomialLexer(GenericPolynomialLexer):
            def __init__(self):
                self.t_SCALAR = factory.inner_factory.scalar_regex
                self.t_VAR = variable
                super(SpecificPolynomialLexer, self).__init__()
        __parsers[factory] = PolynomialParser(factory, SpecificPolynomialLexer())
    return __parsers[factory]
def poly_parse(factory, exp, variable):
    return ParserGenerator(factory, variable)(exp)
