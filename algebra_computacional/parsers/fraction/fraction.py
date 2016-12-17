import ply.lex as lex
import ply.yacc as yacc

class FractionParser(object):
    """docstring for FractionParser"""
    def __init__(self, cl, ring):
        self.ring = ring
        self.fraction_class = cl
    def build_parser(self):
        class InnerParser(object):
            """docstring for InnerParser"""
            def __init__(self, fraction_class, ring):
                super(InnerParser, self).__init__()
                tokens = ('PLUS','MINUS','TIMES','DIVIDE','LPAREN','RPAREN','SCALAR',)
                self.tokens = ('PLUS','MINUS','TIMES','DIVIDE','LPAREN','RPAREN','SCALAR',)
                # Tokens
                self.t_PLUS    = r'\+'
                self.t_MINUS   = r'-'
                self.t_TIMES   = r'\*'
                self.t_DIVIDE  = r'/'
                self.t_LPAREN  = r'\('
                self.t_RPAREN  = r'\)'
                self.t_SCALAR  = ring.scalar_translator
                self.result = None
                # Ignored characters
                self.t_ignore = " \t"
                # Build the lexer
                
                # Parsing rules
                precedence = (
                    ('left','PLUS','MINUS'),
                    ('left','TIMES','DIVIDE'),
                    ('right','UMINUS'),
                    )
                self.lexer = lex.lex(object=self)
                # dictionary of names
                def p_start(t):
                    'start : statement'
                    self.result = t[1]
                def p_statement_expr(t):
                    '''statement : statement PLUS statement  
                                 | statement MINUS statement
                                 | statement TIMES statement '''
                    if t[2] == '+': t[0] = t[1] + t[3]
                    elif t[2] == '-': t[0] = t[1] - t[3]
                    elif t[2] == '*': t[0] = t[1] * t[3]
                def p_statement_div(t):
                    '''statement : statement DIVIDE fraction '''
                    t[0] = t[1] / t[3]
                def p_statement_uni(t):
                    '''statement : fraction '''
                    t[0] = t[1]
                def p_fraction_divide(t):
                    'fraction : fraction DIVIDE fraction'
                    t[0] = t[1] / t[3]
                def p_fraction_scalar(t):
                    'fraction : SCALAR'
                    aux = (t[1],)
                    t[0] = fraction_class(aux)
                def p_fraction_state(t):
                    'fraction : LPAREN statement RPAREN'
                    t[0] = t[2]
                def p_fraction_minus(t):
                    'fraction : MINUS fraction %prec UMINUS'
                    t[0] = -t[2]
                def p_error(t):
                    print("Syntax error at '%s'" % t.value)
                self.parser = yacc.yacc()
            def t_newline(self, t):
                    r'\n+'
                    t.lexer.lineno += t.value.count("\n")
            def t_error(self, t):
                print("Illegal character '%s'" % t.value[0])
                t.lexer.skip(1)
            def parse(self, expression):
                self.parser.parse(expression, lexer=self.lexer)
                return self.result
        return InnerParser(self.fraction_class, self.ring)