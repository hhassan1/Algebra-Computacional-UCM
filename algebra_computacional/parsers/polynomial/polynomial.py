import ply.lex as lex
import ply.yacc as yacc
class PolyParser(object):
    def __init__(self):
        self.parsers = dict()
    def get_parser(self, ring):
        if ring in self.parsers:
            return self.parsers[ring]
        else:
            class InnerParser(object):
                """docstring for InnerParser"""
                tokens = (
                        'VAR','NUMBER',
                        'PLUS','MINUS','TIMES','DIVIDE','POWER',
                        'LPAREN','RPAREN', 'SCALAR',
                        )
                states = (('num','inclusive'),)
                def t_num_NUMBER(self,t):
                    r'\d+'
                    t.value = long(t.value)
                    if self.parenthesis is 0:
                        t.lexer.begin('INITIAL')
                    return t
                def t_LPAREN(self,t):
                    r'\('
                    self.parenthesis = self.parenthesis+1
                    return t
                def t_RPAREN(self,t):
                    r'\)'
                    self.parenthesis = self.parenthesis-1
                    return t
                def t_VAR(self,t):
                    r'x'
                    self.variable = True
                    return t
                def t_PLUS(self,t):
                    r'\+'
                    if self.parenthesis is 0:
                        self.variable = False
                        t.lexer.begin('INITIAL')
                    return t
                def t_MINUS(self,t):
                    r'-'
                    if self.parenthesis is 0:
                        self.variable = False
                        t.lexer.begin('INITIAL')
                    return t
                def t_TIMES(self,t):
                    r'\*'
                    if self.parenthesis is 0:
                        self.variable = False
                        t.lexer.begin('INITIAL')
                    return t
                def t_POWER(self,t):
                    r'\^'
                    if self.parenthesis is 0 and self.variable:
                        self.variable = False
                        t.lexer.begin('num')
                    return t
                def __init__(self):
                    super(InnerParser, self).__init__()
                    self.parenthesis = 0
                    self.variable = False
                    tokens = (
                            'VAR','NUMBER',
                            'PLUS','MINUS','TIMES','DIVIDE','POWER',
                            'LPAREN','RPAREN', 'SCALAR',
                            )
                    #self.t_VAR     = r'x'
                    #self.t_PLUS    = r'\+'
                    #self.t_MINUS   = r'-'
                    #self.t_TIMES   = r'\*'
                    #self.t_POWER   = r'\^'
                    #self.t_LPAREN  = r'\('
                    #self.t_RPAREN  = r'\)'
                    self.t_ignore = " \t"
                    self.t_SCALAR = ring.scalar_translator
                    self.result = dict()
                    self.lexer = lex.lex(object=self)

                    # Parsing rules

                    precedence = (
                        ('left','PLUS','MINUS'),
                        ('left','TIMES','DIVIDE'),
                        ('right','POWER'),
                        ('right','UMINUS'),
                        )

                    # dictionary of names
                    monomials = dict()

                    def p_polynomial_expr(t):
                        'polynomial : statement'

                    def p_statement_expr(t):
                        '''statement : statement PLUS monomial  
                                     | statement minmonomial 
                                     | monomial 
                                     | minmonomial'''

                    def p_ringarith_binop(t):
                        '''ringarith : ringarith PLUS   ringarith
                                 | ringarith MINUS  ringarith
                                 | ringarith TIMES  ringarith
                                 | ringarith DIVIDE ringarith
                                 | ringarith POWER  arith'''
                        if t[2] == '+'  : t[0] = t[1] + t[3]
                        elif t[2] == '-': t[0] = t[1] - t[3]
                        elif t[2] == '*': t[0] = t[1] * t[3]
                        elif t[2] == '/': t[0] = t[1] / t[3]
                        elif t[2] == '^': t[0] = t[1] ** t[3]

                    def p_ringarith_number(t):
                        'ringarith : SCALAR'
                    def p_arith_binop(t):
                        '''arith : arith PLUS   arith
                                 | arith MINUS  arith
                                 | arith TIMES  arith
                                 | arith DIVIDE arith
                                 | arith POWER  arith'''
                        if t[2] == '+'  : t[0] = t[1] + t[3]
                        elif t[2] == '-': t[0] = t[1] - t[3]
                        elif t[2] == '*': t[0] = t[1] * t[3]
                        elif t[2] == '/': t[0] = t[1] / t[3]
                        elif t[2] == '^': t[0] = t[1] ** t[3]

                    def p_arith_number(t):
                        'arith : NUMBER'
                        t[0] = t[1]
                    def p_arith_uminus(t):
                        'arith : MINUS arith %prec UMINUS'
                        t[0] = -t[2]
                    def p_minmonomial_uminus(t):
                        'minmonomial : MINUS monomial %prec UMINUS'
                        aux = -self.result[t[2]]
                        del self.result[t[2]]
                        self.result[t[2]] = aux
                        t[0] = t[2]
                    def p_arith_group(t):
                        'arith : LPAREN arith RPAREN'
                        t[0] = t[2]
                    def p_exparith_group(t):
                        'exparith : LPAREN arith RPAREN'
                        t[0] = t[2]
                    def p_exparith_number(t):
                        'exparith : NUMBER'
                        t[0] = t[1]
                    def p_exparith_one(t):
                        'exparith : '
                        t[0] = 1
                    def p_ringscalar_group(t):
                        'ringscalar : LPAREN ringarith RPAREN'
                        t[0] = t[2]
                    def p_ringscalar_number(t):
                        'ringscalar : SCALAR'
                        t[0] = t[1]
                    def p_ringscalar_one(t):
                        'ringscalar : '
                        t[0] = ring.one()
                    def p_monomial_powered(t):
                        'monomial : ringscalar VAR POWER exparith'
                        t[0] = t[4]
                        if t[4] in self.result:
                            self.result[t[4]] += t[1]
                        else:
                            self.result[t[4]] = t[1]
                    def p_monomial_exp1(t):
                        'monomial : ringscalar VAR'
                        t[0] = 1
                        if 1 in self.result:
                            self.result[1] += t[1]
                        else:
                            self.result[1] = t[1]
                    def p_monomial_scalar(t):
                        'monomial : ringscalar'
                        t[0] = 0
                        if 0 in self.result:
                            self.result[0] += t[1]
                        else:
                            self.result[0] = t[1]

                    def p_error(t):
                        print("Syntax error at '%s'" % t.value)
                    self.parser = yacc.yacc(start='polynomial')

                def t_newline(self,t):
                    r'\n+'
                    t.lexer.lineno += t.value.count("\n")
                    
                def t_error(self,t):
                    print "Illegal character '%s'" % t.value[0]
                    t.lexer.skip(1)
                def __call__(self, expression):
                    self.result.clear()
                    self.parser.parse(expression)
                    aux = {k: v for k, v in self.result.iteritems() if v != ring.zero()}
                    return  aux

            aux = InnerParser()
            self.parsers[ring] = aux
            return aux
parse = PolyParser().get_parser