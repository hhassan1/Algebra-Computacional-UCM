import algebra_computacional.factories.IntegerRing
import ply.lex as lex
import ply.yacc as yacc
class IntegerParser(object):
    """docstring for InnerParser"""
    tokens = ('PLUS', 'NUMBER', 'MINUS','TIMES','POWER','LPAREN','RPAREN')
    def __init__(self):
        super(IntegerParser, self).__init__()
        tokens = ('PLUS', 'NUMBER', 'MINUS','TIMES','POWER','LPAREN','RPAREN')
        self.t_PLUS    = r'\+'
        self.t_MINUS   = r'-'
        self.t_TIMES   = r'\*'
        self.t_POWER   = r'\^'
        self.t_LPAREN  = r'\('
        self.t_RPAREN  = r'\)'
        self.result = None
        self.t_ignore = " \t"
        self.r_term = r'\d+'
        precedence = (
            ('left','PLUS','MINUS'),
            ('left','TIMES'),
            ('right','POWER'),
            ('right','UMINUS'),
            )
        self.lexer = lex.lex(object=self)
        def p_num_number(t):
            'num : NUMBER'
            t[0] = t[1]
        def p_num_uminus(t):
            'num : MINUS num %prec UMINUS'
            t[0] = -t[2]
        def p_num_group(t):
            'num : LPAREN num RPAREN'
            t[0] = t[2]
        def p_error(t):
            print("Syntax error at '%s'" % t.value)
        def p_scalar(t):
            'scalar : num'
            self.result = IntegerRing(t[1])
        def p_num_binop(t):
            '''num : num PLUS   num
                   | num MINUS  num
                   | num TIMES  num
                   | num POWER  num'''
            if t[2] == '+'  : t[0] = t[1] + t[3]
            elif t[2] == '-': t[0] = t[1] - t[3]
            elif t[2] == '*': t[0] = t[1] * t[3]
            elif t[2] == '^': t[0] = t[1] ** t[3]
        self.parser = yacc.yacc(start='scalar')
    def t_NUMBER(self, t):
        r'\d+'
        t.value = long(t.value)
        return t
    def t_newline(self, t):
        r'\n+'
        t.lexer.lineno += t.value.count("\n")
    def t_error(self, t):
        print("Illegal character '%s'" % t.value[0])
        t.lexer.skip(1)
    
    def parse(self, expression):
        self.parser.parse(expression, lexer=self.lexer)
        return self.result
parse = IntegerParser().parse