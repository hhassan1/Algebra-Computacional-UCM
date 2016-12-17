# lextab.py. This file automatically created by PLY (version 3.9). Don't edit!
_tabversion   = '3.8'
_lextokens    = set(('RPAREN', 'DIVIDE', 'POWER', 'NUMBER', 'TIMES', 'PLUS', 'LPAREN', 'VAR', 'MINUS'))
_lexreflags   = 0
_lexliterals  = ''
_lexstateinfo = {'INITIAL': 'inclusive'}
_lexstatere   = {'INITIAL': [('(?P<t_NUMBER>\\d+)|(?P<t_newline>\\n+)|(?P<t_VAR>[a-zA-Z])|(?P<t_LPAREN>\\()|(?P<t_PLUS>\\+)|(?P<t_POWER>\\^)|(?P<t_RPAREN>\\))|(?P<t_TIMES>\\*)|(?P<t_MINUS>-)|(?P<t_DIVIDE>/)', [None, ('t_NUMBER', 'NUMBER'), ('t_newline', 'newline'), (None, 'VAR'), (None, 'LPAREN'), (None, 'PLUS'), (None, 'POWER'), (None, 'RPAREN'), (None, 'TIMES'), (None, 'MINUS'), (None, 'DIVIDE')])]}
_lexstateignore = {'INITIAL': ' \t'}
_lexstateerrorf = {'INITIAL': 't_error'}
_lexstateeoff = {}
