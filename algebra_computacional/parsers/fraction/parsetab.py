
# parsetab.py
# This file is automatically generated. Do not edit.
_tabversion = '3.8'

_lr_method = 'LALR'

_lr_signature = 'E3A7E965D4089A53514685F4C0B6CF1B'
    
_lr_action_items = {'RPAREN':([2,3,12,13,14,15,16,17,18,19,],[-8,-6,19,-10,-7,-5,-4,-2,-3,-9,]),'DIVIDE':([2,3,4,12,13,14,15,16,17,18,19,],[-8,7,8,8,-10,-7,-5,-4,8,8,-9,]),'TIMES':([2,3,4,12,13,14,15,16,17,18,19,],[-8,-6,9,9,-10,-7,-5,-4,9,9,-9,]),'SCALAR':([0,5,6,7,8,9,10,11,],[2,2,2,2,2,2,2,2,]),'PLUS':([2,3,4,12,13,14,15,16,17,18,19,],[-8,-6,10,10,-10,-7,-5,-4,-2,-3,-9,]),'LPAREN':([0,5,6,7,8,9,10,11,],[5,5,5,5,5,5,5,5,]),'MINUS':([0,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,],[6,-8,-6,11,6,6,6,6,6,6,6,11,-10,-7,-5,-4,-2,-3,-9,]),'$end':([1,2,3,4,13,14,15,16,17,18,19,],[0,-8,-6,-1,-10,-7,-5,-4,-2,-3,-9,]),}

_lr_action = {}
for _k, _v in _lr_action_items.items():
   for _x,_y in zip(_v[0],_v[1]):
      if not _x in _lr_action:  _lr_action[_x] = {}
      _lr_action[_x][_k] = _y
del _lr_action_items

_lr_goto_items = {'start':([0,],[1,]),'fraction':([0,5,6,7,8,9,10,11,],[3,3,13,14,15,3,3,3,]),'statement':([0,5,9,10,11,],[4,12,16,17,18,]),}

_lr_goto = {}
for _k, _v in _lr_goto_items.items():
   for _x, _y in zip(_v[0], _v[1]):
       if not _x in _lr_goto: _lr_goto[_x] = {}
       _lr_goto[_x][_k] = _y
del _lr_goto_items
_lr_productions = [
  ("S' -> start","S'",1,None,None,None),
  ('start -> statement','start',1,'p_start','FracParser.py',38),
  ('statement -> statement PLUS statement','statement',3,'p_statement_expr','FracParser.py',41),
  ('statement -> statement MINUS statement','statement',3,'p_statement_expr','FracParser.py',42),
  ('statement -> statement TIMES statement','statement',3,'p_statement_expr','FracParser.py',43),
  ('statement -> statement DIVIDE fraction','statement',3,'p_statement_div','FracParser.py',48),
  ('statement -> fraction','statement',1,'p_statement_uni','FracParser.py',51),
  ('fraction -> fraction DIVIDE fraction','fraction',3,'p_fraction_divide','FracParser.py',54),
  ('fraction -> SCALAR','fraction',1,'p_fraction_scalar','FracParser.py',57),
  ('fraction -> LPAREN statement RPAREN','fraction',3,'p_fraction_state','FracParser.py',61),
  ('fraction -> MINUS fraction','fraction',2,'p_fraction_minus','FracParser.py',64),
]