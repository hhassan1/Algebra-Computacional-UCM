
# parsetab.py
# This file is automatically generated. Do not edit.
_tabversion = '3.10'

_lr_method = 'LALR'

_lr_signature = 'EBD3F2ECB97163A7C81C3FC43E6BFB88'
    
_lr_action_items = {'RPAREN':([1,2,4,5,9,10,12,18,19,20,21,23,24,25,26,27,29,30,31,],[-18,-8,-6,-3,-19,19,-4,-16,-15,-7,-17,-20,-12,-9,-10,-11,-13,-5,-14,]),'DIVIDE':([1,2,4,5,7,9,10,12,18,19,20,21,23,24,25,26,27,29,30,31,],[-18,-8,-6,-3,14,-19,14,-4,-16,-15,-7,-17,-20,-12,14,14,-11,-13,-5,-14,]),'POWER':([1,4,5,9,12,19,21,23,],[-18,11,13,-19,22,28,-17,-20,]),'NUMBER':([0,3,8,11,13,14,15,16,17,19,22,28,],[9,9,9,21,21,9,9,9,9,9,21,21,]),'TIMES':([1,2,4,5,7,9,10,12,18,19,20,21,23,24,25,26,27,29,30,31,],[-18,-8,-6,-3,17,-19,17,-4,-16,-15,-7,-17,-20,-12,17,17,-11,-13,-5,-14,]),'SCALAR':([0,3,8,14,15,16,17,19,],[1,1,1,1,1,1,1,1,]),'PLUS':([1,2,4,5,7,9,10,12,18,19,20,21,23,24,25,26,27,29,30,31,],[-18,-8,-6,-3,15,-19,15,-4,-16,-15,-7,-17,-20,-12,-9,-10,-11,-13,-5,-14,]),'LPAREN':([0,3,8,14,15,16,17,],[3,3,3,3,3,3,3,]),'VAR':([0,1,3,5,8,9,14,15,16,17,19,21,23,],[4,-18,4,12,4,-19,4,4,4,4,4,-17,-20,]),'MINUS':([0,1,2,3,4,5,7,8,9,10,12,14,15,16,17,18,19,20,21,23,24,25,26,27,29,30,31,],[8,-18,-8,8,-6,-3,16,8,-19,16,-4,8,8,8,8,-16,-15,-7,-17,-20,-12,-9,-10,-11,-13,-5,-14,]),'$end':([0,1,2,4,5,6,7,9,12,18,19,20,21,23,24,25,26,27,29,30,31,],[-1,-18,-8,-6,-3,0,-2,-19,-4,-16,-15,-7,-17,-20,-12,-9,-10,-11,-13,-5,-14,]),}

_lr_action = {}
for _k, _v in _lr_action_items.items():
   for _x,_y in zip(_v[0],_v[1]):
      if not _x in _lr_action:  _lr_action[_x] = {}
      _lr_action[_x][_k] = _y
del _lr_action_items

_lr_goto_items = {'polynomial':([0,3,8,14,15,16,17,],[7,10,18,24,25,26,27,]),'monomial':([0,3,8,14,15,16,17,19,],[2,2,2,2,2,2,2,29,]),'scalar':([0,3,8,14,15,16,17,19,],[5,5,5,5,5,5,5,5,]),'number':([11,13,22,28,],[20,23,30,31,]),'statement':([0,],[6,]),}

_lr_goto = {}
for _k, _v in _lr_goto_items.items():
   for _x, _y in zip(_v[0], _v[1]):
       if not _x in _lr_goto: _lr_goto[_x] = {}
       _lr_goto[_x][_k] = _y
del _lr_goto_items
_lr_productions = [
  ("S' -> statement","S'",1,None,None,None),
  ('statement -> <empty>','statement',0,'p_statement_empty','polynomial.py',77),
  ('statement -> polynomial','statement',1,'p_statement_expr','polynomial.py',80),
  ('monomial -> scalar','monomial',1,'p_monomial_scalar','polynomial.py',84),
  ('monomial -> scalar VAR','monomial',2,'p_monomial_single','polynomial.py',87),
  ('monomial -> scalar VAR POWER number','monomial',4,'p_monomial_exp','polynomial.py',90),
  ('monomial -> VAR','monomial',1,'p_monomial_single_one','polynomial.py',93),
  ('monomial -> VAR POWER number','monomial',3,'p_monomial_exp_one','polynomial.py',96),
  ('polynomial -> monomial','polynomial',1,'p_polynomial_single','polynomial.py',100),
  ('polynomial -> polynomial PLUS polynomial','polynomial',3,'p_polynomial_arith','polynomial.py',103),
  ('polynomial -> polynomial MINUS polynomial','polynomial',3,'p_polynomial_arith','polynomial.py',104),
  ('polynomial -> polynomial TIMES polynomial','polynomial',3,'p_polynomial_arith','polynomial.py',105),
  ('polynomial -> polynomial DIVIDE polynomial','polynomial',3,'p_polynomial_arith','polynomial.py',106),
  ('polynomial -> LPAREN polynomial RPAREN monomial','polynomial',4,'p_polynomial_silentproduct','polynomial.py',112),
  ('polynomial -> LPAREN polynomial RPAREN POWER number','polynomial',5,'p_polynomial_power','polynomial.py',115),
  ('polynomial -> LPAREN polynomial RPAREN','polynomial',3,'p_polynomial_parenthesis','polynomial.py',118),
  ('polynomial -> MINUS polynomial','polynomial',2,'p_polynomial_uminus','polynomial.py',121),
  ('number -> NUMBER','number',1,'p_number_single','polynomial.py',125),
  ('scalar -> SCALAR','scalar',1,'p_scalar_single','polynomial.py',129),
  ('scalar -> NUMBER','scalar',1,'p_scalar_single','polynomial.py',130),
  ('scalar -> scalar POWER number','scalar',3,'p_scalar_power','polynomial.py',133),
]
