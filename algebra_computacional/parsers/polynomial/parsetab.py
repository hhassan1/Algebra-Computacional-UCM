
# parsetab.py
# This file is automatically generated. Do not edit.
_tabversion = '3.8'

_lr_method = 'LALR'

_lr_signature = 'EF8B8F8753867C8A996712F85FF45B74'
    
_lr_action_items = {'RPAREN':([9,10,23,24,25,28,29,30,39,40,41,42,43,44,45,46,47,],[-11,15,-9,-10,-17,-8,-6,-7,47,-18,48,-15,-16,-14,-12,-13,-20,]),'DIVIDE':([9,10,23,24,25,28,29,30,39,40,41,42,43,44,45,46,47,],[-11,16,-9,-10,-17,-8,16,16,34,-18,34,-15,-16,-14,34,34,-20,]),'POWER':([9,10,13,23,24,25,28,29,30,39,40,41,42,43,44,45,46,47,],[-11,17,22,17,35,-17,17,17,17,35,-18,35,35,35,35,35,35,-20,]),'NUMBER':([17,22,26,27,33,34,35,36,37,38,],[25,31,25,25,25,25,25,25,25,25,]),'TIMES':([9,10,23,24,25,28,29,30,39,40,41,42,43,44,45,46,47,],[-11,18,-9,-10,-17,-8,18,18,36,-18,36,-15,-16,-14,36,36,-20,]),'SCALAR':([0,1,8,11,16,18,19,20,],[4,9,4,4,9,9,9,9,]),'PLUS':([0,2,3,4,5,7,8,9,10,11,12,13,14,15,21,22,23,24,25,28,29,30,31,32,39,40,41,42,43,44,45,46,47,48,],[-26,-4,-5,-25,11,-29,-26,-11,19,-26,-3,-28,-19,-24,-2,-23,-9,-10,-17,-8,-6,-7,-22,-27,37,-18,37,-15,-16,-14,-12,-13,-20,-21,]),'LPAREN':([0,8,11,17,22,26,27,33,34,35,36,37,38,],[1,1,1,26,33,26,26,26,26,26,26,26,26,]),'VAR':([0,4,7,8,11,15,],[-26,-25,13,-26,-26,-24,]),'MINUS':([0,2,3,4,5,7,8,9,10,11,12,13,14,15,17,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,],[8,-4,-5,-25,8,-29,-26,-11,20,-26,-3,-28,-19,-24,27,-2,-23,-9,-10,-17,27,27,-8,-6,-7,-22,-27,27,27,27,27,27,27,38,-18,38,-15,-16,-14,-12,-13,-20,-21,]),'$end':([0,2,3,4,5,6,7,8,11,12,13,14,15,21,22,31,32,48,],[-26,-4,-5,-25,-1,0,-29,-26,-26,-3,-28,-19,-24,-2,-23,-22,-27,-21,]),}

_lr_action = {}
for _k, _v in _lr_action_items.items():
   for _x,_y in zip(_v[0],_v[1]):
      if not _x in _lr_action:  _lr_action[_x] = {}
      _lr_action[_x][_k] = _y
del _lr_action_items

_lr_goto_items = {'arith':([17,26,27,33,34,35,36,37,38,],[24,39,40,41,42,43,44,45,46,]),'monomial':([0,8,11,],[2,14,21,]),'minmonomial':([0,5,],[3,12,]),'ringscalar':([0,8,11,],[7,7,7,]),'exparith':([22,],[32,]),'statement':([0,],[5,]),'polynomial':([0,],[6,]),'ringarith':([1,16,18,19,20,],[10,23,28,29,30,]),}

_lr_goto = {}
for _k, _v in _lr_goto_items.items():
   for _x, _y in zip(_v[0], _v[1]):
       if not _x in _lr_goto: _lr_goto[_x] = {}
       _lr_goto[_x][_k] = _y
del _lr_goto_items
_lr_productions = [
  ("S' -> polynomial","S'",1,None,None,None),
  ('polynomial -> statement','polynomial',1,'p_polynomial_expr','PolyParser.py',49),
  ('statement -> statement PLUS monomial','statement',3,'p_statement_expr','PolyParser.py',52),
  ('statement -> statement minmonomial','statement',2,'p_statement_expr','PolyParser.py',53),
  ('statement -> monomial','statement',1,'p_statement_expr','PolyParser.py',54),
  ('statement -> minmonomial','statement',1,'p_statement_expr','PolyParser.py',55),
  ('ringarith -> ringarith PLUS ringarith','ringarith',3,'p_ringarith_binop','PolyParser.py',58),
  ('ringarith -> ringarith MINUS ringarith','ringarith',3,'p_ringarith_binop','PolyParser.py',59),
  ('ringarith -> ringarith TIMES ringarith','ringarith',3,'p_ringarith_binop','PolyParser.py',60),
  ('ringarith -> ringarith DIVIDE ringarith','ringarith',3,'p_ringarith_binop','PolyParser.py',61),
  ('ringarith -> ringarith POWER arith','ringarith',3,'p_ringarith_binop','PolyParser.py',62),
  ('ringarith -> SCALAR','ringarith',1,'p_ringarith_number','PolyParser.py',70),
  ('arith -> arith PLUS arith','arith',3,'p_arith_binop','PolyParser.py',72),
  ('arith -> arith MINUS arith','arith',3,'p_arith_binop','PolyParser.py',73),
  ('arith -> arith TIMES arith','arith',3,'p_arith_binop','PolyParser.py',74),
  ('arith -> arith DIVIDE arith','arith',3,'p_arith_binop','PolyParser.py',75),
  ('arith -> arith POWER arith','arith',3,'p_arith_binop','PolyParser.py',76),
  ('arith -> NUMBER','arith',1,'p_arith_number','PolyParser.py',84),
  ('arith -> MINUS arith','arith',2,'p_arith_uminus','PolyParser.py',87),
  ('minmonomial -> MINUS monomial','minmonomial',2,'p_minmonomial_uminus','PolyParser.py',90),
  ('arith -> LPAREN arith RPAREN','arith',3,'p_arith_group','PolyParser.py',96),
  ('exparith -> LPAREN arith RPAREN','exparith',3,'p_exparith_group','PolyParser.py',99),
  ('exparith -> NUMBER','exparith',1,'p_exparith_number','PolyParser.py',102),
  ('exparith -> <empty>','exparith',0,'p_exparith_one','PolyParser.py',105),
  ('ringscalar -> LPAREN ringarith RPAREN','ringscalar',3,'p_ringscalar_group','PolyParser.py',108),
  ('ringscalar -> SCALAR','ringscalar',1,'p_ringscalar_number','PolyParser.py',111),
  ('ringscalar -> <empty>','ringscalar',0,'p_ringscalar_one','PolyParser.py',114),
  ('monomial -> ringscalar VAR POWER exparith','monomial',4,'p_monomial_powered','PolyParser.py',117),
  ('monomial -> ringscalar VAR','monomial',2,'p_monomial_exp1','PolyParser.py',124),
  ('monomial -> ringscalar','monomial',1,'p_monomial_scalar','PolyParser.py',131),
]
