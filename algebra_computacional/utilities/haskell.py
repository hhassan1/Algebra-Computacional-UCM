def scanl(op, neutral, seq):
    result = [neutral]
    for x in seq:
        result.append(op(result[-1],x))
    return result
def scanr(op, neutral, seq):
    result = [neutral]
    for x in reversed(seq):
        result.append(op(result[-1],x))
    return result
def zipWith(operator, left_list, right_list):
    return [ operator(l_operand,r_operand) for (l_operand, r_operand) in zip(left_list, right_list)] 