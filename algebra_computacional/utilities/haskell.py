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
def foldl(op, neutral, seq):
    result = neutral
    for x in seq:
        result = op(result,x)
    return result
def foldr(op, neutral, seq):
    result = neutral
    for x in reversed(seq):
        result = op(result,x)
    return result
def zipWith(operator, left_list, right_list):
    return [ operator(l_operand,r_operand) for (l_operand, r_operand) in zip(left_list, right_list)] 
def zipWithAll(operator, else_operator, left_list, right_list):
    d = len(left_list) - len(right_list)
    if d > 0:
        long_list = left_list
    else:
        long_list = right_list
        d = -d
    return zipWith(operator, left_list, right_list) + map(else_operator, long_list[-d:])
def concat(lists):
    result = []
    for item in reversed(lists):
        result += item
    return result
def eea(a,b):
    R = [a, b]
    S = [1, 0]
    T = [0, 1]
    i = 1
    while R[i] != 0:
        Q = R[i-1] / R[i]
        R.append(R[i-1] - Q*R[i])
        S.append(S[i-1] - Q*S[i])
        T.append(T[i-1] - Q*T[i])
        i = i + 1
    g = R[i-1]
    u, v = S[i-1], T[i-1]
    return g, u, v
def gcd(a, b):
    return eea(a, b)[0]
