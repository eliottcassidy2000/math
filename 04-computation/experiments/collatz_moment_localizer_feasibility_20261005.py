"""Exact finite moment-localizer and tail-LP floor certificates.

No Collatz orbit, ROOT certificate, or weight bank is read. Moment packets
are independent premises about the selected source's actual distribution.
"""
from fractions import Fraction as F
from itertools import product, combinations
import json

C = F(8, 9)


def need(ok, text):
    if not ok:
        raise ValueError(text)


def natural(x):
    need(type(x) is int and x >= 0, "exact natural required")


def exact(x):
    need(type(x) in (int, F), "exact rational required")
    return F(x)


def vector(v):
    need(type(v) is tuple and len(v) > 0, "nonempty exact tuple required")
    return tuple(exact(x) for x in v)


def dot(a, b):
    need(len(a) == len(b), "dimensions disagree")
    return sum((x*y for x, y in zip(a, b)), F(0))


def evaluate(poly, x):
    y = F(0)
    for coefficient in reversed(poly):
        y = y*x+coefficient
    return y


def kernel(m, j):
    natural(m); natural(j)
    t = 1 << abs(m-j)
    return F(4*t, (1+t)**2)


def localizer(moments):
    h = vector(moments)
    need(len(h) >= 2 and len(h) % 2 == 0 and h[0] == 1,
         "normalized moments through an odd degree required")
    need(all(0 <= x <= 1 for x in h), "bounded kernel moments required")
    d = len(h)//2-1
    return tuple(tuple(C*h[i+j]-h[i+j+1] for j in range(d+1))
                 for i in range(d+1))


def quadratic(matrix, a):
    return dot(a, tuple(dot(row, a) for row in matrix))


def psd(matrix):
    """Exact Schur-complement PSD decision, including singular matrices."""
    a = [list(row) for row in matrix]
    while a:
        if any(a[i][i] < 0 for i in range(len(a))):
            return False
        pivot = next((i for i in range(len(a)) if a[i][i] > 0), None)
        if pivot is None:
            return not any(x for row in a for x in row)
        indices = [i for i in range(len(a)) if i != pivot]
        value = a[pivot][pivot]
        a = [[a[i][j]-a[i][pivot]*a[pivot][j]/value for j in indices]
             for i in indices]
    return True


def solve(matrix, rhs):
    """Canonical exact solution: RREF pivot variables, free variables zero."""
    if not matrix:
        need(not rhs, "empty dimensions")
        return ()
    n = len(matrix[0])
    need(len(matrix) == len(rhs), "linear-system dimensions")
    a = [list(row)+[v] for row, v in zip(matrix, rhs)]
    pivots = []
    row = 0
    for col in range(n):
        pivot = next((i for i in range(row, len(a)) if a[i][col]), None)
        if pivot is None:
            continue
        a[row], a[pivot] = a[pivot], a[row]
        scale = a[row][col]
        a[row] = [x/scale for x in a[row]]
        for i in range(len(a)):
            if i != row:
                scale = a[i][col]
                a[i] = [x-scale*y for x, y in zip(a[i], a[row])]
        pivots.append(col)
        row += 1
        if row == len(a):
            break
    need(all(any(line[:n]) or line[n] == 0 for line in a),
         "inconsistent minimizer equations: packet fails necessary condition")
    answer = [F(0)]*n
    for i, col in enumerate(pivots):
        answer[col] = a[i][n]
    return tuple(answer)


def optimize_exact(m, moments):
    """Conditional optimum for P(1)=1; not a realizability validator.

    Every truthful packet has a PSD tangent Gram matrix and solvable normal
    equations. Passing these tests does NOT validate its oracle provenance.
    """
    natural(m)
    matrix = localizer(moments)
    d = len(matrix)-1
    gram = tuple(tuple(matrix[i][j]-matrix[i][0]-matrix[0][j]+matrix[0][0]
                       for j in range(1, d+1)) for i in range(1, d+1))
    linear = tuple(matrix[i][0]-matrix[0][0] for i in range(1, d+1))
    need(psd(gram), "packet fails necessary tangent PSD condition")
    z = solve(gram, tuple(-x for x in linear))
    a = (1-sum(z),)+z
    value = quadratic(matrix, a)
    need(sum(a) == 1, "normalization")
    need(all(dot(matrix[i], a) == dot(matrix[0], a) for i in range(1, d+1)),
         "normal equations")
    return {"source_index": m, "polynomial": a, "minimum": value,
            "conditional_floor": -value/(1-C), "degree": d}


def selector(a):
    a = vector(a)
    need(sum(a) == 1, "P(1)=1 is required")
    square = [F(0)]*(2*len(a)-1)
    for i, x in enumerate(a):
        for j, y in enumerate(a):
            square[i+j] += x*y
    q = [F(0)]*(len(square)+1)
    for i, x in enumerate(square):
        q[i] -= C*x/(1-C)
        q[i+1] += x/(1-C)
    return tuple(q)


def packet_intervals(intervals):
    need(type(intervals) is tuple and len(intervals) > 0, "interval tuple required")
    out = []
    for pair in intervals:
        need(type(pair) is tuple and len(pair) == 2, "exact interval pair")
        lo, hi = map(exact, pair)
        need(0 <= lo <= hi <= 1, "bounded ordered moment interval")
        out.append((lo, hi))
    need(out[0] == (F(1), F(1)), "moment zero is exactly one")
    return tuple(out)


def interval_floor(m, a, intervals):
    natural(m)
    q = selector(a)
    bounds = packet_intervals(intervals)
    need(len(bounds) == len(q), "one interval per selector coefficient")
    lower = sum((c*(lo if c >= 0 else hi) for c, (lo, hi) in zip(q, bounds)), F(0))
    return {"source_index": m, "conditional_floor": lower, "selector": q,
            "coefficient_norm": sum(map(abs, q), F(0))}


def exact_moments(law, maximum):
    return tuple(sum((mass*x**k for x, mass in law), F(0)) for k in range(maximum+1))


def finite_tail_model(m, head, intervals, tail_cap):
    """Outer LP on head masses plus one tail mass; all tail budgets retained."""
    natural(m); natural(head)
    need(head > m, "finite head must include target")
    intervals = packet_intervals(intervals)
    cap = exact(tail_cap)
    need(0 <= cap <= 1, "tail mass cap")
    size = head+1
    rows = [tuple(F(1) for _ in range(size)), tuple(F(-1) for _ in range(size))]
    rhs = [F(1), F(-1)]
    rows.append((F(0),)*head+(F(1),)); rhs.append(cap)
    u = kernel(m, head)
    for k, (lo, hi) in enumerate(intervals[1:], 1):
        values = tuple(kernel(m, j)**k for j in range(head))
        rows.append(values+(F(0),)); rhs.append(hi)
        rows.append(tuple(-x for x in values)+(-u**k,)); rhs.append(-lo)
    objective = tuple(F(int(j == m)) for j in range(size))
    return tuple(rows), tuple(rhs), objective


def verify_lp_floor(rows, rhs, objective, multipliers):
    """Check Ax<=b,x>=0 => objective*x>=-y*b by rational multiplication."""
    rhs = vector(rhs); objective = vector(objective); y = vector(multipliers)
    need(type(rows) is tuple and len(rows) == len(rhs) == len(y), "LP row dimensions")
    rows = tuple(vector(row) for row in rows)
    need(all(len(row) == len(objective) for row in rows), "LP column dimensions")
    need(all(v >= 0 for v in y), "nonnegative dual multipliers required")
    slack = tuple(objective[j]+sum((y[i]*rows[i][j] for i in range(len(rows))), F(0))
                  for j in range(len(objective)))
    need(all(v >= 0 for v in slack), "dual pointwise inequality fails")
    return -dot(y, rhs)


def unique_solution(matrix, rhs):
    """Independent square elimination for the bounded vertex-enumeration audit."""
    n = len(rhs)
    a = [list(row)+[v] for row, v in zip(matrix, rhs)]
    for col in range(n):
        pivot = next((i for i in range(col, n) if a[i][col]), None)
        if pivot is None:
            return None
        a[col], a[pivot] = a[pivot], a[col]
        scale = a[col][col]
        a[col] = [x/scale for x in a[col]]
        for row in range(n):
            if row != col:
                scale = a[row][col]
                a[row] = [x-scale*y for x, y in zip(a[row], a[col])]
    return tuple(a[i][-1] for i in range(n))


def finite_vertices(rows, rhs, dimension):
    """Exact small control only; caller supplies a bounded nonnegative polytope."""
    constraints = tuple(rows)+tuple(tuple(F(-int(i == j)) for j in range(dimension))
                                   for i in range(dimension))
    values = tuple(rhs)+(F(0),)*dimension
    vertices = set()
    for chosen in combinations(range(len(constraints)), dimension):
        x = unique_solution(tuple(constraints[i] for i in chosen),
                            tuple(values[i] for i in chosen))
        if x is not None and all(dot(row, x) <= b for row, b in zip(constraints, values)):
            vertices.add(x)
    return tuple(sorted(vertices))


def tail_lp_selector(m, head, a, intervals, tail_cap):
    """Compile a signed polynomial to a deliberately coarser finite LP dual."""
    q = selector(a)
    need(len(intervals) == len(q), "selector/packet dimensions")
    rows, rhs, objective = finite_tail_model(m, head, intervals, tail_cap)
    y = [F(0)]*len(rows)
    y[1 if q[0] >= 0 else 0] = abs(q[0])
    for k, coefficient in enumerate(q[1:], 1):
        upper, lower = 3+2*(k-1), 4+2*(k-1)
        y[lower if coefficient >= 0 else upper] = abs(coefficient)
    u = kernel(m, head)
    tail_upper = q[0]+sum((coefficient*u**k for k, coefficient in enumerate(q[1:], 1)
                          if coefficient >= 0), F(0))
    y[2] = max(F(0), tail_upper)
    floor = verify_lp_floor(rows, rhs, objective, tuple(y))
    return floor, tuple(y), (rows, rhs, objective)


def main():
    checks = 0
    def check(ok, message):
        nonlocal checks
        need(ok, message); checks += 1

    # Independent exact laws; no Collatz computation or ROOT inputs.
    laws = (
        ((F(1),F(1,100)), (C,F(1,2)), (F(16,25),F(49,100))),
        ((C,F(1,2)), (F(16,25),F(1,2))),
        ((F(1),F(1)),),
        ((F(16,25),F(1)),),
    )
    optimized_rows = []
    for law_id, law in enumerate(laws):
        target_mass = sum((p for x,p in law if x == 1), F(0))
        for degree in range(5):
            moments = exact_moments(law, 2*degree+1)
            result = optimize_exact(0, moments)
            a = result['polynomial']; q = selector(a)
            check(result['conditional_floor'] <= target_mass, "independent law lower bound")
            check(dot(q, moments) == result['conditional_floor'], "matrix versus polynomial")
            check(sum(a) == 1, "canonical normalization")
            for j in range(41):
                check(evaluate(q, kernel(0,j)) <= int(j == 0), "whole-support finite control")
            for small in product((-1,0,1), repeat=degree):
                trial = (F(1-sum(small)),)+tuple(map(F,small))
                check(quadratic(localizer(moments),trial) >= result['minimum'],
                      "independent finite competitor bound")
            optimized_rows.append((law_id,degree,str(result['conditional_floor']),tuple(map(str,a))))
    packet = exact_moments(laws[0],3)
    optimum = optimize_exact(0,packet)
    check(optimum['polynomial'] == (F(-16,9),F(25,9)), "same-order exact optimized vector")
    check(optimum['conditional_floor'] == F(1,100), "exact target recovered")
    monomials = tuple(9*packet[k+1]-8*packet[k] for k in range(3))
    check(all(x < 0 for x in monomials), "all earlier monomial tests fail at same order")

    # A tail of unknown law and explicit mass budget; it is not dropped.
    delta=F(1,100000); upper=kernel(0,3)
    head=(F(1,100),F(1,2),F(49,100)-delta)
    intervals=[(F(1),F(1))]
    for k in range(1,4):
        partial=sum((p*kernel(0,j)**k for j,p in enumerate(head)),F(0))
        intervals.append((partial,partial+delta*upper**k))
    intervals=tuple(intervals)
    bound=interval_floor(0,optimum['polynomial'],intervals)
    lp, y, model=tail_lp_selector(0,3,optimum['polynomial'],intervals,delta)
    check(0 < lp <= bound['conditional_floor'] <= F(1,100), "tail-safe finite LP floor")
    rows,rhs,objective=model
    actual_x=head+(delta,)
    check(all(dot(row,actual_x)<=b for row,b in zip(rows,rhs)), "outer-LP positive witness")
    check(dot(objective,actual_x)==F(1,100), "target is the actual LP objective")
    vertices = finite_vertices(rows,rhs,len(objective))
    check(bool(vertices), "bounded outer polytope is nonempty")
    lp_optimum = min(dot(objective,x) for x in vertices)
    check(lp <= lp_optimum <= F(1,100), "dual versus independent exact vertex optimum")
    zero_vertices = finite_vertices(rows+(objective,),rhs+(F(0),),len(objective))
    check(not zero_vertices, "zero-target relaxed model is infeasible")
    for tail_j in range(3,45):
        law=tuple((kernel(0,j),p) for j,p in enumerate(head))+((kernel(0,tail_j),delta),)
        h=exact_moments(law,3)
        check(all(lo<=v<=hi for (lo,hi),v in zip(intervals,h)), "all declared tail placements")
        check(lp<=dot(selector(optimum['polynomial']),h)<=F(1,100), "independent tail replay")

    # A missing-atom law has a deterministic canonical optimizer but no floor.
    absent=optimize_exact(0,exact_moments(laws[1],7))
    check(absent['conditional_floor']==0,"unique selected witness does not imply positive floor")
    # Continuum ghost: all localizer orders pass, but it is not an h-node law.
    ghost=F(7,10)
    check(F(16,25)<ghost<C,"ghost in forbidden discrete gap")
    check((ghost-F(16,25))*(ghost-C)<0,"discrete gap polynomial rejects ghost")
    for degree in range(5):
        hm=tuple(ghost**k for k in range(2*degree+2))
        check(psd(localizer(hm)),"relaxed continuum law passes localizer PSD")
        check(optimize_exact(0,hm)['conditional_floor']<=0,"no false positive from ghost")

    bad=(lambda: optimize_exact(True,packet),lambda: optimize_exact(0,(1,0.5)),
         lambda: selector((F(0),F(0))),lambda: interval_floor(0,(F(1),),((F(0),F(1)),(F(0),F(1)))),
         lambda: verify_lp_floor(rows,rhs,objective,tuple(-v for v in y)),
         lambda: finite_tail_model(3,3,intervals,delta))
    for call in bad:
        try: call()
        except ValueError: check(True,"typed or false-dual hostile")
        else: raise ValueError("invalid certificate accepted")
    print(json.dumps({'status':'PROVED conditional finite certificate; independent oracle obligation OPEN',
        'checks':checks,'same_order_optimized_floor':str(optimum['conditional_floor']),
        'same_order_monomial_readouts':list(map(str,monomials)),
        'optimized_polynomial':list(map(str,optimum['polynomial'])),
        'unknown_tail_mass':str(delta),'signed_interval_floor':str(bound['conditional_floor']),
        'finite_LP_floor':str(lp),'LP_dual_multipliers':list(map(str,y)),
        'finite_LP_exact_optimum':str(lp_optimum),'finite_LP_vertices':len(vertices),
        'missing_target_canonical_floor':str(absent['conditional_floor']),
        'optimized_exact_laws':optimized_rows,
        'Collatz_root_inputs':0},indent=2))


if __name__ == '__main__':
    main()
