"""Independent no-three-in-line exponential-bound referee.

Uses distinct row-pair boards, cycle rook polynomials, and literal
permutation-prefix averages. No producer imports or source reads.
"""
from collections import Counter, defaultdict
from fractions import Fraction as Q
from itertools import combinations, permutations
from math import comb, factorial
import hashlib
import json
import sys

sys.stdout.reconfigure(newline="\n")
gates = 0


def need(ok, message):
    global gates
    gates += 1
    if not ok:
        raise RuntimeError(message)


def falling(n, k):
    return factorial(n) // factorial(n-k) if 0 <= k <= n else 0


def parts(n, lower=2):
    if n == 0:
        yield ()
    for first in range(lower, n+1):
        for rest in parts(n-first, first):
            yield (first,) + rest


def rook(shape):
    coefficients = [1]
    for half in shape:
        size = 2*half
        local = [1] + [size*comb(size-k, k)//(size-k) for k in range(1, half+1)]
        product = [0]*(len(coefficients)+len(local)-1)
        for i, a in enumerate(coefficients):
            for j, b in enumerate(local):
                product[i+j] += a*b
        coefficients = product
    return coefficients


def mean_from_rooks(shape):
    n, r = sum(shape), rook(shape)
    return sum((Q((-1)**(k+1)*(k-2)*(2*n-k+1)*r[k], (k+1)*falling(n,k))
                for k in range(3,n+1)), Q(0))


def boards(n):
    degree, rows = [0]*n, []
    def visit(row):
        if row == n:
            yield tuple(rows)
            return
        for a,b in combinations([c for c in range(n) if degree[c]<2],2):
            degree[a] += 1
            degree[b] += 1
            if all(2-v <= n-row-1 for v in degree):
                rows.append((a,b))
                yield from visit(row+1)
                rows.pop()
            degree[a] -= 1
            degree[b] -= 1
    yield from visit(0)


def shape_of(rows):
    n = len(rows)
    neighbors = [set() for _ in range(2*n)]
    for r, cs in enumerate(rows):
        for c in cs:
            neighbors[r].add(n+c)
            neighbors[n+c].add(r)
    unseen, shape = set(range(2*n)), []
    while unseen:
        todo, seen = [min(unseen)], set()
        while todo:
            v = todo.pop()
            if v not in seen:
                seen.add(v)
                todo.extend(neighbors[v]-seen)
        unseen -= seen
        shape.append(len(seen)//2)
    return tuple(sorted(shape))


def defect(points):
    return sum(max(v-2,0) for v in Counter(c-r for r,c in points).values())


def board_controls():
    result = []
    for n in range(3,7):
        data = defaultdict(lambda: [0,0,0,[0]*(n+1)])
        for rows in boards(n):
            points = [(r,c) for r,cs in enumerate(rows) for c in cs]
            occupancy = Counter(c-r for r,c in points)
            f = defect(points)
            triple = sum(comb(v,3) for v in occupancy.values())
            need((f == 0) == (triple == 0), "same selected-direction zero event")
            shape = shape_of(rows)
            entry = data[shape]
            entry[0] += 1
            entry[1] += f
            entry[2] += (f == 0)
            for k in range(3,n+1):
                entry[3][k] += sum(comb(v,k) for v in occupancy.values() if v>=k)
            if n <= 5:
                full = sum((b[0]-a[0])*(c[1]-a[1]) == (c[0]-a[0])*(b[1]-a[1])
                           for a,b,c in combinations(points,3))
                need(full != 0 or f == 0, "actual target implies defect zero")
                for a,b in combinations(range(n),2):
                    swap = lambda v: b if v == a else a if v == b else v
                    need(abs(defect([(swap(r),c) for r,c in points])-f)<=4,
                         "literal row transposition range")
                    need(abs(defect([(r,swap(c)) for r,c in points])-f)<=4,
                         "literal column transposition range")
        need(sum(v[0] for v in data.values()) == {3:6,4:90,5:2040,6:67950}[n],
             "complete distinct-board count")
        for shape,(count,total,zero,moments) in sorted(data.items()):
            expected = mean_from_rooks(shape)
            need(Q(total,count) == expected, "independent exact rook mean")
            r = rook(shape)
            for k in range(3,n+1):
                need(Q(moments[k],count) == Q((2*n-k+1)*r[k],(k+1)*falling(n,k)),
                     "summed diagonal rook normalization")
            result.append(dict(n=n,shape=shape,boards=count,mean=str(expected),
                               zero_probability=str(Q(zero,count))))
    hostile = {(r,c) for r in range(5) for c in (r,(r+1)%5)}
    swapped = {(2 if r==0 else 0 if r==2 else r,c) for r,c in hostile}
    need((defect(hostile),defect(swapped)) == (5,1), "range four is attained")
    return result


def skeleton(shape):
    points, offset = [], 0
    for half in shape:
        points.extend((offset+r,offset+c) for r in range(half) for c in (r,(r+1)%half))
        offset += half
    return points


def conditional_range_controls():
    count, largest = 0, Q(0)
    for shape in [(3,),(2,2),(4,)]:
        n = sum(shape)
        points = skeleton(shape)
        prefix_data = defaultdict(lambda: [0,0])
        for rho in permutations(range(n)):
            for sigma in permutations(range(n)):
                f = defect([(rho[r],sigma[c]) for r,c in points])
                code = rho[:-1] + sigma[:-1]
                for j in range(len(code)+1):
                    entry = prefix_data[code[:j]]
                    entry[0] += f
                    entry[1] += 1
        children = defaultdict(list)
        for code,(total,number) in prefix_data.items():
            if code:
                children[code[:-1]].append(Q(total,number))
        for values in children.values():
            width = max(values)-min(values)
            need(width <= 4, "literal conditional mean range, not only increment bound")
            largest = max(largest,width)
            count += 1
        total,number = prefix_data[()]
        need(Q(total,number)==mean_from_rooks(shape), "prefix root mean agrees")
    return dict(nodes=count,maximum_width=str(largest))


def algebra_controls():
    for y in range(101):
        target = max(y-2,0)
        for K in range(2,31):
            partial = sum((-1)**(k+1)*(k-2)*comb(y,k) for k in range(3,min(y,K)+1))
            residual = (Q((-1)**(K+1)*((K-2)*y+2)*comb(y-2,K-1),K)
                        if y >= 2 and y-2 >= K-1 else 0)
            need(partial-target == residual, "closed alternating-truncation residual")
            need(partial <= target if K%2==0 else partial >= target,
                 "pointwise Bonferroni orientation")
    family_count = 0
    for n in range(2,19):
        for shape in parts(n):
            r = rook(shape)
            family_count += 1
            if n>=3:
                need(r[3] == 2*n*(n-2)*(2*n-5)//3, "universal three-rook coefficient")
            if n>=4:
                c4 = shape.count(2)
                need(r[4] == n*(n-3)*(2*n-5)*(2*n-7)//6+c4,
                     "four-rook coefficient retains c4")
                B4 = Q(2*n-5,3)-Q(2*(2*n-3)*r[4],5*falling(n,4))
                lower = Q((2*n-5)*(n*n+5*n-11),15*(n-1)*(n-2))
                lower -= Q(2*(2*n-3)*c4,5*n*(n-1)*(n-2)*(n-3))
                need(B4==lower and B4>=Q(2*n,15), "finite universal lower mean")
                need(mean_from_rooks(shape)>=B4, "alternating exact mean lower bound")
            if n<=12:
                for k in range(2,n+5):
                    q = Q(factorial(k)*r[k],falling(n,k)*2**k) if k<=n else Q(0)
                    if k<=n:
                        need(0<=1-q<=Q(k*(k-1),4*(n-1)), "column-collision union bound")
                    for L in range(n+1):
                        theta = Q(L,n)
                        moment = Q(falling(L,k)*factorial(k)*r[k],falling(n,k)**2) if k<=n else Q(0)
                        deficit = (2*theta)**k-moment
                        upper = 2**k*k*(k-1)*(theta**(k-1)/Q(2*n)+theta**k/Q(4*(n-1)))
                        need(0<=deficit<=upper,"all-length factorial deficit, including k>n")
    for n in range(4,101):
        z = n-4
        need(11*n*n-62*n+78 == 11*z*z+26*z+6, "positive finite lower-bound numerator")
    return family_count


def exponential_constants():
    # Exact Taylor interval for exp(2), with a geometric tail bound.
    last = Q(2**40,factorial(40))
    e2_lower = sum((Q(2**j,factorial(j)) for j in range(41)),Q(0))
    next_term = last*Q(2,41)
    e2_upper = e2_lower + next_term/(1-Q(2,42))
    inv_lower,inv_upper = 1/e2_upper,1/e2_lower
    cminus_upper = Q(13,6)*e2_upper-Q(21,2)*inv_lower+2
    cplus_upper = Q(13,6)*e2_upper+Q(29,2)*inv_upper-2
    need(cminus_upper < 17 and cplus_upper < 16, "exact uniform mean error constants")
    alpha_lower,alpha_upper = 1-5*inv_upper,1-5*inv_lower
    need(alpha_lower>Q(2,15), "strictly positive asymptotic mean")
    for n in range(4,19):
        for shape in parts(n):
            mu = mean_from_rooks(shape)
            need(n*alpha_upper-17 <= mu <= n*alpha_lower+16,
                 "exact mean inside stronger rational interval controls")
    # Hoeffding: two sets of n-1 reveals, each range width four.
    for n in range(4,101):
        ranges_square = 2*(n-1)*4**2
        need(Q(2,ranges_square) == Q(1,16*(n-1)), "conditional-range exponent")
        need(Q(2*n,15)**2/Q(16*(n-1)) == Q(n*n,900*(n-1)), "finite exponent constant")
    return dict(alpha_interval=[str(float(alpha_lower)),str(float(alpha_upper))],
                lower_error_bound=str(float(cminus_upper)),upper_error_bound=str(float(cplus_upper)))


def main():
    family_count = algebra_controls()
    census = board_controls()
    coupling = conditional_range_controls()
    constants = exponential_constants()
    print("STATUS: PASS; independent rook, board, conditional-prefix, and rational constant audit")
    print("BOARD UNIVERSE: all 70086 distinct boards n3..6; all literal transpositions and determinants n<=5")
    for row in census:
        print(json.dumps(row,sort_keys=True))
    print("ROOK CYCLE TYPES n2..18:",family_count)
    print("CONDITIONAL PREFIXES:",json.dumps(coupling,sort_keys=True))
    print("RATIONAL CONSTANT CERTIFICATE:",json.dumps(constants,sort_keys=True))
    print("EXACT FINITE BOUND: P(X=0)<=P(F=0)<=exp[-n^2/(900(n-1))], n>=4")
    print("ASYMPTOTIC RATE: alpha^2/16, alpha=1-5e^-2; no all-slope zero equivalence")
    print("SEMANTIC SHA256:",hashlib.sha256(json.dumps([census,coupling,family_count],sort_keys=True).encode()).hexdigest())
    print("ACTIVE GATES:",gates)


if __name__ == "__main__":
    main()
