"""Exact controls for the two-direction random saturated-board theorem.

No repository imports. Python + SymPy; all gates remain live under -O.
The all-size and uniform-remainder proofs are in the companion report.
"""
from collections import Counter
from fractions import Fraction as Q
from itertools import combinations, permutations
from math import comb, factorial
import hashlib
import json
import sys
import sympy as s

sys.stdout.reconfigure(newline="\n")
GATES = 0


def need(ok, label):
    global GATES
    GATES += 1
    if not ok:
        raise RuntimeError(label)


def fall(x, k):
    answer = 1
    for j in range(k):
        answer *= x-j
    return answer


def grid_lines(n):
    plus = [frozenset((r, r+d) for r in range(n) if 0 <= r+d < n)
            for d in range(1-n, n)]
    minus = [frozenset((r, e-r) for r in range(n) if 0 <= e-r < n)
             for e in range(2*n-1)]
    return plus, minus


def geometry(n, literal=False):
    lengths = [n-abs(d) for d in range(1-n, n)]
    a = [sum((n-abs(d))**2 for d in range(-r, n-r)) for r in range(n)]
    shared = sum((n-abs(c-r))**2*(n-abs(c+r-(n-1)))**2
                 for r in range(n) for c in range(n))
    overlap = 2*sum(v*v for v in a)
    if literal:
        plus, minus = grid_lines(n)
        rowcol = lambda T: ({r for r, c in T}, {c for r, c in T})
        overlap2 = shared2 = 0
        for T in plus:
            for U in minus:
                rt, ct = rowcol(T)
                ru, cu = rowcol(U)
                z = len(T & U)
                need(z <= 1, "distinct directions have at most one common cell")
                w = len(T)**2*len(U)**2
                overlap2 += w*(len(rt & ru)+len(ct & cu))
                shared2 += w*z
        need(overlap == overlap2, "row Fubini formula versus all line pairs")
        need(shared == shared2, "cell Fubini formula versus parity-sensitive line pairs")
    first = sum(t**3 for t in lengths)**2
    need(sum(t**3 for t in lengths) == n*n*(n*n+1)//2,
         "exact diagonal cube sum")
    coefficient = 8*(Q(first,n**8)-Q(overlap,n**7)+Q(shared,n**6))
    return overlap, shared, coefficient


def census(parts):
    n = sum(parts)
    edges, offset = [], 0
    for size in parts:
        for j in range(size):
            edges += [(offset+j, offset+j), (offset+j, offset+(j+1)%size)]
        offset += size
    need(len(set(edges)) == 2*n, "simple degree-two source")
    ps = list(permutations(range(n)))
    hist = Counter()
    joint11 = Q(0)
    witness = None
    plus, minus = grid_lines(n)
    T, U = plus[n-1], minus[n-1]
    R = C = n
    Z = n % 2
    need(len(T & U) == Z, "central crossing parity")
    for rho in ps:
        for sigma in ps:
            B = frozenset((rho[r], sigma[c]) for r,c in edges)
            dp, dm = Counter(c-r for r,c in B), Counter(r+c for r,c in B)
            sp = sum(comb(v,3) for v in dp.values() if v>=3)
            sm = sum(comb(v,3) for v in dm.values() if v>=3)
            X = sum((b[0]-a[0])*(c[1]-a[1]) == (b[1]-a[1])*(c[0]-a[0])
                    for a,b,c in combinations(B,3))
            need(sp+sm <= X, "literal integer determinant contains both directions")
            hist[sp,sm,X] += 1
            joint11 += len(B&T)*len(B&U)
            if sp+sm == 0 and X and witness is None:
                witness = tuple(sorted(B))
    total = factorial(n)**2
    need(sum(hist.values()) == total, "complete independent shore-label universe")
    expectation = lambda f: sum(Q(f(a,b,x)*v,total) for (a,b,x),v in hist.items())
    mean = Q(2*n-5,3)
    need(expectation(lambda a,b,x:a) == mean == expectation(lambda a,b,x:b),
         "both exact means")
    vp = expectation(lambda a,b,x:a*a)-mean*mean
    vm = expectation(lambda a,b,x:b*b)-mean*mean
    cross = expectation(lambda a,b,x:a*b)-mean*mean
    vt = expectation(lambda a,b,x:(a+b)**2)-4*mean*mean
    need(vp == vm and vt == vp+vm+2*cross, "exact covariance decomposition")
    p2, q2 = Q(4*n-6,n*(n-1)**2), Q(2,n*(n-1))
    literal11 = Z*Q(2,n)+(R+C-2*Z)*q2+(n*n-R-C+Z)*p2
    need(joint11/total == literal11, "shared-cell two-edge inclusion identity")
    if Z:
        false_edge_disjoint = (R+C)*q2+(n*n-R-C)*p2
        need(literal11 != false_edge_disjoint, "hostile: crossing cannot be discarded")
    return {"parts":parts,"labelings":total,"mean_each":str(mean),
            "variance_each":str(vp),"cross_covariance":str(cross),
            "variance_sum":str(vt),"P_sum_zero":str(expectation(lambda a,b,x:a+b==0)),
            "P_X_zero":str(expectation(lambda a,b,x:x==0)),
            "zero_sum_positive_X_witness":witness}


def main():
    z, theta, eta, omega, Z = s.symbols('z theta eta omega Z')
    for a in range(1,4):
        for b in range(1,4):
            k = a+b
            alltuples = fall(theta/z,a)*fall(eta/z,b)
            path = a*b*omega/z*(theta/z)**(a-1)*(eta/z)**(b-1)
            repeated = a*b*Z*(theta/z)**(a-1)*(eta/z)**(b-1)
            pk = (2*z)**k*(1+Q(k*(k-1),4)*z)
            pa = (2*z)**a*(1+Q(a*(a-1),4)*z)
            pb = (2*z)**b*(1+Q(b*(b-1),4)*z)
            # Removing repeated tuples from alltuples changes only order z^2;
            # keep it explicitly here to expose that first failed implication.
            expression = s.expand((alltuples-path-repeated)*pk + path*(2*z)**k/2
                                  + repeated*(2*z)**(k-1)-alltuples*pa*pb)
            want = a*b*2**(k-1)*(theta**a*eta**b
                                      +theta**(a-1)*eta**(b-1)*(Z-omega))
            need(s.expand(expression.coeff(z,1)-want)==0,"shared-cell covariance expansion")
    x,y = s.symbols('x y')
    f = s.Rational(1,3)+x-x*x
    ov = 2*s.integrate(f*f,(x,0,1))
    cell = 2*s.integrate(s.integrate((1-x)**2*(1-y)**2,(y,0,1-x)),(x,0,1))
    cross = 8*(s.Rational(1,4)-ov+cell)
    var = s.Rational(80,9)+2*cross
    need(ov == s.Rational(23,45),"overlap integral")
    need(cell == s.Rational(19,90),"shared-cell integral with Jacobian and parity")
    need(cross == -s.Rational(2,5),"cross-direction covariance coefficient")
    need(var == s.Rational(364,45),"two-direction variance")
    need(var/s.Rational(4,3)**2 == s.Rational(91,20),"one-sided zero-event bound")
    need(8*(s.Rational(1,4)-ov) != cross,"omitting crossings changes leading order")
    geometry_rows = []
    for n in range(2,41):
        row = geometry(n,literal=n<=15)
        geometry_rows.append((n,*map(str,row)))
    # Finite geometry identities are checks, not extrapolated proof obligations.
    for n in (60,100,200,400):
        row = geometry(n)
        need(abs(row[2]+Q(2,5)) < Q(20,n),"controlled geometry approach to integral")
        geometry_rows.append((n,*map(str,row)))
    bank = [census(p) for p in ((3,),(4,),(2,2),(5,),(2,3))]
    need(any(row['zero_sum_positive_X_witness'] for row in bank),
         "retained directions do not detect every failure")
    need(bank[1]['cross_covariance'] != bank[2]['cross_covariance'],
         "finite cycle information remains")
    payload = {'geometry':geometry_rows,'census':bank}
    print('STATUS: PASS; exact controls with separate uniform all-size proof')
    print('FIRST-ORDER MIXED FACTORIAL KERNEL: theta^a eta^b + theta^(a-1) eta^(b-1)(Z-omega)')
    print('INTEGRALS: product=1/4; row+column overlap=23/45; occupied crossing=19/90')
    print('COEFFICIENTS: cross=-2/5; variance=364/45; zero-event bound=91/20')
    for row in bank:
        print(json.dumps(row,sort_keys=True))
    print('GEOMETRY FINAL:',geometry_rows[-4:])
    print('SCOPE: random labeled degree-two carriers, uniformly over cycle types; no extremal theorem')
    print('SEMANTIC SHA256:',hashlib.sha256(json.dumps(payload,sort_keys=True).encode()).hexdigest())
    print('ACTIVE GATES:',GATES)


if __name__ == '__main__':
    main()
