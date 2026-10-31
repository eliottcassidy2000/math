#!/usr/bin/env python3
"""Independent exact audit; no imports from either author packet.

Ring coefficients below are in Q(beta).  Planar vectors use basis (1,alpha),
where alpha is a nonreal conjugate, NOT the ring generator beta.
"""
from fractions import Fraction as F
from itertools import product
import json

checks = 0

def need(ok, message):
    global checks
    checks += 1
    if not ok:
        raise RuntimeError(message)

def rem_mod(a, b, p):
    a = [v % p for v in a]
    while len(a) >= len(b):
        v = a[-1]  # divisor is monic
        n = len(a) - len(b)
        for j, x in enumerate(b):
            a[n+j] = (a[n+j]-v*x) % p
        while a and not a[-1]:
            a.pop()
    return a

irreducibility = []
for k, prime in ((3, 3), (4, 2), (5, 5)):
    divisors = 0
    for degree in range(1, k//2+1):
        for coeff in product(range(prime), repeat=degree):
            need(bool(rem_mod([-1]*k+[1], list(coeff)+[1], prime)),
                 'proper modular divisor')
            divisors += 1
    irreducibility.append({'degree': k, 'prime': prime,
                           'all_monic_divisors_checked': divisors})

class Q:
    def __init__(self, x=0):
        if isinstance(x, Q):
            self.c = x.c
        elif isinstance(x, (tuple, list)):
            self.c = tuple(F(v) for v in x)
        else:
            self.c = (F(x), F(0), F(0))
    def __add__(self, x):
        return Q(tuple(a+b for a,b in zip(self.c, Q(x).c)))
    __radd__ = __add__
    def __neg__(self): return Q(tuple(-a for a in self.c))
    def __sub__(self, x): return self+-Q(x)
    def __rsub__(self, x): return Q(x)+-self
    def __mul__(self, x):
        d = [F(0)]*5
        for i,a in enumerate(self.c):
            for j,b in enumerate(Q(x).c): d[i+j] += a*b
        for i in (4,3):
            for j in range(3): d[i-3+j] += d[i]
        return Q(d[:3])
    __rmul__ = __mul__
    def __eq__(self, x): return self.c == Q(x).c
    def __bool__(self): return any(self.c)
    def inverse(self):
        columns = [(self*Q([int(j==i) for j in range(3)])).c for i in range(3)]
        a = [[columns[j][i] for j in range(3)]+[F(i==0)] for i in range(3)]
        for j in range(3):
            pivot = next(i for i in range(j,3) if a[i][j])
            a[j],a[pivot] = a[pivot],a[j]
            v = a[j][j]; a[j] = [x/v for x in a[j]]
            for i in range(3):
                if i != j:
                    v = a[i][j]; a[i] = [x-v*y for x,y in zip(a[i],a[j])]
        return Q([a[i][3] for i in range(3)])
    def __truediv__(self, x): return self*Q(x).inverse()

zero, one, beta = Q(), Q(1), Q((0,1,0))
p, s = beta.inverse(), one-beta
need(p*beta == one, 'beta inverse')
need(s*s-s-one-p == zero, 'quadratic coefficient relation')

def mm(a,b):
    return [[sum((x*y for x,y in zip(row,col)), zero) for col in zip(*b)] for row in a]
def mv(a,v): return tuple(row[0] for row in mm(a,[[x] for x in v]))
def add(a,b): return tuple(x+y for x,y in zip(a,b))
def sub(a,b): return tuple(x-y for x,y in zip(a,b))
def transpose(a): return [list(row) for row in zip(*a)]
def eye(): return [[one,zero],[zero,one]]
def scale(a,q): return [[x*q for x in row] for row in a]

A = [[zero,-p],[one,s]]
G = [[one,s/2],[s/2,p]]
d0, d1 = (zero,zero), (one,zero)
edges = [(0,0,d0),(1,0,d1),(0,1,d0),(2,1,d1),(0,2,d0)]
need(mm(mm(transpose(A),G),A) == scale(G,p), 'exact contraction metric')
translations = [(Q(2),Q(-3)),(Q(-1),Q(4)),(Q(5),Q(7))]
path_count = 0
gain_results = []
for gain in (F(0),F(1),F(-2),F(10),F(7,3)):
    S = [[one,Q(gain)],[zero,one]]
    Si = [[one,Q(-gain)],[zero,one]]
    B = mm(mm(S,A),Si)
    Gh = mm(mm(transpose(Si),G),Si)
    need(mm(mm(transpose(B),Gh),B) == scale(Gh,p), 'transported metric')
    need(mv(S,d1) == d1, 'digit fixed')
    transformed = [(i,j,add(mv(S,d),sub(translations[i],mv(B,translations[j]))))
                   for i,j,d in edges]
    paths = [(j,j,d0,d0) for j in range(3)]
    An, Bn = eye(), eye()
    by_depth = []
    for n in range(7):
        for i,j,v,w in paths:
            need(w == add(mv(S,v), sub(translations[i],mv(Bn,translations[j]))),
                 'affine path/Fourier endpoint gauge')
            path_count += 1
        by_depth.append(len(paths))
        if n == 6: break
        new = []
        for (i,j,d),(_,_,dp) in zip(edges,transformed):
            for source,target,v,w in paths:
                if source == j:
                    new.append((i,target,add(d,mv(A,v)),add(dp,mv(B,w))))
        paths = new
        An, Bn = mm(A,An), mm(B,Bn)
    gain_results.append({'gain':str(gain),'indexed_path_counts':by_depth})

# Test the twisted graph obstruction as an inconsistent linear system.
# Unknowns: four arbitrary common real-linear S entries, then six b coordinates.
# No commutation equations are needed: the digit equations already contradict.
twisted = [(0,0,(-one,one)),(1,0,d0),(0,1,d0),(2,1,d1),(0,2,d0)]
system = []
for (i,j,d),(ip,jp,dp) in zip(edges,twisted):
    need((i,j)==(ip,jp), 'same incidence and edge identification')
    for row in range(2):
        eq = [zero]*11
        for col in range(2):
            eq[2*row+col] = eq[2*row+col]+d[col]
            eq[4+2*j+col] = eq[4+2*j+col]-A[row][col]
        eq[4+2*i+row] = eq[4+2*i+row]+one
        eq[-1] = dp[row]
        system.append(eq)
rank = 0
for col in range(10):
    pivot = next((i for i in range(rank,len(system)) if system[i][col]),None)
    if pivot is None: continue
    system[rank],system[pivot] = system[pivot],system[rank]
    v = system[rank][col]
    system[rank] = [x/v for x in system[rank]]
    for i in range(len(system)):
        if i != rank:
            v = system[i][col]
            system[i] = [x-v*y for x,y in zip(system[i],system[rank])]
    rank += 1
inconsistent = sum(not any(row[:-1]) and bool(row[-1]) for row in system)
need(inconsistent > 0, 'twisted common-affine digit equations inconsistent')

def det3(a):
    return sum(a[0][j]*(a[1][(j+1)%3]*a[2][(j+2)%3]-a[1][(j+2)%3]*a[2][(j+1)%3]) for j in range(3))
def power(x,n):
    answer = one
    for _ in range(n): answer = answer*x
    return answer
def mult_matrix(x):
    columns = [(x*Q([int(j==i) for j in range(3)])).c for i in range(3)]
    return [[columns[j][i] for j in range(3)] for i in range(3)]
def cubic(a,b): return b*b*b-b*b*a-b*a*a-a*a*a

resolvent = []
for a,b in ((0,1),(1,2),(56,103),(1,3)):
    denominator = Q(b)-a*beta
    lam = b*beta/denominator
    need(lam/(one+F(a,b)*lam) == beta, 'Mobius inverse')
    M = mult_matrix(lam)
    norm = det3(M)
    need(norm == F(b**3,cubic(a,b)), 'resolvent norm')
    need(lam.c[0] == F(b*a*a,cubic(a,b)), 'constant coefficient is Tribonacci generating function')
    integral = all(c.denominator == 1 for c in lam.c)
    need(integral == (abs(cubic(a,b)) == 1), 'integral gain controls')
    trace = sum(M[i][i] for i in range(3))
    trace_square = sum(M[i][j]*M[j][i] for i in range(3) for j in range(3))
    char = [-norm,(trace*trace-trace_square)/2,-trace,F(1)]
    value = zero
    for c in reversed(char): value = value*lam+Q(c)
    need(value == zero, 'resolvent characteristic polynomial')
    coord_index = abs(det3([list(row) for row in zip(one.c,lam.c,(lam*lam).c)]))
    need(coord_index == F(b**6,abs(cubic(a,b))**2), 'Mobius power-basis discriminant ratio')
    if (a,b) == (0,1): expected = beta
    elif (a,b) == (1,2): expected = 2*power(beta,4)
    elif (a,b) == (56,103): expected = -103*power(beta,17)
    else: expected = lam
    need(lam == expected, 'signed Pisot power identity')
    resolvent.append({'a':a,'b':b,'norm':str(norm),'integral':integral,
       'power_basis':[str(x) for x in lam.c],
       'monic_polynomial_low_to_high':[str(x) for x in char],
       'power_basis_determinant_absolute':str(coord_index)})
need(F(9,5)**2 > 2, 'first positive branch internal contraction')
need(F(9,5)**17 > 103**2, 'negative branch internal contraction')

# All four cited norm-one pairs and their sign reversed partners are evaluated;
# global completeness comes from the cited theorem, never this finite check.
cited_pairs = [(0,-1),(1,0),(2,1),(-103,-56)]
for b,a in cited_pairs:
    need(cubic(a,b)==1 and cubic(-a,-b)==-1, 'cited Thue pair signs')
positive_b = sorted(set((b,a) for b,a in cited_pairs+ [(-b,-a) for b,a in cited_pairs] if b>0))
need(positive_b == [(1,0),(2,1),(103,56)], 'positive denominator exhaustivity reduction')

nonunit = []
for q,n in ((2,4),(103,17)):
    lam = q*power(beta,n)
    delta = q*q*power(beta.inverse(),n)
    need(delta*lam == Q(q**3), 'nonunit planar measure-balance obstruction')
    need(delta*lam/Q(q**3) == one, 'finite-place normalized measure balance')
    need(det3(mult_matrix(lam)) == q**3, 'nonunit positive signflip norm')
    nonunit.append({'q':q,'n':n,'planar_measure_multiplier':q**3,
                   'finite_place_multiplier':str(F(1,q**3))})

print(json.dumps({'status':'FINITE-EXACT independent corroboration',
    'checks':checks,'irreducibility':irreducibility,
    'gains':gain_results,'indexed_affine_paths_checked':path_count,
    'metric_identity':'B^T Q_h B = beta^(-1) Q_h',
    'twisted_linear_system':{'unknowns':10,'equations':10,
       'coefficient_rank':rank,'inconsistent_zero_rows':inconsistent},
    'resolvent_controls':resolvent,
    'nonunit_measure_controls':nonunit,
    'scope':'No symbolic dynamical conjugacy, Rauzy tiling, or pure-spectrum theorem.'},
    indent=2,sort_keys=True))
