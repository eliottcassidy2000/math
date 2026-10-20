#!/usr/bin/env python3
"""Bounded, independent rational controls; never import or execute the packet verifier.
Run from any directory. Reads the packet, receipts, and inherited q17 certificate.
Prints a JSON receipt; no writes are performed by this program.
"""
from fractions import Fraction as F
from itertools import product, combinations
from collections import Counter
from hashlib import sha256
from pathlib import Path
from math import prod
import json

ROOT = Path(__file__).resolve().parent.parent
PACK = ROOT
checks = Counter()
pins = {}
def need(group, condition):
    if not condition:
        raise RuntimeError('Independent control failed: ' + group)
    checks[group] += 1

def digest(path):
    return sha256(path.read_bytes()).hexdigest()

def rank(rows):
    a = [[F(x) for x in row] for row in rows]
    r = 0
    for c in range(len(a[0])):
        pivot = next((k for k in range(r, len(a)) if a[k][c]), None)
        if pivot is None:
            continue
        a[r], a[pivot] = a[pivot], a[r]
        s = a[r][c]
        a[r] = [x / s for x in a[r]]
        for k in range(len(a)):
            if k != r:
                s = a[k][c]
                a[k] = [x - s*y for x, y in zip(a[k], a[r])]
        r += 1
        if r == len(a):
            break
    return r

def mm(a, b):
    return [[sum(x*y for x,y in zip(row,col)) for col in zip(*b)] for row in a]
def tr(a): return list(map(list, zip(*a)))
def sq(a): return [[x*x for x in row] for row in a]
def ident(n): return [[F(i == j) for j in range(n)] for i in range(n)]

# Receipt integrity, without rerunning the supplied 3,708-check program.
for line in (PACK / 'SHA256SUMS').read_text().splitlines():
    expected, rel = line.split()
    pins[rel] = digest(ROOT / rel)
    need('supplied_hash_matches', pins[rel] == expected)
a = (PACK / 'verification.json').read_bytes()
b = (PACK / 'verification_optimized.json').read_bytes()
receipt = json.loads(a)
need('supplied_receipt_integrity', a == b)
need('supplied_receipt_integrity', receipt['status'] == 'PASS')
need('supplied_receipt_integrity', sum(receipt['checks'].values()) == receipt['exact_check_total'] == 3708)

# Explicit radius-two derivative identity at a rational-square probability law.
sqrt_p = [F(1,3), F(2,3), F(2,3)]
p = [s*s for s in sqrt_p]
u = [F(1), F(-2), F(1)]
need('fisher_normalization', sum(p) == 1 and sum((2*s)**2 for s in sqrt_p) == 4)
need('fisher_normalization', sum(u) == 0 and sum(x*x/q for x,q in zip(u,p)) == sum((x/s)**2 for x,s in zip(u,sqrt_p)))

# Enumerate actual signs, rather than packet bitmask characters.
words = list(product((-1,1), repeat=4))
constant = [1] * len(words)
pair_rows = [[x[i]*x[j] for x in words] for i,j in combinations(range(4),2)]
quartic = [prod(x) for x in words]
star = [min(x.count(-1), x.count(1)) for x in words]
O = [constant] + pair_rows
need('observer_rank', rank(O) == 7)
need('observer_rank', rank(O + [star]) - rank(O) == 1)
need('observer_rank', rank(O + [quartic]) == 8)
# Arbitrary expectation versus individual-character restrictions are genuinely different.
f = [x[0]*x[1] + x[2]*x[3] for x in words]
need('observer_rank', rank([constant,f]) == 2 and rank([constant,pair_rows[0],pair_rows[-1]]) == 3)
plus = [F(1+v,16) for v in quartic]
minus = [F(1-v,16) for v in quartic]
need('maximal_invisible_fiber', sum(plus) == sum(minus) == 1 and all(a*b == 0 for a,b in zip(plus,minus)))
need('maximal_invisible_fiber', all(sum(w*v for w,v in zip(plus,row)) == sum(w*v for w,v in zip(minus,row)) == 0 for row in pair_rows))
need('maximal_invisible_fiber', sum(w*v for w,v in zip(plus,star)) == F(3,2) and sum(w*v for w,v in zip(minus,star)) == 1)
for t in (F(-3,5), F(0), F(5,13)):
    need('fisher_parity_speed', sum(F(1,256)/(F(1,16)*(1+t*q)) for q in quartic) == 1/(1-t*t))

# Klein-four transform, and actual negative probability for the Gram hostile.
reps3 = [(1,) + x for x in product((-1,1), repeat=2)]
H = [[1, x[0]*x[1], x[0]*x[2], x[1]*x[2]] for x in reps3]
need('triangle_character_transform', mm(H,tr(H)) == [[4*int(i==j) for j in range(4)] for i in range(4)])
for a,b,c in [(F(-1,2),)*3, (F(-1,3),)*3, (F(1),)*3, (F(0),)*3, (F(1,5),F(-2,5),F(3,5))]:
    mu = [1,a,b,c]
    masses = [sum(x*y for x,y in zip(row,mu))/4 for row in H]
    M = [[1,a,b,c],[a,1,c,b],[b,c,1,a],[c,b,a,1]]
    need('triangle_character_transform', mm(mm(tr(H), [[masses[i]*int(i==j) for j in range(4)] for i in range(4)]),H) == M)
need('triangle_hostile', 1 + 2*F(-1,2)**3 - 3*F(-1,2)**2 == 0 and (1+3*F(-1,2))/4 == F(-1,8))
need('triangle_q17_local', [1+3*F(-1,3)] + [1-F(-1,3)]*3 == [0,F(4,3),F(4,3),F(4,3)])

# All scalar moments can occur in a proper principal test.
labels = [frozenset()] + [frozenset(x) for x in combinations(range(4),2)] + [frozenset(range(4))]
proper_products = {a ^ b for a in labels[:-1] for b in labels[:-1]}
need('proper_row_all_scalars', proper_products == set(labels))
for m in (3,4,5,6):
    n = 2**(m-1)
    c = F(-1,n-2)
    need('missing_row_spectrum', all(1-c > 0 and 1+(r-1)*c >= 0 for r in range(1,n)))
    need('missing_row_spectrum', 1+(n-1)*c == F(-1,n-2))
    need('missing_row_masses', (1+(n-1)*c)/n == F(-1,n*(n-2)) and (1-c)/n == F(n-1,n*(n-2)))
    eps = F(1,n*(n-2))
    # At threshold the (N-1)-row all-ones eigenvalue is zero; above it it is negative.
    def proper_eigen(e):
        moment = -(1+n*e)/F(n-1)
        return 1+(n-2)*moment
    need('missing_row_sharp_threshold', proper_eigen(eps) == 0 and proper_eigen(eps+F(1,n*n)) < 0)
balanced = [F(1,6) if x.count(-1) == 2 else F(0) for x in words]
mix = [(F(1,16)+w)/2 for w in balanced]
need('pair_vs_full_packet', sum(mix) == 1 and all(sum(w*v for w,v in zip(mix,row)) == F(-1,6) for row in pair_rows))
need('pair_vs_full_packet', sum(w*v for w,v in zip(mix,quartic)) == F(1,2))

# Separate all-sign enumeration for the frozen packet's added parity boundary.
for m in range(2,8):
    balanced_words = [x for x in product((-1,1),repeat=m) if abs(sum(x)) == m%2]
    endpoint = F(-1,m if m%2 else m-1)
    need('equicorrelation_boundary', all(F(sum(x[i]*x[j] for x in balanced_words),len(balanced_words)) == endpoint for i,j in combinations(range(m),2)))
    need('equicorrelation_boundary', F(sum(sum(x)**2 for x in balanced_words),len(balanced_words)) == m%2)
    need('equicorrelation_boundary', endpoint + F(1,m-1) == (F(1,m*(m-1)) if m%2 else 0))

# Rebuild only the small graph, not either refutation tree.
vertices = list(product(range(17),repeat=2))
adj = [set() for _ in vertices]
for i,(x,y) in enumerate(vertices):
    for j in range(i):
        u,v = vertices[j]
        if ((x-u)**2 - 3*(y-v)**2) % 17 == 1:
            adj[i].add(j); adj[j].add(i)
edges = [(i,j) for i in range(289) for j in adj[i] if i<j]
triangles = {(i,j,k) for i,j in edges for k in adj[i] & adj[j] if j<k}
need('q17_graph_rebuild', len(edges) == 2601 and len(triangles) == 1734 and all(len(s) == 18 for s in adj))
need('q17_graph_rebuild', all(len(adj[i] & adj[j]) == 2 for i,j in edges))
qbase = ROOT.parent / 'plane_binary_transfer_20261020'
if not qbase.exists():
    qbase = ROOT.parent / 'math/05-knowledge/results/plane_binary_transfer_20261020'
cert = qbase / 'addendum-four-defects/results/q17_no_two_defects_certificate.json.gz'
pins[str(cert.relative_to(ROOT.parent))] = digest(cert)
need('q17_dependency_identity_only', pins[str(cert.relative_to(ROOT.parent))] == '251cf994140ed97daf9fc963dab8544c3c6018c51cad8f39f095edd3087e23fe')
need('q17_distance_arithmetic', 2601-2*1730 - 2601*F(-1,3) == 8 and 51**2 == 2601)

# Independent rational orthogonal example: compose a reflection with itself.
U = [[F(3,5), F(4,5)],[F(4,5), F(-3,5)]]
V = [[F(5,13), F(12,13)], [F(-12,13), F(5,13)]]
need('walsh_composition', mm(U,tr(U)) == ident(2) and mm(V,tr(V)) == ident(2))
coherent, classical = sq(mm(V,U)), mm(sq(V),sq(U))
D = [[coherent[i][j]-classical[i][j] for j in range(2)] for i in range(2)]
need('walsh_composition', any(x for row in D for x in row))
need('walsh_composition', all(D[i][j] == 2*V[i][0]*V[i][1]*U[0][j]*U[1][j] for i,j in product(range(2),repeat=2)))
need('walsh_composition', all(sum(row) == 0 for row in D) and all(sum(col) == 0 for col in zip(*D)))
# Order-four Walsh needs no irrational entries.
W = [[F(1,2),F(1,2),F(1,2),F(1,2)], [F(1,2),F(-1,2),F(1,2),F(-1,2)], [F(1,2),F(1,2),F(-1,2),F(-1,2)], [F(1,2),F(-1,2),F(-1,2),F(1,2)]]
need('walsh_involution_gap', mm(W,W) == ident(4) and sq(W) == [[F(1,4)]*4 for _ in range(4)])
need('walsh_involution_gap', mm(sq(W),sq(W)) == sq(W) and sum(abs(F(i==0)-F(1,4)) for i in range(4))/2 == F(3,4))

# Cube order duality versus polar combinatorics are distinct finite controls.
grid = list(product((F(0),F(1,2),F(1)), repeat=3))
le = lambda x,y: all(a<=b for a,b in zip(x,y))
complement = lambda x: tuple(1-a for a in x)
need('cube_order_duality', all(complement(complement(x)) == x for x in grid))
need('cube_order_duality', all(le(x,y) == le(complement(y),complement(x)) for x in grid for y in grid))
need('cube_polar_vertex_count', 2**3 != 2*3)

print(json.dumps({'status':'PASS','exact_control_total':sum(checks.values()),'checks':dict(sorted(checks.items())), 'reviewed_hashes':pins, 'scope':'Independent bounded exact controls. Supplied 3708-check verifier and inherited q17 proof trees were not rerun. Q17 graph rebuilt; certificate hash identified only. Proofs and type distinctions separately reviewed.'}, indent=2,sort_keys=True))
