"""Independent BFS / Euler-log audit of the formal no-three-in-line defects.

No producer imports, cyclic-run decomposition, or new geometric event census.
The geometric functional is evaluated only in its retained weights eight/nine.
Every verification uses an optimization-live explicit gate.
"""
from collections import Counter, defaultdict
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
from itertools import combinations, permutations
from math import comb, factorial
from pathlib import Path
import json
import sys

sys.stdout.reconfigure(newline='\n')
ROOT = Path(__file__).resolve().parents[1]
DEGREE = 9
GATES = 0


def need(ok, message):
    global GATES
    GATES += 1
    if not ok:
        raise RuntimeError(message)


@lru_cache(None)
def degree(m):
    return sum(token//3 for token in m)


def add(P, Q, scale=1):
    result = defaultdict(F, P)
    for m, c in Q.items():
        result[m] += scale*c
    return {m: c for m, c in result.items() if c}


def mul(P, Q):
    result = defaultdict(F)
    for a, x in P.items():
        for b, y in Q.items():
            if degree(a)+degree(b) <= DEGREE:
                result[tuple(sorted(a+b))] += x*y
    return {m: c for m, c in result.items() if c}


def scale(P, c):
    return {m: c*x for m, x in P.items() if c*x}


def part(P, weight):
    return {m: c for m, c in P.items() if degree(m) == weight}


def bfs_bank(L):
    """Enumerate edge combinations and classify graph components by BFS."""
    result = Counter()
    for k in range(min(L, DEGREE)+1):
        for selected in combinations(range(L), k):
            adjacent = [[] for _ in range(L)]
            for edge in selected:
                u, v = edge, (edge+1) % L
                adjacent[u].append(v)
                adjacent[v].append(u)
            unseen = {v for v in range(L) if adjacent[v]}
            tokens = []
            while unseen:
                start = min(unseen)
                component, todo = set(), [start]
                while todo:
                    u = todo.pop()
                    if u in component:
                        continue
                    component.add(u)
                    todo.extend(v for v in adjacent[u] if v not in component)
                unseen.difference_update(component)
                edges = sum(len(adjacent[u]) for u in component)//2
                left = sum(u % 2 for u in component)
                right = len(component)-left
                cycle = edges == len(component)
                need(cycle == (edges % 2 == 0 and left == right), "bipartite component type")
                token = 3*edges+(1 if left > right else 2 if right > left else 0)
                tokens.append(token)
            m = tuple(sorted(tokens))
            need(degree(m) == k, "BFS edge-weight grading")
            result[m] += 1
    need(sum(result.values()) == sum(comb(L, k) for k in range(min(L, DEGREE)+1)),
         "complete bounded edge-subset universe")
    return {m: F(c) for m, c in result.items()}


def euler_log(P):
    """Solve P*Euler(log P)=Euler(P) by homogeneous degree recursion."""
    need(P.get(()) == 1, "log constant normalization")
    pieces = [part(P, k) for k in range(DEGREE+1)]
    logs = [{} for _ in range(DEGREE+1)]
    for m in range(1, DEGREE+1):
        R = pieces[m].copy()
        for i in range(1, m):
            R = add(R, mul(pieces[i], logs[m-i]), -F(m-i, m))
        logs[m] = R
    result = {}
    for Pm in logs:
        result = add(result, Pm)
    return result


def falling(n, k):
    return factorial(n)//factorial(n-k) if k <= n else 0


def grid_copies(signature):
    left = sum(x[1] for x in signature)
    right = sum(x[2] for x in signature)
    multiplicities = Counter(tuple(x) for x in signature)
    automorphisms = 1
    for (edges, L, R, cycle), copies in multiplicities.items():
        local = edges if cycle else 2 if edges % 2 == 0 else 1
        automorphisms *= local**copies*factorial(copies)
    numerator = falling(8, left)*falling(8, right)
    need(numerator % automorphisms == 0, "shore-preserving automorphism quotient")
    return numerator//automorphisms


def det5(M):
    value = 0
    for perm in permutations(range(5)):
        inversions = sum(perm[i] > perm[j] for i in range(5) for j in range(i+1, 5))
        term = (-1)**inversions
        for i in range(5):
            term *= M[i][perm[i]]
        value += term
    return value


def main():
    print('source_sha256', sha256(Path(__file__).read_bytes()).hexdigest())
    banks = {L: bfs_bank(L) for L in (4, 6, 8, 10, 12, 16)}
    logs = {L: euler_log(P) for L, P in banks.items()}
    bulk = scale(logs[10], F(1, 10))
    for L in (12, 16):
        need(scale(logs[L], F(1, L)) == bulk, 'independent long-cycle bulk agreement')
    need(all(not(token % 3 == 0 and token//3 % 2 == 0) for m in bulk for token in m),
         'bulk contains only path variables')
    need(part(bulk, 1) == {(3,): F(1)}, 'Euler-log first bulk term')
    defects = {L: add(logs[L], bulk, -L) for L in (4, 6, 8)}
    for L, P in defects.items():
        need(all(degree(m) >= L for m in P), 'cycle defect begins at circumference')
    D, D5 = part(defects[4], 4), part(defects[4], 5)
    p1 = {(3,): F(1)}
    expected_D = {(12,): 1, (3, 3, 3, 3): 1, (3, 3, 7): -2, (3, 3, 8): -2,
                  (3, 9): 4, (7, 7): 1, (8, 8): 1, (13,): -2, (14,): -2}
    T = {(3, 3, 9): F(1), (3, 7, 8): F(-1), (3, 12): F(1),
         (3, 13): F(-1), (3, 14): F(-1), (7, 9): F(1), (8, 9): F(1), (15,): F(-1)}
    need(D == expected_D, 'explicit homogeneous D4')
    need(D5 == add(scale(mul(p1, D), -8), scale(T, 4)), 'D5=-8p1D4+4T')
    square = mul(D, D)
    Fcopy = add(scale(square, F(1, 2)), scale(mul(D, T), 4))
    direct_F = scale(add(add(mul(mul(banks[4], banks[4]), mul(banks[4], banks[4])),
                                   mul(banks[4], banks[12]), -4), banks[16], 3), F(1, 12))
    need(Fcopy == direct_F, 'whole quadratic polynomial vs literal skeleton contrast')
    d8, d9 = part(defects[8], 8), part(defects[8], 9)
    Ecopy = add(add(d8, d9), mul(p1, d8), 16)
    need(Ecopy == scale(add(mul(banks[8], banks[8]), banks[16], -1), F(1, 2)),
         'whole eight-cycle polynomial vs literal skeleton contrast')
    for P in (D, D5, d8, d9):
        need(sum(P.values()) == 0, 'edge-count specialization kills defects')
    path = ROOT/'05-knowledge/results/overnight2_20260906_no3line_third_certificates.json'
    need(sha256(path.read_bytes()).hexdigest() ==
         'f2c566ac1b2bcedb530af72f3a290479841e7db87844b331558bed68c93ba727', 'frozen geometry identity')
    data = json.loads(path.read_text())
    geometry = {}
    for row in data['profiles']:
        sig = row['signature']
        m = tuple(sorted(3*e+(1 if l > r else 2 if r > l else 0) for e, l, r, cyc in sig))
        need(m not in geometry and degree(m) in (8, 9), 'retained profile uniqueness and scope')
        copies = grid_copies(sig)
        need(copies == row['grid_copies'], 'independent injective grid-copy normalization')
        geometry[m] = F(row['unordered_event_triples'], copies)
        need(Ecopy.get(m, 0) == F(row['profile_coefficients_A_B_D_E_F'][3]), 'eight coefficient per profile')
        need(Fcopy.get(m, 0) == F(row['profile_coefficients_A_B_D_E_F'][4]), 'quadratic coefficient per profile')
        need([6*geometry[m]*Ecopy.get(m, 0), 6*geometry[m]*Fcopy.get(m, 0)] ==
             list(map(F, row['ordered_weighted_contributions_E_F'])), 'ordered whole-event normalization')
    need(len(geometry) == 150, 'complete inherited retained bank')
    def W(P):
        need(all(degree(m) in (8, 9) for m in P), 'functional evaluated only in retained weights')
        return sum((c*geometry.get(m, F(0)) for m, c in P.items()), F(0))
    E, FF = 6*W(Ecopy), 6*W(Fcopy)
    need((E, FF) == (F(172483, 529200), F(11881, 50400)), 'inherited full third-moment coefficients')
    need((3*W(square), 24*W(mul(D, T))) == (F(456371, 2116800), F(42631, 2116800)),
         'signed square and cross contribution split')
    Q = {(3, 3, 7): F(1), (7, 7): F(-1)}
    need(all(degree(m) == 4 for m in Q), 'homogeneous weight-four hostile')
    Q2 = mul(Q, Q)
    need(all(m in geometry and degree(m) == 8 and all(t % 3 for t in m if t//3 % 2 == 0)
             for m in Q2), 'square supported on realized weight-eight forests')
    need(W(Q2) == F(-397, 529200), 'restricted homogeneous square is negative')
    need([geometry[m] for m in ((3, 3, 3, 3, 7, 7), (3, 3, 7, 7, 7), (7, 7, 7, 7))] ==
         [F(16747, 1058400), F(2749, 235200), F(1, 147)], 'three forest weights')
    basis = [[1, 0, 0, 0, 0], [1, 1, 0, 0, 1], [1, 0, 1, 0, 0],
             [1, 0, 0, 2, 0], [1, 4, 0, 0, 16]]
    need(det5(basis) == 24, 'valid n8 cycle basis full rank')
    multisets = ((0, 0, 2), (2, 0, 1), (4, 0, 0))
    need(all(multisets[0][j]-2*multisets[1][j]+multisets[2][j] == 0 for j in range(3)),
         'direct contrast kills every component-additive function')
    need(8*FF == F(11881, 6300), 'third-cumulant nonadditivity contrast')
    print('BFS_CYCLE_LENGTHS', sorted(banks))
    print('LOG_METHOD Euler homogeneous derivative recurrence; bulk10=bulk12=bulk16 through9')
    print('DEFECT_COUNTS D4,D5,D8,D9', len(D), len(D5), len(d8), len(d9))
    print('WHOLE_POLYNOMIAL_COUNTS eight,quadratic', len(Ecopy), len(Fcopy))
    print('INHERITED_PROFILES', len(geometry), 'ORDERED_M3_COEFFICIENTS', str(E), str(FF))
    print('SQUARE_CROSS_SPLIT', str(3*W(square)), str(24*W(mul(D, T))))
    print('HOMOGENEOUS_WEIGHT4_FOREST_SQUARE', str(W(Q2)))
    print('THIRD_CUMULANT_COMPONENT_CONTRAST', str(8*FF), 'BASIS_MINOR', det5(basis))
    print('GEOMETRY_SHA256', sha256(path.read_bytes()).hexdigest())
    print('active_gates', GATES)
    print('PASS independent graph enumeration, formal identity controls, signed square and cumulant obstruction')


if __name__ == '__main__':
    main()
