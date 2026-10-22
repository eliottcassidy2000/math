# AUDIT (task 4, Section 1 numbers): the excess-1,2,3 core lists cores_e1/2/3.json.
#  * own validity checks: connected, min degree >= 3 (loop counts 2), E = V + excess, pairwise non-isomorphic (own
#    canonical form: min over all vertex permutations of the sorted edge multiset);
#  * own automorphism count |Aut| = |vertex automorphisms| * prod_(a<b) m_ab! * prod_a m_aa! 2^(m_aa)  vs the stored 'aut';
#  * counts by V, sum 1/|G| (claim 1661/240), connected cubic part (claim 1105/1152);
#  * configuration-model completeness identity: [x^V y^E] exp(sum_(cores, excess<=3) x^V y^E/|Aut|) =
#    (1/V!) sum over ordered degree sequences (all >= 3, sum 2E) of (2E-1)!!/prod d_v!, for all 1 <= E - V <= 3.
import json, itertools
from fractions import Fraction
from math import factorial
from collections import Counter, defaultdict
HERE = '/Users/e/Documents/GitHub/math/.claude/worktrees/gifted-chatterjee-298788/04-computation/experiments/sierpinski_ocf_20261022'

def mult(V, edges):
    m = Counter()
    for a, b in edges: m[(min(a, b), max(a, b))] += 1
    return m

def own_aut(V, edges):
    m = mult(V, edges)
    va = 0
    for p in itertools.permutations(range(V)):
        if all(m.get((min(p[a], p[b]), max(p[a], p[b])), 0) == c for (a, b), c in m.items()): va += 1
    f = 1
    for (a, b), c in m.items():
        f *= factorial(c) * (2 ** c if a == b else 1)
    return va * f

def own_canon(V, edges):
    return min(tuple(sorted((min(p[a], p[b]), max(p[a], p[b])) for a, b in edges)) for p in itertools.permutations(range(V)))

def connected(V, edges):
    adj = defaultdict(set)
    for a, b in edges: adj[a].add(b); adj[b].add(a)
    seen = {0}; st = [0]
    while st:
        x = st.pop()
        for y in adj[x]:
            if y not in seen: seen.add(y); st.append(y)
    return len(seen) == V

gf = defaultdict(Fraction)       # (V, E) -> sum 1/|Aut| over connected cores
for e in (1, 2, 3):
    cs = json.load(open(f'{HERE}/cores_e{e}.json'))
    byV = Counter(c['V'] for c in cs)
    ok = True; canons = set(); s = Fraction(0); cub = Fraction(0)
    for c in cs:
        V, edges = c['V'], [tuple(x) for x in c['edges']]
        deg = Counter()
        for a, b in edges: deg[a] += 1; deg[b] += 1
        ok &= connected(V, edges) and len(edges) == V + e and all(deg[v] >= 3 for v in range(V))
        canons.add(own_canon(V, edges))
        a = own_aut(V, edges)
        ok &= (a == c['aut'])
        s += Fraction(1, a); gf[(V, len(edges))] += Fraction(1, a)
        if all(deg[v] == 3 for v in range(V)): cub += Fraction(1, a)
    print(f"excess {e}: {len(cs)} cores, by V {dict(sorted(byV.items()))}; valid + own |Aut| == stored aut: {ok}; "
          f"pairwise non-isomorphic: {len(canons) == len(cs)}; sum 1/|Aut| = {s}; connected cubic part = {cub}")

# exponential formula: all (not necessarily connected) min-degree->=3 multigraphs with total excess <= 3
def exp_gf(gf, emax=3):
    # generating function in x (vertices) and excess k = E - V; coefficients Fraction
    C = defaultdict(Fraction)
    for (V, E), w in gf.items(): C[(V, E - V)] += w
    # exp via sum over multisets: use exp(f) = sum f^k/k!, truncated at excess <= emax (each component has excess >= 1)
    res = defaultdict(Fraction); res[(0, 0)] = Fraction(1)
    term = {(0, 0): Fraction(1)}
    for k in range(1, emax + 1):
        new = defaultdict(Fraction)
        for (v1, x1), w1 in term.items():
            for (v2, x2), w2 in C.items():
                if x1 + x2 <= emax: new[(v1 + v2, x1 + x2)] += w1 * w2
        term = {key: w / k for key, w in new.items()}
        for key, w in term.items(): res[key] += w
    return res

def config_count(V, E):
    tot = Fraction(0)
    def rec(i, rem, acc):
        nonlocal tot
        if i == V:
            if rem == 0: tot += acc
            return
        for d in range(3, rem - 3 * (V - i - 1) + 1):
            rec(i + 1, rem - d, acc / factorial(d))
    rec(0, 2 * E, Fraction(1))
    dfact = 1
    for k in range(1, 2 * E, 2): dfact *= k
    return tot * dfact / factorial(V)

ex = exp_gf(gf)
allok = True
for V in range(1, 7):
    for x in range(1, 4):
        E = V + x
        lhs = ex.get((V, x), Fraction(0)); rhs = config_count(V, E)
        allok &= (lhs == rhs)
        print(f"   (V, E) = ({V}, {E}): exp(cores) = {lhs}, configuration model = {rhs}, equal: {lhs == rhs}")
print("configuration-model identity at every (V, E) with 1 <= E - V <= 3:", allok)
# JKLP check: all cubic multigraphs on 2r vertices: e_r = (6r-1)!!/((2r)! 6^(2r)); connected part at r = 3
def df(m):
    r = 1
    for k in range(1, m + 1, 2): r *= k
    return r
e1, e2, e3 = [Fraction(df(6 * r - 1), factorial(2 * r) * 6 ** (2 * r)) for r in (1, 2, 3)]
e2c = e2 - e1 ** 2 / 2
e3c = e3 - e1 * e2c - e1 ** 3 / 6
print(f"JKLP e_1, e_2, e_3 = {e1}, {e2}, {e3}; connected cubic excess-3 part = {e3c}")
