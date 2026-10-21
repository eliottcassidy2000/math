# AUDIT a20 (independent auditor, 2026-10-21): HYP-9357 design-theoretic claims.
#  (1) Paley designs D(q) = translates of the nonzero squares of F_q (q = 3 mod 4) for q in {7,11,19,23,27,31,43,47}:
#      |Aut(D(q))| (points/blocks kept separate) and 2-transitivity (point stabiliser transitive on the other points),
#      by nauty's dreadnaut on the point-block incidence graph;
#  (2) Singer difference sets of PG(m-1,2) in Z_{2^m-1} (m = 3, 5, 7): full multiplier groups vs those of QR_q, and the
#      block-triple intersection distribution (a design invariant) of D(31) vs PG(4,2) and D(127) vs PG(6,2);
#  (3) the Fano plane: |Aut| = 168 and the point stabiliser has the element-order profile of S_4.
import itertools, subprocess, re
from collections import Counter

def field(q):
    """return (elements, add, mul) for F_q, q prime or 27"""
    if q == 27:
        els = [(a, b, c) for a in range(3) for b in range(3) for c in range(3)]
        def add(x, y): return tuple((u + v) % 3 for u, v in zip(x, y))
        def mul(x, y):
            # polynomials a + b t + c t^2 modulo t^3 = t + 1
            r = [0] * 5
            for i in range(3):
                for j in range(3):
                    r[i + j] += x[i] * y[j]
            for d in (4, 3):
                c = r[d]; r[d] = 0
                r[d - 3] += c; r[d - 2] += c   # t^d = t^(d-3) * (t + 1)
            return tuple(v % 3 for v in r[:3])
        return els, add, mul
    els = list(range(q))
    return els, (lambda x, y: (x + y) % q), (lambda x, y: x * y % q)

def paley_blocks(q):
    els, add, mul = field(q)
    zero = els[0]
    sq = {mul(x, x) for x in els if x != zero}
    idx = {x: i for i, x in enumerate(els)}
    return [frozenset(idx[add(s, b)] for s in sq) for b in els]

def singer_blocks(m):
    q = 2 ** m - 1
    poly = {3: 0b1011, 5: 0b100101, 7: 0b10000011}[m]   # x^3+x+1, x^5+x^2+1, x^7+x+1 (primitive)
    powers = []; x = 1
    for i in range(q):
        powers.append(x); x <<= 1
        if x >> m & 1: x ^= poly
    def tr(v):  # absolute trace F_{2^m} -> F_2: sum of conjugates
        t = 0; y = v
        for _ in range(m):
            t ^= y
            # square y
            r = 0
            for i in range(m):
                if y >> i & 1:
                    r ^= powers[(2 * i) % q] if True else 0
            y = r
        return t
    D = [i for i in range(q) if tr(powers[i]) == 0]
    return D, [frozenset((d + b) % q for d in D) for b in range(q)]

def dread(q, blocks, stab_point=None):
    n = 2 * q
    edges = []
    for b, B in enumerate(blocks):
        for p in B:
            edges.append((p, q + b))
    adj = {v: [] for v in range(n)}
    for a, b in edges:
        adj[a].append(b)
    g = "; ".join(f"{v}:" + " ".join(map(str, adj[v])) for v in range(n) if adj[v])
    if stab_point is None:
        part = f"f=[0:{q-1}|{q}:{n-1}]"
    else:
        part = f"f=[{stab_point}|" + ",".join(str(v) for v in range(q) if v != stab_point) + f"|{q}:{n-1}]"
    inp = f"n={n} g {g}. {part} x o q\n"
    out = subprocess.run(["dreadnaut"], input=inp, capture_output=True, text=True).stdout
    grp = re.search(r"grpsize=([0-9.e+]+)", out).group(1)
    orb_text = out.split("cpu time")[-1]
    return grp, orb_text

def orbit_of(orb_text, v, q):
    # parse 'a:b (k); c d e (k);' style orbit listing
    text = orb_text.replace("\n", " ")
    for chunk in text.split(";"):
        chunk = re.sub(r"\(\d+\)", "", chunk).strip()
        if not chunk:
            continue
        members = set()
        for tok in chunk.split():
            if ":" in tok:
                a, b = map(int, tok.split(":")); members |= set(range(a, b + 1))
            elif tok.isdigit():
                members.add(int(tok))
        if v in members:
            return members
    return None

print("(1) Paley designs")
for q in [7, 11, 19, 23, 27, 31, 43, 47]:
    B = paley_blocks(q)
    grp, _ = dread(q, B)
    grp0, orb0 = dread(q, B, stab_point=0)
    o1 = orbit_of(orb0, 1, q)
    pts_orbit = {v for v in (o1 or set()) if v < q}
    print(f"  q={q}: |Aut(D(q))| = {grp}; stabiliser of point 0 has the other points in one orbit "
          f"(2-transitive): {pts_orbit == set(range(1, q))}")

print("(2) Singer vs QR multiplier groups and triple-intersection invariants")
def multipliers(D, q):
    D = set(D); res = []
    for a in range(1, q):
        aD = {a * d % q for d in D}
        if any(aD == {(d + g) % q for d in D} for g in range(q)):
            res.append(a)
    return res
def triple_profile(blocks):
    return Counter(len(a & b & c) for a, b, c in itertools.combinations(blocks, 3))
for m in [3, 5, 7]:
    q = 2 ** m - 1
    D, SB = singer_blocks(m)
    QR = sorted({x * x % q for x in range(1, q)})
    # check Singer set is a difference set with lambda = (q-3)/4
    diffs = Counter((a - b) % q for a in D for b in D if a != b)
    assert set(diffs.values()) == {(q - 3) // 4}, "not a difference set"
    MS = multipliers(D, q); MQ = multipliers(QR, q)
    print(f"  q={q}: Singer set size {len(D)}; |M(Singer)| = {len(MS)} {MS if len(MS) < 10 else ''}; |M(QR)| = {len(MQ)}")
    if q <= 31:
        print(f"     triple-intersection profile Singer {sorted(triple_profile(SB).items())}")
        print(f"     triple-intersection profile Paley  {sorted(triple_profile(paley_blocks(q)).items())}")
    else:
        # q = 127: compare profiles on all block triples through block 0 (still a design invariant up to transitivity)
        def prof0(bl):
            return Counter(len(bl[0] & b & c) for b, c in itertools.combinations(bl[1:], 2))
        print(f"     triples through block 0: Singer {sorted(prof0(SB).items())}")
        print(f"     triples through block 0: Paley  {sorted(prof0(paley_blocks(q)).items())}")

print("(3) Fano plane automorphisms")
fano = [frozenset((d + b) % 7 for d in (1, 2, 4)) for b in range(7)]
fs = set(fano)
auts = [p for p in itertools.permutations(range(7)) if all(frozenset(p[x] for x in L) in fs for L in fano)]
stab = [p for p in auts if p[0] == 0]
def order(p):
    k = 1; x = list(p)
    while x != list(range(7)):
        x = [p[i] for i in x]; k += 1
    return k
print(f"  |Aut(Fano)| = {len(auts)}; |Stab(0)| = {len(stab)}; element orders in Stab(0): {sorted(Counter(order(p) for p in stab).items())} "
      f"(S_4 has [(1,1),(2,9),(3,8),(4,6)])")
