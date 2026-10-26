# Stress test of the legal-path reconstruction beyond the exhaustive range (n = 9..13), by bitmask dynamic
# programming instead of path enumeration.  Random tournaments are almost always 2-strong (no missed arcs, nothing
# to test), so we build tournaments with cut vertices on purpose: random blocks C_1 => ... => C_k, a cut vertex u
# with random arcs to the blocks (forced u -> C_1 and C_k -> u, so T is strong), optionally nested twice (the whole
# thing becomes a block, or the cut-vertex construction is applied again with it as a block), then a random
# relabelling.  For each tournament:
#   (R)  decode from path data alone (used arcs, cut positions, block endpoint pairs) == T;
#   (S)  every missed arc u -> v inside a strong T: T - u non-strong with v terminal, T - v non-strong with u initial;
#   (F)  for random nonempty sets D of missed arcs, T_D has a Hamilton path through a reversed arc of D.
import sys, random
sys.path.insert(0, '../joint_port_bisimulation_20261026')
from rooted import is_strong, relabel
from ham_recovery import components

def bits(x):
    while x:
        b = x & -x; yield b.bit_length() - 1; x ^= b

def path_data(adj, mask):
    """F[S]: possible last vertices of a path covering exactly S (S within mask); B[S]: possible first vertices."""
    n = len(adj)
    inn = [sum(1 << u for u in range(n) if (adj[u] >> v) & 1) for v in range(n)]
    F = {}; B = {}
    subs = [S for S in range(1 << n) if S & ~mask == 0 and S]
    subs.sort(key=lambda S: bin(S).count('1'))
    for v in bits(mask): F[1 << v] = 1 << v; B[1 << v] = 1 << v
    for S in subs:
        f = F.get(S, 0)
        for v in bits(f):
            for w in bits(adj[v] & mask & ~S): F[S | 1 << w] = F.get(S | 1 << w, 0) | 1 << w
        b = B.get(S, 0)
        for v in bits(b):
            for w in bits(inn[v] & mask & ~S): B[S | 1 << w] = B.get(S | 1 << w, 0) | 1 << w
    return F, B, subs

def legal_data(adj):
    n = len(adj); full = (1 << n) - 1
    F, B, subs = path_data(adj, full)
    used = set(); prefixes = {}
    for S in subs:
        if S == full: continue
        rest = B.get(full ^ S, 0)
        hit = False
        for a in bits(F.get(S, 0)):
            for b in bits(adj[a] & rest): used.add((a, b)); hit = True
        if hit: prefixes.setdefault(bin(S).count('1'), set()).add(S)
    return used, prefixes

def block_ends(adj, C):
    """(x, y) such that T[C] has a Hamilton path from x to y"""
    n = len(adj); ends = set()
    for x in bits(C):
        R = {1 << x: 1 << x}
        order = sorted((S for S in range(1 << n) if S & ~C == 0 and (S >> x) & 1), key=lambda S: bin(S).count('1'))
        for S in order:
            for v in bits(R.get(S, 0)):
                for w in bits(adj[v] & C & ~S): R[S | 1 << w] = R.get(S | 1 << w, 0) | 1 << w
        for y in bits(R.get(C, 0)):
            if y != x: ends.add((x, y))
    return ends

def decode(adj):
    n = len(adj)
    used, prefixes = legal_data(adj)
    cuts = [m for m in range(1, n) if len(prefixes.get(m, ())) == 1]
    chain = [0] + [next(iter(prefixes[m])) for m in cuts] + [(1 << n) - 1]
    blocks = [chain[i + 1] & ~chain[i] for i in range(len(chain) - 1)]
    bi = {v: i for i, C in enumerate(blocks) for v in bits(C)}
    ends = {i: block_ends(adj, C) if bin(C).count('1') >= 3 else set() for i, C in enumerate(blocks)}
    rows = [0] * n
    for x in range(n):
        for y in range(n):
            if x == y: continue
            if (x, y) in used: rows[x] |= 1 << y; continue
            if (y, x) in used: continue
            if bi[x] != bi[y]:
                if bi[x] < bi[y]: rows[x] |= 1 << y
                continue
            fx, fy = (x, y) in ends[bi[x]], (y, x) in ends[bi[x]]
            if fx == fy: return None
            if fx: rows[x] |= 1 << y
    return tuple(rows), used

def random_block(rows, verts, rng):
    for i, a in enumerate(verts):
        for b in verts[i + 1:]:
            if rng.random() < 0.5: rows[a] |= 1 << b
            else: rows[b] |= 1 << a

def build(rng, n_target, depth=0):
    """blocks => chain, plus a cut vertex u; blocks are random, or (recursively) earlier builds.
    Half the time u has a sparse in-neighbourhood (a few vertices of the terminal block, never the penultimate
    block), which is what produces arcs missed by every Hamilton path."""
    sizes = []
    while sum(abs(s) for s in sizes) < n_target - 1:
        s = rng.choice([1, 1, 3, 3, 4, 5, -5, -6]) if depth == 0 else rng.choice([1, 3, 3, 4])
        sizes.append(s)
    if len(sizes) < 2: sizes.append(1)
    n = sum(abs(s) for s in sizes) + 1; rows = [0] * n; blocks = []; start = 0
    for s in sizes:
        verts = list(range(start, start + abs(s))); blocks.append(verts); start += abs(s)
        if s < 0:                              # recursive block: a strong build with its own cut vertex
            inner = build(rng, abs(s), depth + 1)
            assert len(inner) >= abs(s)
            inner = inner[:abs(s)] if len(inner) == abs(s) else None
            if inner is None:
                random_block(rows, verts, rng)
            else:
                for i, a in enumerate(verts):
                    for j, b in enumerate(verts):
                        if (inner[i] >> j) & 1: rows[a] |= 1 << b
        else:
            random_block(rows, verts, rng)
    for i, Bk in enumerate(blocks):
        for C in blocks[i + 1:]:
            for a in Bk:
                for b in C: rows[a] |= 1 << b
    u = n - 1
    if rng.random() < 0.5:
        for x in range(n - 1): rows[u] |= 1 << x
        k = len(blocks)
        ins = rng.sample(blocks[-1], rng.randint(1, min(2, len(blocks[-1]))))
        ins += [x for i in range(k - 2) for x in blocks[i] if rng.random() < 0.3 and i > 0]
        for x in ins: rows[u] &= ~(1 << x); rows[x] |= 1 << u
    else:
        for x in range(n - 1):
            if rng.random() < 0.5: rows[u] |= 1 << x
            else: rows[x] |= 1 << u
        a = rng.choice(blocks[0]); rows[a] &= ~(1 << u); rows[u] |= 1 << a
        b = rng.choice(blocks[-1]); rows[u] &= ~(1 << b); rows[b] |= 1 << u
    return rows

def sub(adj, keep):
    idx = {v: i for i, v in enumerate(keep)}
    return tuple(sum(1 << idx[w] for w in keep if (adj[v] >> w) & 1) for v in keep)

def comps_minus(adj, x):
    n = len(adj); keep = [y for y in range(n) if y != x]
    return [[keep[i] for i in c] for c in components(sub(adj, keep))]

rng = random.Random(20261026)
trials = int(sys.argv[1]) if len(sys.argv) > 1 else 200
tot = dict(cases=0, decoded=0, strong=0, with_missed=0, missed=0, S_bad=0, flipsets=0, F_bad=0, nmax=0)
for t in range(trials):
    rows = build(rng, rng.randint(8, 11))
    if rng.random() < 0.3:                    # nest: append 1-2 vertices dominated by the whole build
        extra = rng.randint(1, 2); n0 = len(rows)
        rows = rows + [0] * extra
        for e in range(n0, n0 + extra):
            for a in range(n0): rows[a] |= 1 << e
        for e in range(n0, n0 + extra):
            for f in range(e + 1, n0 + extra): rows[e] |= 1 << f
    n = len(rows); perm = list(range(n)); rng.shuffle(perm)
    adj = relabel(tuple(rows), perm)
    tot['nmax'] = max(tot['nmax'], n)
    res = decode(adj)
    tot['cases'] += 1
    if res is None: continue
    dec, used = res
    tot['decoded'] += dec == adj
    missed = [(u, v) for u in range(n) for v in range(n) if (adj[u] >> v) & 1 and (u, v) not in used]
    tot['missed'] += len(missed); tot['with_missed'] += bool(missed)
    if is_strong(adj):
        tot['strong'] += 1
        for u, v in missed:
            cu, cv = comps_minus(adj, u), comps_minus(adj, v)
            tot['S_bad'] += not (len(cu) >= 2 and v in cu[-1] and len(cv) >= 2 and u in cv[0])
    for _ in range(3 if missed else 0):
        D = rng.sample(missed, rng.randint(1, len(missed)))
        new = list(adj)
        for u, v in D: new[u] &= ~(1 << v); new[v] |= 1 << u
        used2, _ = legal_data(tuple(new))
        tot['flipsets'] += 1
        tot['F_bad'] += not any((v, u) in used2 for u, v in D)
print(tot)
