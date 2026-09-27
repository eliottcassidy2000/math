"""procgen_cnpos_20260926_core.py -- first-return reduction of C_n to a residual problem, lift, verifier,
and two solvers for the residual problem.

Session collatz-procgen-20260922, lane "cnpos" (2026-09-26).

C_n: vertices 1..n, x ~ y iff x + y is a power of 2 or of 3 (targets >= 3).

Top targets.  P = largest power of 2 <= n, Q = largest power of 3 <= n, {T1 < T2} = {2P, 3Q}.
For n in a window (and n != 3^a) we have n < T1, T2 <= 2n, and m = T2 - T1 is odd.

Reduction (Lemma F_j of the note).  Put j = T1 - 1 - n >= 0 and r = T2 - 1 - n = m + j.
    alpha = M_{T1} on [1, n]   (pairs x + y = T1),
    beta  = M_{T2} on [1, n]   (pairs x + y = T2; its domain is the zone [r + 1, n]),
    residual R = [1, r] = F u W with F = [1, j] (free: no alpha, no beta) and W = [j + 1, j + m]
    (alpha, no beta).  phi(x) = (T1 - x) mod m, represented in W, is an involution of W with one fixed point u0.
For a set mu of C-edges inside R such that every W-vertex has exactly one mu-edge and every F-vertex two,
except one end e (one fewer), the graph alpha u beta u mu is a Hamiltonian path of C_n iff
phi u mu is a Hamiltonian path of R from u0 to e.  (Proof in the note; checked by `lift` + `verify_path`.)

Nothing here is trusted: every path produced is re-verified edge by edge by `verify_path`, which uses only
the definition of C_n.
"""
import sys
from array import array

EXP = '/tmp/math-wt-collatz-procgen-20260922/04-computation/experiments'
if EXP not in sys.path:
    sys.path.insert(0, EXP)

TG = []
_v = 4
while _v < 1 << 62:
    TG.append(_v)
    _v *= 2
_v = 3
while _v < 1 << 62:
    TG.append(_v)
    _v *= 3
TG.sort()
TGSET = frozenset(TG)


def is_target(t):
    return t in TGSET


def top_targets(n):
    P = 1
    while P * 2 <= n:
        P *= 2
    Q = 1
    while Q * 3 <= n:
        Q *= 3
    return P, Q, min(2 * P, 3 * Q), max(2 * P, 3 * Q)


def window(a):
    """Window W_a of THM-4505 (C2): (lo, hi, kind)."""
    B = 3 ** (a - 1)
    pw = [2 ** k for k in range(1, 4 * a + 8) if B < 2 ** k < 3 * B]
    if len(pw) == 2:
        p0, p1 = pw
        return max(3 * B - p0, p1 - B), 3 * B - 1, 'two-power'
    return 2 * B, 3 * B, 'one-power'


def reduced_data(n):
    """Parameters of the reduction at n (n in a window, n not a power of 3)."""
    P, Q, T1, T2 = top_targets(n)
    assert T1 > n and T2 <= 2 * n, (n, T1, T2)
    m = T2 - T1
    assert m % 2 == 1
    j = T1 - 1 - n
    r = m + j
    # fixed point of phi in W: 2 u0 = T1 (mod m)
    u0 = (T1 * ((m + 1) // 2)) % m
    while u0 < j + 1:
        u0 += m
    while u0 > j + m:
        u0 -= m
    f = T1 // 2 if T1 % 2 == 0 else T2 // 2      # the degree-one vertex of alpha u beta (a power of 2)
    return dict(n=n, P=P, Q=Q, T1=T1, T2=T2, m=m, j=j, r=r, u0=u0, f=f)


def phi_of(d, x):
    """phi(x) for x in W = [j+1, j+m]."""
    j, m = d['j'], d['m']
    y = (d['T1'] - x) % m
    while y < j + 1:
        y += m
    while y > j + m:
        y -= m
    return y


def lift(d, mu_adj=None, MU=None):
    """alpha u beta u mu as a vertex sequence from its degree-one vertex f.
    mu_adj: dict x -> list of mu-partners (x in R), or MU: array with MU[x] = the mu-partner (0 if none), x <= r.
    Returns array('i') of length <= n (not verified here)."""
    n, T1, T2, f, r = d['n'], d['T1'], d['T2'], d['f'], d['r']

    def nbrs(x):
        out = []
        y = T1 - x
        if 1 <= y <= n and y != x:
            out.append(y)
        y = T2 - x
        if 1 <= y <= n and y != x:
            out.append(y)
        if MU is not None:
            if x <= r and MU[x]:
                out.append(MU[x])
        else:
            out.extend(mu_adj.get(x, ()))
        return out

    seq = array('i', [f])
    prev, cur = 0, f
    steps = 0
    while steps <= n:
        nx = [y for y in nbrs(cur) if y != prev]
        if not nx:
            break
        if len(nx) > 1:
            raise ValueError('degree > 2 at %d' % cur)
        prev, cur = cur, nx[0]
        seq.append(cur)
        steps += 1
    return seq


def verify_path(n, seq):
    """Hamiltonian path of C_n?  Own check: a permutation of 1..n with every consecutive sum a target."""
    if len(seq) != n:
        return False
    seen = bytearray(n + 1)
    for x in seq:
        if x < 1 or x > n or seen[x]:
            return False
        seen[x] = 1
    for i in range(n - 1):
        if (seq[i] + seq[i + 1]) not in TGSET:
            return False
    return True


# ------------------------------------------------------------------------------------------------
# residual problem, exact search (uses the sumgraph lane's propagation/DPLL engine, read-only)
# ------------------------------------------------------------------------------------------------
def c_edges(lo, hi, members=None):
    """C-edges {x, y} with lo <= x < y <= hi (both in `members` if given)."""
    E = []
    for t in TG:
        if t > 2 * hi:
            break
        a = max(lo, t - hi)
        b = (t - 1) // 2
        for x in range(a, b + 1):
            y = t - x
            if members is None or (x in members and y in members):
                E.append((x, y, t))
    return E


def solve_residual_dpll(d, end_set=None, node_limit=200000):
    """Exact search on the residual problem of reduced_data d.  Returns mu_adj or None/'ABORT'."""
    import procgen_sumgraph_20260926_solver as SV
    r, j, m, u0 = d['r'], d['j'], d['m'], d['u0']
    E = [(x, y, t) for (x, y, t) in c_edges(1, r)]
    idx = {(x, y): i for i, (x, y, t) in enumerate(E)}
    phi_ids = []
    for x in range(j + 1, j + m + 1):
        y = phi_of(d, x)
        if x < y:
            if (x, y) in idx:
                phi_ids.append(idx[(x, y)])
            else:
                E.append((x, y, -1))
                phi_ids.append(len(E) - 1)
    Z = r + 1
    for x in range(1, r + 1):
        if x != u0 and (end_set is None or x in end_set):
            E.append((x, Z, -2))
    S = SV.State(r + 1, E, [u0, Z])
    q = []
    try:
        for i in phi_ids:
            S._force(i, 'BR', None, q)
        S.propagate(q + list(range(1, r + 2)))
    except SV.Contra:
        return None
    counter = [0]
    res = SV.portfolio(S, node_limit, counter)
    if res != 'FOUND':
        return 'ABORT' if res == 'ABORT' else None
    mu = {}
    for i in range(len(E)):
        if S.st[i] == SV.FO and E[i][2] > 0:
            x, y, t = E[i]
            mu.setdefault(x, []).append(y)
            mu.setdefault(y, []).append(x)
    # remove the phi pairs that happen to be C-edges (they are alpha-walk pairs, not mu edges)
    phiset = set(phi_ids)
    for i in phiset:
        x, y, t = E[i]
        if t > 0:
            mu[x].remove(y)
            mu[y].remove(x)
    return mu


# ------------------------------------------------------------------------------------------------
# residual problem with j = 0 (RP): renormalisation by zone steps (Lemma Z), exact search at the bottom
# ------------------------------------------------------------------------------------------------
class GP:
    """V (bytearray membership over [0, N]), psi (array: involution of V, psi[s] = s), s start, e forced end or 0."""

    def __init__(self, N, inV, psi, s, e=0, size=None):
        self.N, self.inV, self.psi, self.s, self.e = N, inV, psi, s, e
        self.size = sum(inV) if size is None else size

    def verts(self):
        inV = self.inV
        return [x for x in range(1, self.N + 1) if inV[x]]


def rp_gp(m, K):
    N = m
    inV = bytearray(b'\x01') * (N + 1)
    inV[0] = 0
    psi = array('i', [0]) * (N + 1)
    s = 0
    for x in range(1, m + 1):
        y = (K - x) % m
        y = m if y == 0 else y
        psi[x] = y
        if y == x:
            s = x
    return GP(N, inV, psi, s, 0, m)


def zone_step(g, t, extra):
    """mu = M_t on Z (pairs inside V avoiding the extra pairs' vertices and the forced end) plus `extra` (dict).
    Z is kept implicitly (bytearray inZ; partner = extra.get(y, t - y)).
    Returns (GP', (t, extra, inZ, |Z|)) or None (a walk that never leaves Z = a cycle)."""
    inV, psi, s, e = g.inV, g.psi, g.s, g.e
    N = g.N
    inZ = bytearray(N + 1)
    nz = 0
    for x in extra:
        inZ[x] = 1
        nz += 1
    lo = max(1, t - N)
    for x in range(lo, min(N, t - 1) + 1):
        if inV[x] and not inZ[x] and x != e:
            y = t - x
            if y != x and y <= N and inV[y] and y not in extra and y != e:
                inZ[x] = 1
                nz += 1
    if nz == 0:
        return None
    inV2 = bytearray(inV)
    for x in range(1, N + 1):
        if inZ[x]:
            inV2[x] = 0
    psi2 = array('i', [0]) * (N + 1)
    sp = 0
    visited = 0
    for x in range(1, N + 1):
        if not inV2[x] or psi2[x]:
            continue
        y = psi[x]
        if y == x:
            psi2[x] = x
            sp = x
            continue
        hit = False
        while inZ[y]:
            z = extra[y] if y in extra else t - y
            visited += 2
            if z == s:
                hit = True
                break
            y = psi[z]
        if hit:
            psi2[x] = x
            sp = x
        else:
            psi2[x] = y
            psi2[y] = x
    if visited != nz:
        return None
    return GP(N, inV2, psi2, sp, e, g.size - nz), (t, extra, inZ, nz)


def _partners(x, inV, N):
    out = []
    for t in TG:
        if t <= x:
            continue
        y = t - x
        if y > N:
            break
        if y != x and inV[y]:
            out.append((t, y))
    return out


def _isolated(g, cands):
    return [x for x in cands if g.inV[x] and x != g.e and not _partners(x, g.inV, g.N)]


def gp_dpll(g, node_limit=200000):
    import procgen_sumgraph_20260926_solver as SV
    V = g.verts()
    idx = {x: i + 1 for i, x in enumerate(V)}
    k = len(V)
    Vs = set(V)
    E = []
    seen = {}
    for x in V:
        for t in TG:
            if t <= x:
                continue
            y = t - x
            if y > V[-1]:
                break
            if y > x and y in Vs:
                seen[(x, y)] = len(E)
                E.append((idx[x], idx[y], t))
    phi_ids = []
    for x in V:
        y = g.psi[x]
        if x < y:
            if (x, y) in seen:
                phi_ids.append(seen[(x, y)])
            else:
                E.append((idx[x], idx[y], -1))
                phi_ids.append(len(E) - 1)
    Zv = k + 1
    for x in V:
        if x != g.s and (g.e == 0 or x == g.e):
            E.append((idx[x], Zv, -2))
    S = SV.State(k + 1, E, [idx[g.s], Zv])
    q = []
    try:
        for i in phi_ids:
            S._force(i, 'BR', None, q)
        S.propagate(q + list(range(1, k + 2)))
    except SV.Contra:
        return None
    counter = [0]
    if SV.portfolio(S, node_limit, counter) != 'FOUND':
        return None
    seq = S.extract_path()
    if seq[0] != idx[g.s]:
        seq = seq[::-1]
    seq = [V[i - 1] for i in seq[:-1]]
    mu = {}
    for i in range(0, len(seq) - 1, 2):
        mu[seq[i]] = seq[i + 1]
        mu[seq[i + 1]] = seq[i]
    return mu


def gp_solve(g, MU, small=2500, budget=None, repair_depth=6):
    """Zone steps from the top (target of the largest vertex, full zone; the power-of-2 fixed point of an even
    target is paired through another target, and vertices left without partners are re-paired along a short
    repair chain), exact search below `small`.  Heuristic.  On success writes the matching into MU
    (array, MU[x] = partner, 0 for the end) and returns True."""
    if budget is None:
        budget = [300]
    if g.size <= small:
        mu = gp_dpll(g)
        if mu is None:
            return False
        for x, y in mu.items():
            MU[x] = y
        return True
    N, inV = g.N, g.inV
    xmax = N
    while xmax > 0 and (not inV[xmax] or xmax == g.e):
        xmax -= 1
    cands = [t for t in TG if xmax < t < 2 * xmax and t - xmax <= N and inV[t - xmax] and t - xmax != g.e]
    for t in cands:
        f = t // 2
        opts = [{}]
        if t % 2 == 0 and f <= N and inV[f] and f != g.e:
            ws = [y for (tt, y) in _partners(f, inV, N) if tt != t and y != g.e]
            opts = [{f: w, w: f} for w in ws] + ([None] if g.e == 0 else [])
        for extra in opts:
            if budget[0] <= 0:
                return False
            budget[0] -= 1
            newe = g.e
            ex = dict(extra) if extra is not None else {}
            if extra is None:
                newe = f            # f becomes the end
            res = _zone_with_repair(g, t, ex, newe, repair_depth)
            if res is None:
                continue
            g2, zinfo = res
            if gp_solve(g2, MU, small, budget, repair_depth):
                t_, ex_, inZ, nz = zinfo
                for x in range(1, N + 1):
                    if inZ[x]:
                        MU[x] = ex_[x] if x in ex_ else t_ - x
                return True
            del g2, zinfo, res
    return False


def _zone_with_repair(g, t, ex, newe, depth):
    """Zone step; if it leaves vertices without partners, re-pair one of them with a zone vertex through
    another target (the zone vertex's M_t partner returns to V) and retry, up to `depth` times."""
    gg = GP(g.N, g.inV, g.psi, g.s, newe, g.size)
    for _ in range(depth + 1):
        r = zone_step(gg, t, ex)
        if r is None:
            return None
        g2, zinfo = r
        inZ = zinfo[2]
        cands = set()
        for x in TG:
            if x > g.N:
                break
            cands.add(x)
        for x in ex:
            y = t - x
            if 1 <= y <= g.N:
                cands.add(y)
        cands.add(t // 2)
        # leftovers above the new residual interval (earlier extras, repaired mirrors): few vertices
        lo_zone = max(1, t - g.N)
        cands.update(x for x in range(lo_zone, g.N + 1) if g2.inV[x])
        iso = _isolated(g2, sorted(c for c in cands if 1 <= c <= g.N))
        if not iso:
            if g2.e and not g2.inV[g2.e]:
                return None
            return g2, zinfo
        if len(iso) == 1 and g2.e == 0 and iso[0] != g2.s:
            g2.e = iso[0]
            return g2, zinfo
        y = iso[0]
        opts = [z for (tt, z) in _partners(y, g.inV, g.N) if tt != t and inZ[z] and z not in ex]
        if not opts:
            return None
        z = opts[0]
        ex[y] = z
        ex[z] = y
    return None


def check_gp(g0, MU):
    """psi u mu is a Hamiltonian path of V from s (mu first, MU[x] = 0 for the end).  Own walk."""
    s = g0.s
    seen = bytearray(g0.N + 1)
    seen[s] = 1
    cnt = 1
    cur = s
    while MU[cur]:
        y = MU[cur]
        if seen[y] or MU[y] != cur:
            return False
        seen[y] = 1
        z = g0.psi[y]
        if z == y or seen[z]:
            return False
        seen[z] = 1
        cnt += 2
        cur = z
    return cnt == g0.size


# ------------------------------------------------------------------------------------------------
# exact reduction by contracting the forced top structure (Theorem F', all regimes)
# ------------------------------------------------------------------------------------------------
def forced_structure(n):
    """Forced edges of C_n: both edges of every top vertex (Lemma T) and the edge of the leaf P.
    Returns (P, forced adjacency dict, fdeg bytearray)."""
    P, Q, T1, T2 = top_targets(n)
    M = max(P, Q)
    fadj = {}

    def add(x, y):
        fadj.setdefault(x, set()).add(y)
        fadj.setdefault(y, set()).add(x)
    for x in range(M + 1, n + 1):
        for t in (2 * P, 3 * Q):
            y = t - x
            if 1 <= y <= n and y != x:
                add(x, y)
    # the leaf P (THM-4505 C1): its unique neighbour
    nb = [t - P for t in TG if P < t <= P + n and t - P != P]
    assert len(nb) == 1, (n, nb)
    add(P, nb[0])
    return P, fadj


def solve_exact_reduced(n, node_limit=200000):
    """Contract the forced chains; exact search on the reduced graph H_n.  Returns path (array) or None/'ABORT'.
    Covers the Hamiltonian paths of C_n whose free end is not saturated by forced edges."""
    import procgen_sumgraph_20260926_solver as SV
    P, fadj = forced_structure(n)
    fdeg = {v: len(s) for v, s in fadj.items()}
    if any(d > 2 for d in fdeg.values()):
        return None
    # walk chains
    seen = set()
    virt = []
    start = None
    for v in fadj:
        if fdeg[v] == 1 and v not in seen:
            prev, cur = None, v
            seen.add(cur)
            while True:
                nx = [w for w in fadj[cur] if w != prev]
                if not nx:
                    break
                prev, cur = cur, nx[0]
                if cur in seen:
                    return None           # forced cycle
                seen.add(cur)
            u, w = v, cur
            if P in (u, w):
                start = w if u == P else u
            else:
                virt.append((u, w))
    VH = [v for v in range(1, n + 1) if fdeg.get(v, 0) <= 1 and v != P and not (fdeg.get(v, 0) == 1 and v == start and False)]
    VH = [v for v in VH if not (fdeg.get(v, 0) == 2)]
    idx = {v: i + 1 for i, v in enumerate(VH)}
    k = len(VH)
    E = []
    seenE = set()
    for x in VH:
        for t in TG:
            if t <= x:
                continue
            y = t - x
            if y > n:
                break
            if y > x and y in idx and y not in fadj.get(x, ()):
                E.append((idx[x], idx[y], t))
                seenE.add((x, y))
    vid = []
    for (u, w) in virt:
        E.append((idx[u], idx[w], -1))
        vid.append(len(E) - 1)
    Z = k + 1
    for x in VH:
        if x != start:
            E.append((idx[x], Z, -2))
    S = SV.State(k + 1, E, [idx[start], Z])
    q = []
    try:
        for i in vid:
            S._force(i, 'BR', None, q)
        S.propagate(q + list(range(1, k + 2)))
    except SV.Contra:
        return None
    counter = [0]
    res = SV.portfolio(S, node_limit, counter)
    if res != 'FOUND':
        return 'ABORT' if res == 'ABORT' else None
    chosen = {}
    for i, (a_, b_, t) in enumerate(E):
        if S.st[i] == SV.FO and t > 0:
            x, y = VH[a_ - 1], VH[b_ - 1]
            chosen.setdefault(x, []).append(y)
            chosen.setdefault(y, []).append(x)
    # lift: forced edges + chosen edges, walk from P
    seq = array('i', [P])
    prev, cur = 0, P
    for _ in range(n):
        nx = [w for w in list(fadj.get(cur, ())) + chosen.get(cur, []) if w != prev]
        if not nx:
            break
        if len(nx) > 1:
            return None
        prev, cur = cur, nx[0]
        seq.append(cur)
    return seq
