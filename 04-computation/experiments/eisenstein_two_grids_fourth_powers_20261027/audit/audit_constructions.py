"""Audit of THM-4824 (A1)-(A3) as written, plus the auditor's extensions (A_n support reduction, E7/E6 index-2 glue).

Each construction is a function decompose(x) -> None (x is a unit leaf) or a split (y, z) with y + z = x.  The verifier
unfolds the decomposition from a target x (as a DAG over values: identical values have identical subtrees), checks
  * every leaf is a minimal vector (unit charge),
  * every internal node has y + z = x and beta(y, z) <= 0 (so the merge costs beta_+ = 0),
  * the unfolding is finite (no value is re-entered while in progress),
so the unfolded tree is a zero-cost merge tree with root x.

Lattices (exact integer models):
  Z[i]   : Z^2, Q = a^2 + b^2                                       units +-1, +-i
  D_n    : Z^n, sum even, Q = x.x/2                                 roots +-e_i +-e_j
  A_n    : Z^(n+1), sum 0, Q = x.x/2                                roots e_i - e_j
  E8     : doubled coords y = 2x, one parity, sum = 0 mod 4, Q = y.y/8
  E7     : doubled coords, one parity, sum 0, Q = y.y/8           (= A7 u (A7 + (1/2^4, -1/2^4)))
  E6     : doubled (a_1..a_6; t), one parity, sum a = 0, Q = (a.a + 2t^2)/8   (= A5+A1 u glue)
"""
import random, sys

sys.setrecursionlimit(100000)


def verify(target, decompose, is_unit, beta, add, sub, zero):
    state = {}          # value -> 'busy' | 'done'
    stack = [(target, 0)]
    nodes = leaves = 0
    while stack:
        v, phase = stack.pop()
        if phase == 0:
            if state.get(v) == 'done':
                continue
            if state.get(v) == 'busy':
                raise AssertionError(f"cycle at {v}")
            if is_unit(v):
                state[v] = 'done'; leaves += 1
                continue
            sp = decompose(v)
            assert sp is not None, f"non-unit {v} has no split"
            y, z = sp
            assert add(y, z) == v, f"bad split {v} -> {y} + {z}"
            assert y != v and z != v, f"trivial split at {v}"
            b = beta(y, z)
            assert b <= 0, f"paid merge at {v}: {y} + {z}, beta = {b}"
            state[v] = 'busy'
            stack.append((v, 1))
            stack.append((y, 0)); stack.append((z, 0))
        else:
            state[v] = 'done'; nodes += 1
    return nodes, leaves


# ---------------------------------------------------------------- square grid, literally as THM-4824 (A1)
def Zi_decompose(x):
    a, b = x
    if (a, b) == (0, 0):
        return ((1, 0), (-1, 0))
    if b == 0:
        if a < 0:
            y, z = Zi_decompose((-a, 0)); return ((-y[0], -y[1]), (-z[0], -z[1]))
        if a == 1:
            return None
        if a % 2 == 0:
            k = a // 2
            return ((k, k), (k, -k))                    # k(1+i) + k(1-i), beta = 0
        k = (a - 1) // 2
        return ((k + 1, k + 1), (k, -(k + 1)))          # (k+1)(1+i) + (k-(k+1)i), beta = -2(k+1)
    if a == 0:                                          # imaginary: rotate the real construction by i
        sp = Zi_decompose((b, 0))
        if sp is None:
            return None
        rot = lambda p: (-p[1], p[0])
        return (rot(sp[0]), rot(sp[1]))
    return ((a, 0), (0, b))                             # perpendicular merge a + bi


Zi = dict(is_unit=lambda v: v[0]*v[0] + v[1]*v[1] == 1,
          beta=lambda y, z: 2*(y[0]*z[0] + y[1]*z[1]),
          add=lambda y, z: (y[0]+z[0], y[1]+z[1]), sub=None, zero=(0, 0))


# ---------------------------------------------------------------- D_n, literally as THM-4824 (A2) (n >= 3)
def tadd(y, z): return tuple(a + b for a, b in zip(y, z))
def tsub(y, z): return tuple(a - b for a, b in zip(y, z))
def unitvec(n, i, s=1): return tuple(s if k == i else 0 for k in range(n))


def D_decompose(x):
    n = len(x)
    if sum(a*a for a in x) == 2:
        return None
    if not any(x):
        r = unitvec(n, 0); r = tadd(r, unitvec(n, 1)); return (r, tuple(-a for a in r))
    odd = [i for i in range(n) if x[i] % 2]
    evn = [i for i in range(n) if x[i] and x[i] % 2 == 0]
    blocks = [(odd[k], odd[k+1]) for k in range(0, len(odd), 2)] + [(i,) for i in evn]
    def vec(block):
        v = [0]*n
        for i in block: v[i] = x[i]
        return tuple(v)
    if len(blocks) >= 2:                                 # disjoint supports: perpendicular
        y = vec(blocks[0]); return (y, tsub(x, y))
    blk = blocks[0]
    if len(blk) == 1:                                    # even multiple j e_a = (j/2)(e_a+e_c) + (j/2)(e_a-e_c)
        i = blk[0]; c = 0 if i else 1; h = x[i] // 2
        y = [0]*n; y[i] = h; y[c] = h
        z = [0]*n; z[i] = h; z[c] = -h
        return (tuple(y), tuple(z))
    i, j = blk
    if abs(x[i]) == abs(x[j]):                           # odd root multiple j(e_a +- e_b), j >= 3, third index c
        c = next(k for k in range(n) if k not in (i, j))
        y = [0]*n; y[i] = x[i]; y[c] = 1
        z = [0]*n; z[j] = x[j]; z[c] = -1
        return (tuple(y), tuple(z))
    p, q = (x[i] + x[j]) // 2, (x[i] - x[j]) // 2        # pair of odd coordinates: perpendicular root multiples
    y = [0]*n; y[i] = p; y[j] = p
    z = [0]*n; z[i] = q; z[j] = -q
    return (tuple(y), tuple(z))


def Dn_ops(n):
    return dict(is_unit=lambda v: sum(a*a for a in v) == 2,
                beta=lambda y, z: sum(a*b for a, b in zip(y, z)),        # beta = 2<y,z> = y.z for Q = x.x/2
                add=tadd, sub=tsub, zero=tuple([0]*n))


# ---------------------------------------------------------------- E8, literally as THM-4824 (A3), doubled coordinates
def E8_decompose(y):
    if sum(a*a for a in y) == 8:
        return None
    if all(a % 2 == 0 for a in y):                       # D8 point
        sp = D_decompose(tuple(a // 2 for a in y))
        return None if sp is None else (tuple(2*a for a in sp[0]), tuple(2*a for a in sp[1]))
    h = [-1 if a > 0 else 1 for a in y]                  # h' = -1/2 sign(x), doubled
    if sum(1 for a in h if a < 0) % 2:
        j = min(range(8), key=lambda k: abs(y[k])); h[j] = -h[j]
    h = tuple(h)
    return (tsub(y, h), h)


def E8_beta(y, z):
    s = sum(a*b for a, b in zip(y, z))
    assert s % 4 == 0, "non-integral beta"
    return s // 4


E8_ops = dict(is_unit=lambda v: sum(a*a for a in v) == 8, beta=E8_beta, add=tadd, sub=tsub, zero=(0,)*8)


# ---------------------------------------------------------------- auditor: A_n (n >= 3) by support reduction
def square_grid_embed(k, alpha, h):
    """Decomposition of k*alpha in the square grid Z alpha + Z h (alpha, h perpendicular roots)."""
    def emb(p): return tuple(p[0]*a + p[1]*b for a, b in zip(alpha, h))
    return emb


def zero_sum_split(x):
    """A proper nonempty subset S of the support with sum_{S} x = 0 (then x = x|S + x|S^c, disjoint supports, beta = 0)."""
    from itertools import combinations
    supp = [i for i in range(len(x)) if x[i]]
    for r in range(2, len(supp) - 1):
        for S in combinations(supp, r):
            if sum(x[i] for i in S) == 0:
                return S
    return None


def A_decompose_factory(m):
    """A_{m-1} in Z^m, m >= 4.  Rules: (3) perpendicular split along a zero-sum partition of the support;
    (4) obtuse support reduction y = x_i (e_i - e_t), x_i, x_t of the same sign, beta = -x_i (x_i + x_t) < 0;
    (5) root multiple k (e_i - e_j): square-grid (A1) construction in Z(e_i - e_j) + Z(e_c - e_d)."""
    def dec(x):
        if sum(a*a for a in x) == 2:
            return None
        if not any(x):
            r = [0]*m; r[0] = 1; r[1] = -1; r = tuple(r); return (r, tuple(-a for a in r))
        S = zero_sum_split(x)
        if S is not None:
            y = tuple(x[i] if i in S else 0 for i in range(m)); return (y, tsub(x, y))
        supp = [i for i in range(m) if x[i]]
        if len(supp) >= 3:
            pos = [i for i in supp if x[i] > 0]; neg = [i for i in supp if x[i] < 0]
            i, t = (pos[0], pos[1]) if len(pos) >= 2 else (neg[0], neg[1])
            y = [0]*m; y[i] = x[i]; y[t] = -x[i]
            y = tuple(y); return (y, tsub(x, y))
        i, j = supp
        k = x[i]
        c, d = [l for l in range(m) if l not in (i, j)][:2]
        alpha = [0]*m; alpha[i] = 1; alpha[j] = -1
        h = [0]*m; h[c] = 1; h[d] = -1
        emb = lambda p: tuple(p[0]*a + p[1]*b for a, b in zip(alpha, h))
        sp = Zi_decompose((k, 0))
        return (emb(sp[0]), emb(sp[1]))
    return dec


def A_decompose_general(m):
    """Full A-lattice decomposition: points supported anywhere (including off-axis square-grid nodes)."""
    base = A_decompose_factory(m)
    def dec(x):
        if sum(a*a for a in x) == 2:
            return None
        return base(x)
    return dec


An_ops = Dn_ops  # same form x.x/2


# ---------------------------------------------------------------- auditor: E7 = A7 u (A7 + glue), doubled coordinates
def E7_decompose(y):
    if sum(a*a for a in y) == 8:
        return None
    if all(a % 2 == 0 for a in y):
        sp = A_decompose_factory(8)(tuple(a // 2 for a in y))
        return None if sp is None else (tuple(2*a for a in sp[0]), tuple(2*a for a in sp[1]))
    order = sorted(range(8), key=lambda k: y[k])
    h = [0]*8
    for k in order[:4]: h[k] = 1                         # +1/2 on the four smallest, -1/2 on the four largest
    for k in order[4:]: h[k] = -1
    h = tuple(h)
    return (tsub(y, h), h)


# ---------------------------------------------------------------- auditor: E6 = (A5 + A1) u glue, doubled (a; t)
def E6_beta(y, z):
    s = sum(a*b for a, b in zip(y[:6], z[:6])) + 2*y[6]*z[6]
    assert s % 4 == 0
    return s // 4


def E6_decompose(y):
    if sum(a*a for a in y[:6]) + 2*y[6]*y[6] == 8:
        return None
    if not any(y):
        r = (2, -2, 0, 0, 0, 0, 0); return (r, tuple(-a for a in r))
    if all(a % 2 == 0 for a in y):                       # M = A5 + A1
        a, t = y[:6], y[6]
        if any(a) and t:
            return (a + (0,), (0,)*6 + (t,))             # perpendicular
        if t == 0:
            sp = A_decompose_factory(6)(tuple(v // 2 for v in a))
            return None if sp is None else (tuple(2*v for v in sp[0]) + (0,), tuple(2*v for v in sp[1]) + (0,))
        # t*alpha, alpha = (0;1) [doubled (0^6; 2)], perpendicular root e1 - e2 [doubled (2,-2,0,0,0,0; 0)]
        k = t // 2
        alpha = (0,)*6 + (2,); h = (2, -2, 0, 0, 0, 0, 0)
        emb = lambda p: tuple(p[0]*u + p[1]*v for u, v in zip(alpha, h))
        sp = Zi_decompose((k, 0))
        return (emb(sp[0]), emb(sp[1]))
    h0 = (1, 1, 1, -1, -1, -1, 1)                        # a glue root: (6 + 2)/8 = 1
    h = h0 if E6_beta(y, h0) + 2 <= 0 else tuple(-a for a in h0)   # pick sign with x.h <= 0
    return (tsub(y, h), h)


E6_ops = dict(is_unit=lambda v: sum(a*a for a in v[:6]) + 2*v[6]*v[6] == 8, beta=E6_beta, add=tadd, sub=tsub, zero=(0,)*7)
E7_ops = dict(E8_ops)


def E6_decompose_full(y):
    """Square-grid nodes p0*alpha + p1*h of a t*alpha decomposition are general M points; route them through E6_decompose
    is fine except that a node (a; t) with both parts nonzero splits perpendicularly, and pure-alpha / pure-A5 parts
    recurse.  The only extra case: nodes of the form p0*alpha + p1*h are (2p1, -2p1, 0,0,0,0; 2p0) -> handled above."""
    return E6_decompose(y)


# ---------------------------------------------------------------- driver
def rand_vec(n, R, cond):
    while True:
        v = tuple(random.randint(-R, R) for _ in range(n))
        if cond(v): return v


def main():
    random.seed(20261027)
    tot = 0
    # (A1) Z[i]
    for a in range(-60, 61):
        for b in range(-60, 61):
            verify((a, b), Zi_decompose, **Zi); tot += 1
    for _ in range(200):
        verify(rand_vec(2, 5000, lambda v: True), Zi_decompose, **Zi); tot += 1
    print(f"(A1) Z[i]: zero-cost trees verified for all |a|,|b| <= 60 and 200 random points with |coords| <= 5000")
    # (A2) D_n literal
    for n in range(3, 9):
        ops = Dn_ops(n)
        for _ in range(400):
            v = rand_vec(n, 60, lambda v: sum(v) % 2 == 0)
            verify(v, D_decompose, **ops)
        for j in range(1, 80):                           # root multiples j(e_a +- e_b) for every sign
            for s in (1, -1):
                v = [0]*n; v[0] = j; v[1] = s*j
                verify(tuple(v), D_decompose, **ops)
        print(f"(A2) D{n}: literal construction verified on 400 random points (|coords| <= 60) and all j(e1 +- e2), j < 80")
    # (A3) E8 literal (doubled coordinates)
    worst = 0.0
    for _ in range(3000):
        par = random.randint(0, 1)
        v = rand_vec(8, 40, lambda v: all(a % 2 == par for a in v) and sum(v) % 4 == 0)
        verify(v, E8_decompose, **E8_ops)
        if par == 1:
            y, h = E8_decompose(v)
            b = E8_ops['beta'](y, h)
            bound = -(3/8) * sum(abs(a) / 2 for a in v) - 2
            assert b <= bound + 1e-9, (v, b, bound)
            worst = max(worst, b - bound)
    print(f"(A3) E8: literal glue construction verified on 3000 random points; beta <= -(3/8) sum|x_i| - 2 holds (max gap {worst:.3f})")
    # auditor: A_n
    for m in range(4, 11):
        ops = Dn_ops(m)
        dec = A_decompose_factory(m)
        for _ in range(300):
            v = list(rand_vec(m - 1, 50, lambda v: True)); v.append(-sum(v)); v = tuple(v)
            verify(v, dec, **ops)
        print(f"(aud) A{m-1}: support reduction + square-grid root multiples verified on 300 random points (|coords| <= 50)")
    # auditor: E7, E6
    for _ in range(2000):
        par = random.randint(0, 1)
        v = rand_vec(8, 40, lambda v: all(a % 2 == par for a in v) and sum(v) == 0)
        verify(v, E7_decompose, **E7_ops)
    print("(aud) E7: index-2 glue over A7 verified on 2000 random points")
    for _ in range(2000):
        par = random.randint(0, 1)
        v = rand_vec(7, 40, lambda v: all(a % 2 == par for a in v) and sum(v[:6]) == 0)
        verify(v, E6_decompose, **E6_ops)
    print("(aud) E6: index-2 glue over A5+A1 verified on 2000 random points")
    # auditor: reducible A2+A1, model Z^3 (sum 0) + Z t, Q = (x.x)/2 + t^2 ; doubled-free integer form
    def A2A1_beta(y, z): return sum(a*b for a, b in zip(y[:3], z[:3])) + 2*y[3]*z[3]
    def A2A1_unit(v): return sum(a*a for a in v[:3]) + 2*v[3]*v[3] == 2
    def A2A1_dec(x):
        if A2A1_unit(x): return None
        a, t = x[:3], x[3]
        if not any(x): return ((1, -1, 0, 0), (-1, 1, 0, 0))
        if any(a) and t: return (a + (0,), (0, 0, 0, t))
        alpha = (0, 0, 0, 1)
        if t:                                            # t*alpha in square grid with root e1 - e2
            h = (1, -1, 0, 0)
        else:
            supp = [i for i in range(3) if a[i]]
            if len(supp) == 3:
                pos = [i for i in supp if a[i] > 0]; neg = [i for i in supp if a[i] < 0]
                i, s_ = (pos[0], pos[1]) if len(pos) >= 2 else (neg[0], neg[1])
                y = [0, 0, 0, 0]; y[i] = a[i]; y[s_] = -a[i]; y = tuple(y)
                return (y, tsub(x, y))
            i, j = supp
            k = a[i]
            alpha = [0, 0, 0, 0]; alpha[i] = 1; alpha[j] = -1; alpha = tuple(alpha)
            h = (0, 0, 0, 1)
            emb = lambda p: tuple(p[0]*u + p[1]*v for u, v in zip(alpha, h))
            sp = Zi_decompose((k, 0)); return (emb(sp[0]), emb(sp[1]))
        emb = lambda p: tuple(p[0]*u + p[1]*v for u, v in zip(alpha, h))
        sp = Zi_decompose((t, 0)); return (emb(sp[0]), emb(sp[1]))
    ops = dict(is_unit=A2A1_unit, beta=A2A1_beta, add=tadd, sub=tsub, zero=(0, 0, 0, 0))
    for _ in range(2000):
        v = list(rand_vec(2, 60, lambda v: True)); v.append(-sum(v)); v.append(random.randint(-60, 60))
        verify(tuple(v), A2A1_dec, **ops)
    print("(aud) A2+A1: free via perpendicular-root square grids + support reduction, verified on 2000 random points")


if __name__ == '__main__':
    main()
