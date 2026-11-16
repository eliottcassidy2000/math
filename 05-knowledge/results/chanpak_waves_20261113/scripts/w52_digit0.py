"""Digit-0 vertices: B(T,t) = sum over paths P starting at t of H(T - V(P)) (HYP-9453's B), digit 0 <=> B(T,t) = H(T).
Test on strong tournaments of order <= 7: is every digit-0 vertex a near-join tip (out-degree 1, its out-neighbour t0
before a transitive tail: T - t = T0 => M_k with t0 in T0), and does every near-join tip have digit 0?"""
import itertools, sys
sys.path.insert(0, '.')
from w45_hk import ham_split, sub
def classes(n):
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]; m = len(pairs)
    idx = {pr: t for t, pr in enumerate(pairs)}; perms = list(itertools.permutations(range(n)))
    seen = bytearray(1 << m)
    def relabel(bits, p):
        b = 0
        for t, (i, j) in enumerate(pairs):
            x = bits >> t & 1; a, c = p[i], p[j]
            if a < c: b |= x << idx[(a, c)]
            else: b |= (1 - x) << idx[(c, a)]
        return b
    for bits in range(1 << m):
        if seen[bits]: continue
        for p in perms: seen[relabel(bits, p)] = 1
        out = [0] * n
        for t, (i, j) in enumerate(pairs):
            if bits >> t & 1: out[i] |= 1 << j
            else: out[j] |= 1 << i
        yield out
def strong(n, out):
    def reach(fwd):
        r = 1; st = [0]
        while st:
            u = st.pop()
            nb = out[u] if fwd else sum(1 << w for w in range(n) if out[w] >> u & 1)
            t = nb & ~r
            while t:
                w = (t & -t).bit_length() - 1; t &= t - 1; r |= 1 << w; st.append(w)
        return r
    return reach(True) == (1 << n) - 1 and reach(False) == (1 << n) - 1
def B(n, out, t):
    tot = 0
    def ext(path, used):
        nonlocal tot
        rest = [v for v in range(n) if not used >> v & 1]
        tot += ham_split(*sub(out, rest))[0] if rest else 1
        for w in range(n):
            if not used >> w & 1 and out[path[-1]] >> w & 1: ext(path + [w], used | 1 << w)
    ext([t], 1 << t); return tot
def near_join_tip(n, out, t):
    if bin(out[t]).count('1') != 1: return False
    t0 = out[t].bit_length() - 1
    rest = [v for v in range(n) if v != t]
    # T - t must be T0 => M_k (transitive tail of k >= 1 vertices, all beaten by T0), with t0 in T0
    for k in range(1, len(rest)):
        for tail in itertools.combinations(rest, k):
            head = [v for v in rest if v not in tail]
            if t0 in tail: continue
            if not all(out[h] >> s & 1 for h in head for s in tail): continue
            # tail transitive
            if all(sum(1 for b in tail if b != a and out[a] >> b & 1) in range(k) for a in tail) and \
               sorted(sum(1 for b in tail if b != a and out[a] >> b & 1) for a in tail) == list(range(k)):
                return True
    return False
stats = {}
for n in range(3, 8):
    d0 = d0_tip = tips = tips_d0 = 0; counts_d0 = set(); counts_all = set()
    for out in classes(n):
        if not strong(n, out): continue
        H = ham_split(n, out)[0]; counts_all.add(H)
        for t in range(n):
            b = B(n, out, t); tip = near_join_tip(n, out, t)
            if b == H:
                d0 += 1; counts_d0.add(H); d0_tip += tip
            if tip:
                tips += 1; tips_d0 += (b == H)
    print(f"n={n}: digit-0 vertices {d0} (near-join tips among them {d0_tip}); near-join tips {tips} (digit-0 among them {tips_d0});"
          f" strong counts without a digit-0 vertex at this order: {sorted(counts_all - counts_d0)[:12]}")
