"""
lrc_graph_interpretation_ladder_opus_S141.py  (clean rewrite; v1 had window-offset + wrap
bugs, both fixed and documented below)

THE GRAPH-THEORETIC LADDER FOR LRC (owner: a purely graph-theoretical interpretation).

Objects: for finite S subset Z_{>0}: the distance graph G_S = Cay(Z, +-S);
  M(S) = sup_t min_s ||s t||;  mu(S) = Motzkin density (max independence density of G_S);
  chi_f(G_S) = 1/mu(S);  chi_c = circular chromatic number.

LADDER:  LRC(14)  ==>  GRAPH-14 [chi_c(G_S) <= 14 for all |S|=13]  ==>  MOTZKIN-14
[mu(S) >= 1/14].  Both implications via the witness coloring (L1).  The converse gaps:
LINEARIZATION (chi_c = 1/M ?) and mu = M ?.  Empirical finding of this session: the
periodic Motzkin optimum EQUALS M on every set tested (raising the collapse possibility).

RUNS: (1) mu (periodic, N<=240, exact transfer DP) vs M (exact) on small sets;
(2) repo extremal 13-families: exact witness reconstruction + certified independent set of
density exactly M (so mu >= M >= 1/13 > 1/14 on all of them);
(3) chi_c probes: search for circular colorings strictly better than 1/M on tiny S.
"""
from fractions import Fraction as F
from math import gcd
import time, itertools

def M_exact(S):
    S = sorted(S)
    cands = set()
    for i in range(len(S)):
        cands.add(2 * S[i])
        for j in range(i, len(S)):
            cands.add(S[i] + S[j])
            if S[j] > S[i]:
                cands.add(S[j] - S[i])
    best = F(0)
    for q in sorted(cands):
        if q < 2: continue
        for a in range(1, q):
            if gcd(a, q) != 1: continue
            v = min(min((s * a) % q, q - (s * a) % q) for s in S)
            if F(v, q) > best:
                best = F(v, q)
    return best

def witness(S):
    """(a, q, v) with min_s ||s a/q|| = v/q = M(S)."""
    M = M_exact(S)
    cands = set()
    for i in range(len(S)):
        cands.add(2 * S[i])
        for j in range(i, len(S)):
            cands.add(S[i] + S[j])
            if S[j] > S[i]: cands.add(S[j] - S[i])
    for q in sorted(cands):
        if q < 2: continue
        for a in range(1, q):
            if gcd(a, q) != 1: continue
            v = min(min((s * a) % q, q - (s * a) % q) for s in S)
            if F(v, q) == M:
                return a, q, v
    return None

def motzkin_periodic(S, N):
    """Exact max density of a period-N set avoiding differences in S (cyclic)."""
    w = max(S)
    if N <= 26:
        best = 0
        for mask in range(1 << N):
            ok = True
            for i in range(N):
                if not (mask >> i) & 1: continue
                for ss in S:
                    if (mask >> ((i + ss) % N)) & 1:
                        ok = False; break
                if not ok: break
            if ok:
                c = bin(mask).count("1")
                if c > best: best = c
        return F(best, N)
    if N <= 2 * w:
        return F(0, 1)
    # window: bit j of mask = membership at position p-w+j (p = next position).
    # appending 1 at p forbidden iff membership at p-s <=> mask bit (w-s).
    Smask = 0
    for ss in S: Smask |= (1 << (w - ss))
    states = []
    def gen(mask, j):
        if j == w:
            states.append(mask); return
        gen(mask, j + 1)
        for ss in S:
            jj = j - ss
            if jj >= 0 and (mask >> jj) & 1:
                break
        else:
            gen(mask | (1 << j), j + 1)
    gen(0, 0)
    sidx = {m: i for i, m in enumerate(states)}
    trans = []
    for m in states:
        t0 = sidx.get(m >> 1)
        nm1 = (m >> 1) | (1 << (w - 1))
        t1 = sidx.get(nm1) if (m & Smask) == 0 else None
        trans.append((t0, t1))
    NEG = -10**9
    best_overall = -1
    for start_i, hmask in enumerate(states):
        cur = [NEG] * len(states)
        cur[start_i] = bin(hmask).count("1")
        for _ in range(N - w):
            nxt = [NEG] * len(states)
            for si, val in enumerate(cur):
                if val == NEG: continue
                t0, t1 = trans[si]
                if t0 is not None and val > nxt[t0]: nxt[t0] = val
                if t1 is not None and val + 1 > nxt[t1]: nxt[t1] = val + 1
            cur = nxt
        for fi, fmask in enumerate(states):
            if cur[fi] == NEG: continue
            ok = True
            for ss in S:
                for x in range(N - ss, N):
                    jt = x - (N - w)
                    if jt < 0 or jt >= w: continue
                    if (fmask >> jt) & 1:
                        xh = x + ss - N
                        if 0 <= xh < w and (hmask >> xh) & 1:
                            ok = False; break
                if not ok: break
            if ok and cur[fi] > best_overall:
                best_overall = cur[fi]
    return F(max(best_overall, 0), N)

def motzkin_scan(S, Nmax=240):
    best = F(0); bestN = None
    for N in range(2, Nmax + 1):
        v = motzkin_periodic(S, N)
        if v > best:
            best, bestN = v, N
    return best, bestN

def main():
    t0 = time.time()
    print("=" * 100)
    print("(1) periodic Motzkin (exact, N <= 240) vs M (exact)")
    print("=" * 100)
    tests = {
        "AP {1..4}": [1, 2, 3, 4], "{1,2}": [1, 2], "{2,3}": [2, 3],
        "{1,4,5}": [1, 4, 5], "{2,3,5}": [2, 3, 5], "{1,3,4}": [1, 3, 4],
        "{3,5,8}": [3, 5, 8], "{1,2,5,6}": [1, 2, 5, 6],
    }
    for name, S in tests.items():
        M = M_exact(S)
        mu, N = motzkin_scan(S, 240)
        rel = "mu = M" if mu == M else ("mu > M  <-- FRACTIONAL SLACK" if mu > M else "*** mu < M: BUG ***")
        print(f"  {name:<12} M = {str(M):>6} = {float(M):.4f}   mu(<=240) = {str(mu):>6} = {float(mu):.4f} (N={N})   {rel}")

    print()
    print("=" * 100)
    print("(2) repo extremal 13-families: witness reconstruction + certified independent set")
    print("=" * 100)
    fams = {
        "GW {1..11,13,24}": list(range(1, 12)) + [13, 24],
        "prim-sat 2AP+13": [2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 13],
        "parity record": [2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 11, 13],
        "deep well {1..12,182}": list(range(1, 13)) + [182],
    }
    for name, S in fams.items():
        M = M_exact(S)
        a, q, v = witness(S)
        A = [x for x in range(q) if (a * x) % q < v]
        okA = True
        for x in A:
            for s2 in S:
                if (x + s2) % q in A: okA = False
        dens = F(len(A), q)
        print(f"  {name:<22} M = {str(M):>7} = {float(M):.5f}  witness (a,q,v)=({a},{q},{v})"
              f"  indep set density = {dens} (certified: {okA})"
              f"  {'>= 1/14 OK' if dens >= F(1, 14) else '*** below ***'}")

    print()
    print("=" * 100)
    print("(3) chi_c probes on tiny S: any circular coloring strictly better than 1/M?")
    print("=" * 100)
    for S in ([1, 2], [2, 3], [1, 4, 5]):
        M = M_exact(S)
        found = None
        for p in range(2, 14):
            for qq in range(1, p // 2 + 1):
                if gcd(p, qq) != 1: continue
                if F(p, qq) >= 1 / M: continue
                for N in range(2, 13):
                    if p ** N > 2 * 10**6: continue
                    for coloring in itertools.product(range(p), repeat=N):
                        good = True
                        for x in range(N):
                            for s in S:
                                d = (coloring[(x + s) % N] - coloring[x]) % p
                                if min(d, p - d) < qq:
                                    good = False; break
                            if not good: break
                        if good:
                            found = (p, qq, N); break
                    if found: break
                if found: break
            if found: break
        print(f"  S={S}: 1/M = {float(1/M):.3f}; coloring with p/q < 1/M: {found}"
              f"  {'<-- SEPARATION!' if found else '(none: consistent with chi_c = 1/M)'}")

    print(f"\n[{time.time()-t0:.0f}s]")

if __name__ == "__main__":
    main()
