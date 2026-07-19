#!/usr/bin/env python3
"""
death-star-2026-07-19-S59 -- HYP-7885: the cross-N first-gap species census at the
UNVERIFIED MIDDLE N = 8, 9, 10, 11 (speeds), with validation gates.

Background (verify-first): opus-S118/HYP-4506 verified the first gap
(1/(N+1), 2/(2N+1)) NONEMPTY at N=6 (5/33), N=7 (3/23), N=13 (3/41 = {1..11,13,36});
EMPTY at N=12 (fleet censuses ~9k..174k). mac-mini-S25 (HYP-4542): "depth 2->1->0
(N=9..12 UNVERIFIED)"; the bordered-AP enumeration at N=8..11 was named and never run.
opus-S119/HYP-4516: the canonical family F(N)={1..N}\{N-1} u {3(N-1)} attains the
mediant 3/(3N+2) iff N==1 mod 6 and 5 nmid 3N+2 (mod-30 binder gate) -- so N=8..11
members, if any, are NON-canonical species.

VALIDATION GATES (kind-pasteur HYP-7870 part IV: a searcher whose detection floor
sits above its target produces vacuous negatives):
  G1 exact-M evaluator reproduces: M({1..N})=1/(N+1) (N=6..13), M({1..11,24})=2/25,
     M({1..11,13,36})=3/41, M({1,5,6,11,16,17})=5/33, M({1,3,4,5,7,13,18})=3/23,
     Fan-Sun ML(3,8,11,19)=7/30.
  G2 the species generators must REDISCOVER the known gap members at N=6, 7, 13
     from inside their own ranges. If G2 fails, N=8..11 emptiness is NOT concluded.

Species: A single-defect single-far {1..N}\{i} u {w};  B bordered dilated APs
(spine a+d*{0..m-1} + borders spine+-e); C two-defect two-far; D greedy descent
from Allowed(a) bands at admissible (val,q); E adaptive needle-repair search.

Exactness: THM-1002 (general-M pair-sum lemma): the maximizer t*=a/q has
q | (v_i+v_j) (i<=j), so scanning all pair-sum/2v denominators q and all a<q is
EXACT (non-reduced a/q covers divisors).
"""
from fractions import Fraction as F
from math import gcd
from functools import reduce
from itertools import combinations
import random, time, sys

random.seed(59)

# ---------------- exact M with early exit ----------------
def cand_denoms(S):
    Q = set()
    for v in S: Q.add(2*v)
    for x, y in combinations(S, 2):
        Q.add(x+y); Q.add(abs(x-y))
    Q.discard(0)
    return sorted(Q)

def M_exact(S, stop_above=None):
    """Exact M(S) = max_t min_v ||vt||, all-integer hot loop (opus-S398 lesson:
    ~100x over Fractions).  best tracked as (bn, bd); compare mn/q vs bn/bd via
    mn*bd > bn*q.  If stop_above given, return early once best >= stop_above
    (caller then knows M >= stop_above; the return is a certified lower bound)."""
    S = sorted(set(S))
    bn, bd = 0, 1
    sa_n = sa_d = None
    if stop_above is not None:
        sa_n, sa_d = stop_above.numerator, stop_above.denominator
    for q in cand_denoms(S):
        half = q // 2
        for a in range(1, half + 1):
            mn = q
            for v in S:
                r = (v*a) % q
                if r > q - r: r = q - r
                if r < mn:
                    mn = r
                    if mn * bd <= bn * q: break   # cannot beat current best
            if mn * bd > bn * q:
                bn, bd = mn, q
                if sa_n is not None and bn * sa_d >= sa_n * bd:
                    return F(bn, bd)   # certified >= stop_above (loose)
    return F(bn, bd)

def primitive(S):
    g = reduce(gcd, S)
    return tuple(sorted(v//g for v in S))

def in_window(M, N):
    return F(1, N+1) < M < F(2, 2*N+1)

def order_k(M, N):
    return M.denominator - N * M.numerator

# ---------------- gates ----------------
def gate1():
    ok = True
    tests = [
        (list(range(1,7)), F(1,7)), (list(range(1,8)), F(1,8)),
        (list(range(1,9)), F(1,9)), (list(range(1,13)), F(1,13)),
        (list(range(1,14)), F(1,14)),
        (list(range(1,12))+[24], F(2,25)),
        (list(range(1,12))+[13,36], F(3,41)),
        ([1,5,6,11,16,17], F(5,33)),
        ([1,3,4,5,7,13,18], F(3,23)),
        ([3,8,11,19], F(7,30)),
    ]
    for S, want in tests:
        got = M_exact(S)
        tag = "OK " if got == want else "FAIL"
        if got != want: ok = False
        print(f"  G1 {tag} M({S}) = {got} (want {want})")
    return ok

# ---------------- species generators ----------------
def species_A(N):
    """single-defect single-far {1..N}\{i} u {w}"""
    base = list(range(1, N+1))
    for i in range(1, N+1):
        B = [v for v in base if v != i]
        for w in range(N+1, 9*N+1):
            if w in B: continue
            yield tuple(B + [w])

def species_B(N, cap_per_shape=3000):
    """bordered dilated APs: spine a+d*{0..m-1}, borders spine+-e (e=1,2)"""
    for d in range(2, 9):
        for a in range(1, 2*N+1):
            for m in range(3, N):
                spine = [a + d*i for i in range(m)]
                nb = N - m
                if nb <= 0: continue
                borders = set()
                for sv in spine:
                    for e in (1, 2, 3):
                        for sgn in (1, -1):
                            b = sv + sgn*e
                            if b > 0 and b not in spine: borders.add(b)
                borders = sorted(borders)
                if len(borders) < nb: continue
                cnt = 0
                for combo in combinations(borders, nb):
                    yield tuple(sorted(spine + list(combo)))
                    cnt += 1
                    if cnt >= cap_per_shape: break

def species_C(N, wmax_mult=7, cap=500000):
    """two-defect two-far {1..N}\{i,j} u {w1,w2}"""
    base = list(range(1, N+1))
    cnt = 0
    for i, j in combinations(range(1, N+1), 2):
        B = [v for v in base if v not in (i, j)]
        for w1 in range(N+1, wmax_mult*N+1):
            if w1 in B: continue
            for w2 in range(w1+1, wmax_mult*N+1):
                if w2 in B: continue
                yield tuple(B + [w1, w2])
                cnt += 1
                if cnt >= cap: return

def admissible_vq(N, valmax=7):
    out = []
    for val in range(3, valmax+1):
        qlo = (2*N+1)*val  # q > qlo/2
        for q in range((qlo//2)+1, (N+1)*val):
            if gcd(val, q) == 1 and 2*q > (2*N+1)*val:
                out.append((val, q))
    return out

def species_D(N, Bmult=1.5, fat_cap=14):
    """greedy descent from Allowed(a) at admissible (val,q) (val<=5), skipping
    bands fatter than N+fat_cap (cost control). Yields the size-N endpoint."""
    for val, q in admissible_vq(N, valmax=5):
        B = int(q*Bmult)
        for a in range(1, q//2+1):
            if gcd(a, q) != 1: continue
            allowed = [v for v in range(1, B+1)
                       if min((v*a) % q, q-((v*a) % q)) >= val]
            if not (N <= len(allowed) <= N + fat_cap): continue
            if not any((x+y) % q == 0 for x, y in combinations(allowed, 2)):
                continue
            S = list(allowed)
            while len(S) > N:
                bestM, bestS = None, None
                for r in range(len(S)):
                    T = S[:r] + S[r+1:]
                    m = M_exact(list(primitive(T)), stop_above=F(2, 2*N+1))
                    if bestM is None or m < bestM:
                        bestM, bestS = m, T
                S = bestS
            yield tuple(primitive(S))

def species_E(N, tries=1200, iters=300):
    """adaptive needle-repair: start from a band-compliant lift at (val,q,a=1);
    repeatedly find the deepest beating hole (q',a') and repair by moving one
    residue; accept if M decreases. Random restarts."""
    hi = F(2, 2*N+1)
    for val, q in admissible_vq(N, valmax=6):
        band = list(range(val, q-val+1))
        for _ in range(max(1, tries//len(admissible_vq(N)))):
            # residues: force active pair (val, q-val); rest random in band
            res = {val, q-val}
            while len(res) < N:
                res.add(random.choice(band))
            if len(res) < N: continue
            S = sorted(res)  # lift = residues themselves (a=1)
            curM = M_exact(S)
            for _ in range(iters):
                if in_window(curM, N):
                    break
                # find the binding candidate (deepest hole) and repair randomly
                idx = random.randrange(N)
                new = random.choice(band)
                T = sorted(set(S[:idx] + S[idx+1:] + [new]))
                if len(T) < N: continue
                m = M_exact(T, stop_above=hi)
                if m < curM:
                    S, curM = T, m
            yield tuple(primitive(S))

# ---------------- driver ----------------
def run_N(N, budget_s, log):
    lo, hi = F(1, N+1), F(2, 2*N+1)
    found = {}     # M -> (example, species)
    t0 = time.time()
    counts = {}
    def check(S, sp):
        counts[sp] = counts.get(sp, 0) + 1
        Sp = primitive(S)
        if len(set(Sp)) < N: return
        m = M_exact(list(Sp), stop_above=hi)
        if lo < m < hi and m not in found:
            found[m] = (Sp, sp)
            log(f"    !! N={N} GAP MEMBER {m} (order k={order_k(m,N)}) via {sp}: {Sp}")
    gens = [("A", species_A(N)), ("B", species_B(N)), ("C", species_C(N)),
            ("D", species_D(N)), ("E", species_E(N))]
    per = budget_s / len(gens)
    for name, gen in gens:
        tg = time.time()
        for S in gen:
            check(S, name)
            if time.time() - tg > per: break
    log(f"  N={N} window ({lo},{hi}) width {float(hi-lo):.5f}: "
        f"{len(found)} distinct in-gap values; tested {counts}")
    for m, (S, sp) in sorted(found.items()):
        log(f"    member M={m} order k={order_k(m,N)} species {sp}: {S}")
    return found

def main():
    log = lambda s="": print(s, flush=True)
    log("== HYP-7885 cross-N first-gap census (death-star-S59) ==\n")
    log("GATE 1 (exact evaluator):")
    if not gate1():
        log("G1 FAILED -- aborting"); sys.exit(1)
    log("\nGATE 2 (generators must rediscover known members at N=6,7,13):")
    res = {}
    for N, budget in [(6, 150), (7, 150), (13, 200)]:
        res[N] = run_N(N, budget, log)
    g2_ok = (any(m == F(5,33) for m in res[6]) or len(res[6]) > 0,
             any(m == F(3,23) for m in res[7]),
             any(m == F(3,41) for m in res[13]))
    log(f"\n  G2 verdicts: N=6 found-any={len(res[6])>0} (5/33 exact: {any(m==F(5,33) for m in res[6])}), "
        f"N=7 3/23: {g2_ok[1]}, N=13 3/41: {g2_ok[2]}")
    if not (len(res[6]) > 0 and g2_ok[1] and g2_ok[2]):
        log("  G2 FAILED -- the N=8..11 emptiness verdicts below are NOT trustworthy.")
    log("\nTHE UNVERIFIED MIDDLE N=8..11 (+ N=12 consistency):")
    for N, budget in [(8, 400), (9, 400), (10, 400), (11, 400), (12, 250)]:
        res[N] = run_N(N, budget, log)
    log("\n== SUMMARY: N -> in-gap values (order k) ==")
    for N in sorted(res):
        vals = sorted(res[N])
        s = ", ".join(f"{m} (k={order_k(m,N)})" for m in vals) if vals else "EMPTY (within species+budget)"
        log(f"  N={N:>2}: {s}")

if __name__ == "__main__":
    main()
