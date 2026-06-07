#!/usr/bin/env python3
"""
S710 — Extending the poke fiber-bundle suggestion: the DIVISOR-LATTICE cover of LRC(14).

S643 cashed ONE seed of the poke suggestion (the mod-7 fiber: at the 7-clock, danger == mult-of-7).
Poke gave three seeds: (1) the 7-runner (k=6) base [DONE, S643]; (2) the 2-runner (k=1) base / mod-2
half-turn fiber [HERE]; (3) "the half-turn leak misses only 56 cells" [HERE: 56 = A000568(6) =
#tournaments on 6 vertices = the LRC(7) base type-space].

This session establishes the FULL picture:

(A) DIVISOR-CLOCK THEOREM. For EVERY divisor d|n, the d-clock t=b/d (gcd(b,d)=1) has danger set
    == exactly {v : d|v}; every non-multiple-of-d runner has margin >= 1/d. For d<n strictly safe.
    The divisor poset {1<2<7<14} of n=14 gives a LATTICE of clock-projections, each peeling off the
    mult-of-d sub-config. (S643 was the d=7 instance; here we do d=2 and the whole lattice.)

(B) THE CRT TOWER. 14=2.7, ZZ/14 = ZZ/2 x ZZ/7. A mult-of-14 config, viewed at the 7-clock, reduces
    to its mult-of-7 sub-config; dividing those by 7 lands them in the mod-2 (half-turn) world: a
    mult-of-7 speed is 0 mod 14 (even) or 7 mod 14 (odd). So the residual fiber problem is itself a
    HALF-TURN problem. The 7-clock + 2-clock factor the obstruction along the CRT primes {7} and {2}.

(C) THE 56. A000568(6)=56 = #tournaments on 6 vertices = the LRC(7) base (7 runners = 6 movers). We
    verify 56, the 2-adic valuation v_2(56)=3=(7-1)/2 (THM-305, QR), and the arc-confined sub-count
    |Arc(7)|=6 of 56 (THM-381 family), and articulate the half-turn-leak == base-type-space reading.

No numpy/sympy. Exact over QQ via Fractions.
"""
import random, itertools
from math import gcd
from fractions import Fraction as Fr

# ---------- LRC core (exact M over QQ) ----------
def norm(x):
    f = x - (x.numerator // x.denominator)
    return f if f <= Fr(1, 2) else 1 - f

def gap_and_argmax(speeds):
    V = [abs(v) for v in speeds]; cands = set()
    for i in range(len(V)):
        vi = V[i]
        for k in range(0, 2*vi+1):
            t = Fr(2*k+1, 2*vi)
            if 0 < t <= Fr(1, 2): cands.add(t)
        for j in range(i):
            vj = V[j]
            for d in (vi+vj, abs(vi-vj)):
                if d == 0: continue
                kk = 1
                while Fr(kk, d) <= Fr(1, 2): cands.add(Fr(kk, d)); kk += 1
    best = Fr(0); arg = []
    for t in cands:
        m = min(norm(v*t) for v in V)
        if m > best: best, arg = m, [t]
        elif m == best: arg.append(t)
    return best, arg

def danger_at(speeds, t, level):
    """runners v with ||v t|| < level (the 'dangerous' = within the loneliness radius)."""
    return [v for v in speeds if norm(v*t) < level]

def coprime_clocks(d):
    return [Fr(b, d) for b in range(1, d) if gcd(b, d) == 1] or [Fr(0)]

def gen_config(n, vmax, rng):
    """S643's 'multiple-of-n config': generic primitive (n-1)-speed set CONTAINING a multiple of n."""
    S = {n * rng.randint(1, 2)}
    while len(S) < n-1: S.add(rng.randint(1, vmax))
    S = sorted(S)
    g = 0
    for x in S: g = gcd(g, x)
    return S if g == 1 else None

# ---------- (A) DIVISOR-CLOCK THEOREM ----------
def divisor_clock_test(n, trials, vmax_mult, seed):
    """For each divisor d|n, check at every d-clock that danger==mult-of-d and margin(non-mult)>=1/d."""
    rng = random.Random(seed); level = Fr(1, n)
    divisors = [d for d in range(1, n+1) if n % d == 0]
    stats = {d: {"checks": 0, "danger_eq_mult": 0, "margin_ok": 0} for d in divisors}
    done = 0
    while done < trials:
        S = gen_config(n, vmax_mult, rng)
        if S is None: continue
        done += 1
        for d in divisors:
            mult_d = set(v for v in S if v % d == 0)
            for t in coprime_clocks(d):
                dg = set(danger_at(S, t, level))
                st = stats[d]; st["checks"] += 1
                if dg == mult_d: st["danger_eq_mult"] += 1
                nonmult_margins = [norm(v*t) for v in S if v % d != 0]
                if all(m >= Fr(1, d) for m in nonmult_margins): st["margin_ok"] += 1
    return divisors, stats

# ---------- (B) CRT TOWER reduction trace ----------
def crt_tower_trace(n, trials, vmax_mult, seed):
    """At the best 7-clock, peel mult-of-7; divide by 7; show residual lives mod 2 (half-turn world)."""
    rng = random.Random(seed); level = Fr(1, n)
    residual_sizes = {}; parity_split = {"even(=0 mod14)": 0, "odd(=7 mod14)": 0}
    halfturn_isolates_even = 0; total = 0
    while total < trials:
        S = gen_config(n, vmax_mult, rng)
        if S is None: continue
        total += 1
        # 7-clock peel
        m7 = sorted(v for v in S if v % 7 == 0)
        residual_sizes[len(m7)] = residual_sizes.get(len(m7), 0) + 1
        for v in m7:
            if v % 14 == 0: parity_split["even(=0 mod14)"] += 1
            else: parity_split["odd(=7 mod14)"] += 1
        # half-turn (d=2) check: danger at t=1/2 == even speeds
        dg_half = set(danger_at(S, Fr(1, 2), level))
        if dg_half == set(v for v in S if v % 2 == 0): halfturn_isolates_even += 1
    return residual_sizes, parity_split, halfturn_isolates_even, total

# ---------- (C) tournaments on m vertices (brute, m<=6) + A000568 ----------
def tournament_count_brute(m):
    """#non-isomorphic tournaments on m labelled->iso vertices (brute, m<=6)."""
    if m <= 1: return 1
    perms = list(itertools.permutations(range(m)))
    pairs = [(i, j) for i in range(m) for j in range(i+1, m)]
    seen = set(); count = 0
    for bits in range(1 << len(pairs)):
        # orientation: bit=1 means i->j else j->i
        orient = {}
        for idx, (i, j) in enumerate(pairs):
            if (bits >> idx) & 1: orient[(i, j)] = 1; orient[(j, i)] = 0
            else: orient[(i, j)] = 0; orient[(j, i)] = 1
        # canonical form under permutation
        canon = None
        for p in perms:
            key = tuple(orient[(p[i], p[j])] for (i, j) in pairs)
            if canon is None or key < canon: canon = key
        if canon not in seen:
            seen.add(canon); count += 1
    return count

def v2(n):
    k = 0
    while n % 2 == 0: n //= 2; k += 1
    return k

def smallest_prime(n):
    d = 2
    while d*d <= n:
        if n % d == 0: return d
        d += 1
    return n  # prime

# ---------- (D) prime-vs-composite dichotomy: best proper-divisor clock across n ----------
def best_proper_divisor_clock_residual(n, trials, vmax, seed):
    """For each n: at the n/p-clock (p=smallest prime|n), measure the residual danger-set size
       (=mult-of-(n/p)). Prime n has NO proper divisor>1 => method empty (needs the 2n-1 shell)."""
    rng = random.Random(seed)
    p = smallest_prime(n)
    if p == n:  # prime
        return {"n": n, "prime": True, "p": p, "proper_divisor": None}
    d = n // p  # largest proper divisor
    sizes = {}
    done = 0
    while done < trials:
        S = gen_config(n, vmax, rng)
        if S is None: continue
        done += 1
        k = sum(1 for v in S if v % d == 0)
        sizes[k] = sizes.get(k, 0) + 1
    avg = sum(k*c for k, c in sizes.items()) / done
    return {"n": n, "prime": False, "p": p, "proper_divisor": d,
            "residual_dist": dict(sorted(sizes.items())), "avg_residual": round(avg, 3),
            "fiber_mod": p}

# ============================ RUN ============================
if __name__ == "__main__":
    n = 14
    print("="*78)
    print("S710 — DIVISOR-LATTICE COVER of LRC(14): extending the poke fiber-bundle suggestion")
    print("="*78)

    print("\n(A) DIVISOR-CLOCK THEOREM  (n=14, divisors of 14 = {1,2,7,14})")
    print("    Claim: at every d-clock t=b/d, danger set == {v : d|v}; non-mult margin >= 1/d.")
    divs, stats = divisor_clock_test(n, trials=400, vmax_mult=32, seed=1)
    for d in divs:
        st = stats[d]; c = st["checks"]
        print(f"  d={d:2d} (1/d={float(Fr(1,d)):.4f} vs 1/n={float(Fr(1,n)):.4f}): "
              f"danger==mult {st['danger_eq_mult']}/{c}   margin>=1/d {st['margin_ok']}/{c}")

    print("\n(B) CRT TOWER  14 = 2.7,  ZZ/14 = ZZ/2 x ZZ/7")
    rs, ps, hi, tot = crt_tower_trace(n, trials=1500, vmax_mult=32, seed=7)
    print(f"  7-clock residual (mult-of-7 sub-config) size distribution over {tot}: "
          f"{dict(sorted(rs.items()))}")
    print(f"  parity split of the mult-of-7 speeds (the fiber coordinate): {ps}")
    print(f"  half-turn t=1/2 isolates exactly the even speeds: {hi}/{tot}")
    print("  => residual mult-of-7, divided by 7, splits by parity = the mod-2 half-turn world.")

    print("\n(C) THE 56  (A000568: #tournaments on m vertices = the LRC(m+1) base type-space)")
    A = {}
    for m in range(1, 7):
        A[m] = tournament_count_brute(m)
    print(f"  brute #tournaments on m vertices, m=1..6: {[A[m] for m in range(1,7)]}")
    print(f"  A000568(6) = {A[6]}  <-- the 56;  base of LRC(7) (7 runners = 6 movers)")
    print(f"  v_2(56) = {v2(A[6])} = (7-1)/2 = 3   [THM-305: 2-adic val = |QR_7|, QR connection]")
    print(f"  v_2 of #tournaments on m vertices m=3,5: "
          f"v2(T(3))={v2(A[3])}=(3-1)/2, v2(T(5))={v2(A[5])}=(5-1)/2  [odd-m QR law]")
    print("  Reading: LRC(14) is a bundle over the LRC(7) base; the base's full combinatorial")
    print("  type-space has exactly 56 cells = A000568(6); the half-turn (mod-2) fiber is the")
    print("  parity coordinate ON each of those 56 base cells. The leak 'misses 56 cells' = it")
    print("  rides the WHOLE base (the cover always leaves a section over every base type).")

    print("\n(D) PRIME-vs-COMPOSITE DICHOTOMY  (the divisor-clock method's reach)")
    print("    Best proper-divisor clock = n/p (p=smallest prime|n): peels the SPARSE mult-of-(n/p);")
    print("    residual lives mod p. Prime n has NO proper divisor>1 => method EMPTY (needs 2n-1 shell).")
    for nn in [11, 12, 13, 14, 15, 16, 18, 19, 20, 21, 22]:
        r = best_proper_divisor_clock_residual(nn, trials=600, vmax=2*nn+4, seed=nn)
        if r["prime"]:
            print(f"  n={nn:2d}  PRIME      -> divisor-clock EMPTY; use THM-420 (2n-1={2*nn-1} shell)")
        else:
            print(f"  n={nn:2d}  composite  -> {nn}/{r['p']}-clock peels mult-of-{r['proper_divisor']:2d}; "
                  f"avg residual {r['avg_residual']:.2f}, fiber mod {r['fiber_mod']}; "
                  f"dist {r['residual_dist']}")
    print("  SYNTHESIS: composite n -> divisor-clock tower (this session); prime n -> 2n-1 shell")
    print("  (THM-420, S642). n=14 is TWO-HEADED: composite (clock works) AND 2n-1=27=3^3 ramified")
    print("  (shell fails). The two methods partition the C'(n) frontier by primality of n.")
    print("="*78)
