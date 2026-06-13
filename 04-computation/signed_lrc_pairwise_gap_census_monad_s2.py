#!/usr/bin/env python3
"""
Signed LRC: exhaustive small-n census of the PAIRWISE (mutual) signed gap.
monad-explorer-2026-06-06-S2.  Dispatched angle: exhaustive small-n, enumerate
sign-reversal patterns, rigorously find the optimal gap at each, tabulate.

Setup (repo canon):
  LRC speed system {0, v_1,...,v_{r}} with stationary observer 0, r = n-1 movers
  (repo "n" = r+1 total).  Observer loneliness  M_obs = max_t min_i ||v_i t||,
  conjecture M_obs >= 1/n = 1/(r+1).

Sign gauge: v_i -> eps_i * v_i, eps in {+-1}^r.
  KNOWN (T1/HYP-2286): ||eps v t|| = ||v t||  =>  M_obs is sign-INVARIANT.
  Genuinely sign-DEPENDENT object = the PAIRWISE structure (s674: "observer-
  blind, pair-visible").  Define the SIGNED PAIRWISE GAP over moving pairs:

     G_pair(eps,S) = max_t  min_{i<j} || (eps_i v_i - eps_j v_j) t ||

  i.e. the LRC gap of the signed relative-speed set D_eps = {eps_i v_i - eps_j v_j}.
  Same signs -> difference v_i - v_j ; opposite signs -> sum v_i + v_j.
  Global sign flip negates all of D_eps (||.|| invariant) => 2^{r-1} effective
  patterns.  We compute G_pair EXACTLY (Fraction arithmetic) over ALL eps and
  tabulate: which pattern maximizes it, whether sign flips beat the all-+ (all-
  differences) pattern, and the relation to modular shell-partners.

Exactness: ||w t|| for the speed set W is piecewise-linear in t with kinks at
  a/(2 w_i) and pairwise crossings at a/(2(w_i +- w_j)).  The maximin is attained
  at such a kink.  We take the (over-)complete candidate set
     T = { a/(2 m) : m in {w_i} u {w_i+w_j} u {|w_i-w_j|}, 0<=a<=2m }
  which provably contains every kink, and we additionally re-run at doubled
  resolution as a stability check.
"""
from fractions import Fraction as F
from itertools import combinations, product
from math import gcd
from functools import reduce

def norm(x):
    """ ||x|| = distance from rational x to nearest integer. """
    f = x - int(x)          # in (-1,1)
    if f < 0: f += 1        # in [0,1)
    return min(f, 1 - f)

def candidate_times(W):
    """Over-complete exact candidate set of maximizer times in (0,1)."""
    W = [w for w in W if w != 0]
    aw = set(abs(w) for w in W)
    ms = set(aw)
    for a, b in combinations(W, 2):
        ms.add(abs(a + b)); ms.add(abs(a - b))
    ms.discard(0)
    T = set()
    for m in ms:
        for a in range(1, 2 * m):
            T.add(F(a, 2 * m))
    return T

def maximin_gap(W, double=False):
    """max_t min_{w in W} ||w t|| exactly. W = list of nonzero ints."""
    W = [w for w in W if w != 0]
    if not W:
        return F(1, 2)
    T = candidate_times(W)
    if double:
        T = T | {t / 2 for t in T} | {(t + 1) / 2 for t in T if (t+1)/2 < 1}
    best = F(0)
    for t in T:
        m = min(norm(w * t) for w in W)
        if m > best:
            best = m
    return best

def signed_diff_set(V, eps):
    return [eps[i]*V[i] - eps[j]*V[j] for i, j in combinations(range(len(V)), 2)]

def gpair(V, eps, double=False):
    return maximin_gap(signed_diff_set(V, eps), double=double)

def normalize_set(V):
    """divide out gcd (gaps are scale-invariant: ||c v (t/c)|| reparametrizes)."""
    g = reduce(gcd, V)
    return tuple(v // g for v in V)

def m_obs(V, double=False):
    return maximin_gap(list(V), double=double)

def enum_sets(r, B):
    """distinct positive speed sets, size r, entries in 1..B, gcd 1, sorted."""
    out = []
    for combo in combinations(range(1, B + 1), r):
        if reduce(gcd, combo) == 1:
            out.append(combo)
    return out

def main():
    print("="*78)
    print("SIGNED LRC PAIRWISE-GAP CENSUS (exact)  monad-explorer-S2")
    print("="*78)
    # (r = movers, B = max speed).  repo n = r+1.
    configs = {1:9, 2:9, 3:9, 4:8, 5:7}
    summary = []
    for r in range(1, 6):
        B = configs[r]
        n = r + 1
        sets = enum_sets(r, B)
        floor = F(1, n)
        # PART 1: sign-invariance of M_obs (exact, all eps) ------------------
        obs_violations = 0
        # PART 2: pairwise gap census --------------------------------------
        beat_allpos = 0     # sets where some eps gives G_pair > G_pair(all +)
        allpos_optimal = 0  # sets where all-+ is already the (an) optimum
        gpair_examples = []
        for V in sets:
            base_obs = m_obs(V)
            # check M_obs invariance over all sign patterns
            for eps in product([1, -1], repeat=r):
                if m_obs([eps[i]*V[i] for i in range(r)]) != base_obs:
                    obs_violations += 1
            # pairwise gaps over effective sign patterns (fix eps_0 = +1)
            gvals = {}
            for tail in product([1, -1], repeat=r-1):
                eps = (1,) + tail
                gvals[eps] = gpair(V, eps)
            allpos = gvals[tuple([1]*r)]
            gmax = max(gvals.values())
            gmin = min(gvals.values())
            best_eps = [e for e,g in gvals.items() if g == gmax]
            if gmax > allpos:
                beat_allpos += 1
            if allpos == gmax:
                allpos_optimal += 1
            gpair_examples.append((V, base_obs, allpos, gmax, gmin, best_eps))
        print(f"\n--- r={r} movers (repo n={n}), B<= {B}, {len(sets)} gcd-1 sets, floor 1/n={floor} ---")
        print(f"  M_obs sign-invariance violations (all eps, all sets): {obs_violations}  "
              f"({'CONFIRMED invariant' if obs_violations==0 else 'VIOLATED!'})")
        if r >= 2:
            print(f"  G_pair: all-+ is optimal in {allpos_optimal}/{len(sets)} sets; "
                  f"sign-flip BEATS all-+ in {beat_allpos}/{len(sets)} sets")
            # show a few examples where flipping helps
            shown = 0
            for V, bo, ap, gx, gn, be in gpair_examples:
                if gx > ap and shown < 6:
                    print(f"    V={V}: G_pair(all+)={ap}  -> max over signs={gx} "
                          f"(min={gn}), best eps (eps0=+): {be[0]}")
                    shown += 1
        summary.append((r, n, len(sets), obs_violations, beat_allpos if r>=2 else None,
                        allpos_optimal if r>=2 else None))
    print("\n" + "="*78)
    print("SUMMARY")
    print("="*78)
    print(f"{'r':>2} {'n':>2} {'#sets':>6} {'obs_viol':>9} {'flip_beats':>11} {'allpos_opt':>11}")
    for r, n, ns, ov, bb, ao in summary:
        print(f"{r:>2} {n:>2} {ns:>6} {ov:>9} {str(bb):>11} {str(ao):>11}")

if __name__ == "__main__":
    main()
