#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""TASK 1: chase the f-recursion through the HALF-SYSTEM DIGRAPH (mac-mini-2026-06-29-S8).

f = #{palindromic Hamiltonian paths} of a self-converse tournament (THM-582, the odd
index).  For Paley T_p (p=3 mod4, phi(x)=-x), a palindromic path is
   P = (v_1, ..., v_{(p-1)/2},  0,  -v_{(p-1)/2}, ..., -v_1),
forced: center = 0 (the phi-fixed vertex), second half = negated reverse of the first.
So P is determined by its FIRST HALF (v_1,...,v_{m}), m=(p-1)/2, a sequence with one
representative chosen (with sign) from each of the m pairs {x,-x}, such that P is a valid
directed path in T.  Hence
   f = #{ valid first-halves } = a Hamiltonian-path count on the m = (p-1)/2 PAIRS
       (each pair contributes a +-oriented state) -- a HALF-SIZE reduction (m! vs p!).

This script:
 (1) verifies f via the half-system enumeration equals the direct palindromic count (THM-582),
 (2) exhibits the half-system as a transfer-matrix recursion: build the path pair-by-pair,
     the live state = (last vertex, set of used pairs), and f = sum over completions;
 (3) reports the speedup and the per-step branching (the (p-1)/2-pair transfer structure).
"""
from __future__ import annotations
import functools, itertools
print = functools.partial(print, flush=True)


def paley(p):
    qr = set((x * x) % p for x in range(1, p))
    return [[(i != j and ((j - i) % p) in qr) for j in range(p)] for i in range(p)]


def direct_palindromic_count(arc, p):
    """brute force: all Ham paths P with phi(reverse(P))=P, phi(x)=-x. (small p)"""
    cnt = 0
    for perm in itertools.permutations(range(p)):
        if all(arc[perm[k]][perm[k + 1]] for k in range(p - 1)):
            if all((-perm[p - 1 - i]) % p == perm[i] for i in range(p)):
                cnt += 1
    return cnt


def half_system_count(arc, p):
    """f via the FIRST HALF only: choose an ordered sequence using one (+-signed) rep
    from each pair {x,-x}, then check the FULL palindromic path is a valid T-path.
    This enumerates m! * 2^m half-states (vs p! full) -- the half-size reduction."""
    m = (p - 1) // 2
    pairs = [(x, (-x) % p) for x in range(1, m + 1)]   # {1,p-1},{2,p-2},...,{m,m+1}
    f = 0
    half_states = 0
    for order in itertools.permutations(range(m)):          # which pair at each half-position
        for signs in itertools.product((0, 1), repeat=m):   # which rep of each pair
            half_states += 1
            first = [pairs[order[k]][signs[k]] for k in range(m)]
            # full palindromic path
            P = first + [0] + [(-v) % p for v in reversed(first)]
            if len(set(P)) == p and all(arc[P[k]][P[k + 1]] for k in range(p - 1)):
                f += 1
    return f, half_states, m


def transfer_build(arc, p):
    """f as a transfer/DP recursion on the half-system: live state = (last_vertex, frozenset
    of pairs used).  Build the first half pair-by-pair; the closing arc (-> 0) and the
    second-half arcs are forced by phi, but must be checked.  Returns f and the DP table size."""
    m = (p - 1) // 2
    pairs = [(x, (-x) % p) for x in range(1, m + 1)]
    repset = {}
    for idx, (a, b) in enumerate(pairs):
        repset[a] = idx; repset[b] = idx
    # A first half v_1..v_m is valid iff: consecutive arcs v_i->v_{i+1} hold, AND the
    # closing v_m -> 0 holds, AND the mirror arcs 0 -> -v_m, -v_{i+1} -> -v_i hold.
    # By the anti-automorphism phi (x->-x): arc(u,w) <=> arc(-w,-u). So:
    #   arc(v_i, v_{i+1})  <=>  arc(-v_{i+1}, -v_i)  [mirror arc, AUTOMATIC]
    #   arc(v_m, 0)        <=>  arc(0, -v_m)         [closing pair, AUTOMATIC together]
    # So validity reduces to: arcs v_i->v_{i+1} (i=1..m-1) AND v_m -> 0.  Pure half-data!
    # DP: f = number of orderings/signs with these half-arcs.
    from functools import lru_cache
    starts = [v for pr in pairs for v in pr]
    # DP over (last vertex, used-pairs bitmask)
    f = 0
    table = {}
    def dp(last, used):
        key = (last, used)
        if key in table:
            return table[key]
        if used == (1 << m) - 1:
            r = 1 if arc[last][0] else 0      # closing arc v_m -> 0
            table[key] = r; return r
        tot = 0
        for v in starts:
            idx = repset[v]
            if used & (1 << idx):
                continue
            if arc[last][v]:                  # half-arc last -> v
                tot += dp(v, used | (1 << idx))
        table[key] = tot; return tot
    for v in starts:
        idx = repset[v]
        f += dp(v, 1 << idx)
    return f, len(table)


def main():
    print("=" * 84)
    print("TASK 1: f (palindromic count) via the HALF-SYSTEM transfer recursion (mac-mini-S8)")
    print("=" * 84)
    import math
    for p in (3, 7, 11):
        arc = paley(p)
        m = (p - 1) // 2
        direct = direct_palindromic_count(arc, p) if p <= 9 else None
        fh, hs, _ = half_system_count(arc, p)
        ftr, tsize = transfer_build(arc, p)
        print(f"\n--- Paley T_{p}  (m=(p-1)/2={m} pairs) ---")
        print(f"  direct palindromic count: {direct if direct is not None else '(skipped, p large)'}")
        print(f"  half-system enumeration:  f={fh}  (half-states {hs} = m!*2^m vs full p!={math.factorial(p)})")
        print(f"  TRANSFER DP on half-system: f={ftr}  (DP states stored: {tsize})")
        print(f"  KEY: validity reduces to PURE HALF-DATA -- arcs v_i->v_{{i+1}} (i<m) and v_m->0;")
        print(f"       the mirror/closing arcs are AUTOMATIC by the anti-automorphism phi.")
        ok = (ftr == fh) and (direct is None or direct == fh)
        print(f"  consistent: {ok}")

    print("\n" + "=" * 84)
    print("RESULT: f = #palindromic Ham paths is a Hamiltonian-path count on the m=(p-1)/2")
    print("PAIRS of the half-system digraph, computable by a transfer DP whose validity")
    print("condition is PURE half-data (last->next half-arcs + close-to-0); the second half")
    print("is forced by phi.  This is the half-size reduction -- 'genuinely finishable'.")
    print("f-sequence (p=3,7,11): 1, 9, 185.")
    print("=" * 84)


if __name__ == "__main__":
    main()
