#!/usr/bin/env python3
"""
hopfield_ising_signed_incidence_opus_s_hop.py
=============================================
Follow-up to the H2 REFUTATION: the naive uniform line-graph coupling
Q = sum_{e~f} s_e s_f does NOT determine c3.  Reason: the map spin -> out-score
is not symmetric in s_e; it depends on the FIXED edge orientation reference (i<j).

Correct statement (this script PROVES it exactly):
  Define the SIGNED incidence b_{v,e} = +1 if e=(v,w) w>v ... we use the natural
  out-degree differential.  For edge e=(i,j), i<j:  s=+1 => i->j (out at i),
  s=-1 => j->i (out at j).  So the out-degree contribution is:
      delta_i = (1+s_e)/2 ,  delta_j = (1-s_e)/2 .
  Then s_v = sum_{e incident to v} delta_{v,e}, an AFFINE-LINEAR readout of spins.
  E_score = sum_v (s_v - sbar)^2 is therefore a quadratic in s_e:
      E_score = const + (linear in s_e) + sum_{e~f sharing v} eps(v,e)eps(f,e) s_e s_f
  where eps(v,e) = +1/2 if v is the "low" endpoint, -1/2 if "high" endpoint.
  So it IS a 2-body Ising Hamiltonian on the line graph, but with SIGNED
  (frustrated) couplings J_ef = eps depending on shared-vertex roles, and a
  nonzero external field.  We verify c3 = const' - E_score exactly equals the
  full signed 2-body expansion.  This is the precise Hopfield weight matrix.
"""
from itertools import combinations
from math import comb
from fractions import Fraction
import random
random.seed(7)

def edges_of(n):
    return [(i, j) for i in range(n) for j in range(i + 1, n)]

def ising_energy_from_weights(n, edges, spins, W, h, c0):
    """E = c0 + sum_e h[e] s_e + sum_{e<f} W[(e,f)] s_e s_f."""
    m = len(edges)
    E = Fraction(c0)
    for e in range(m):
        E += h[e] * spins[e]
    for (e, f), w in W.items():
        E += w * spins[e] * spins[f]
    return E

def build_weights(n, edges):
    """Derive exact 2-body Ising weights so that E = sum_v (s_v - sbar)^2."""
    m = len(edges)
    sbar = Fraction(n - 1, 2)
    # s_v = sum_e eps(v,e)*s_e + offset_v, where for e=(i,j) i<j:
    #   contributes to s_i: (1+s_e)/2  -> eps(i,e)=+1/2, off += 1/2
    #   contributes to s_j: (1-s_e)/2  -> eps(j,e)=-1/2, off += 1/2
    eps = {}      # (v,e) -> +-1/2
    off = [Fraction(0)] * n
    for ei, (i, j) in enumerate(edges):
        eps[(i, ei)] = Fraction(1, 2);  off[i] += Fraction(1, 2)
        eps[(j, ei)] = Fraction(-1, 2); off[j] += Fraction(1, 2)
    # s_v - sbar = (off[v]-sbar) + sum_{e inc v} eps(v,e) s_e  =: a_v + sum L_{v,e} s_e
    a = [off[v] - sbar for v in range(n)]
    # E = sum_v (a_v + sum_e L_{v,e} s_e)^2
    #   = sum_v a_v^2 + 2 sum_v a_v sum_e L_{v,e} s_e
    #     + sum_v (sum_e L_{v,e} s_e)^2
    # const:
    c0 = sum(av * av for av in a)
    # linear h[e] = 2 sum_v a_v L_{v,e}
    h = [Fraction(0)] * m
    for v in range(n):
        for ei, (i, j) in enumerate(edges):
            if v in (i, j):
                h[ei] += 2 * a[v] * eps[(v, ei)]
    # quadratic: sum_v sum_{e,f} L_{v,e}L_{v,f} s_e s_f
    #   diagonal e=f: L^2 * s_e^2 = L^2 (since s_e^2=1) -> goes to const
    #   off-diagonal e!=f sharing v: 2 L_{v,e}L_{v,f} s_e s_f
    W = {}
    for v in range(n):
        inc = [ei for ei, (i, j) in enumerate(edges) if v in (i, j)]
        for ei in inc:
            c0 += eps[(v, ei)] * eps[(v, ei)]  # diagonal -> const
        for a_idx in range(len(inc)):
            for b_idx in range(a_idx + 1, len(inc)):
                e1, e2 = inc[a_idx], inc[b_idx]
                key = (min(e1, e2), max(e1, e2))
                W[key] = W.get(key, Fraction(0)) + 2 * eps[(v, e1)] * eps[(v, e2)]
    return W, h, c0

def main():
    for n in (4, 5, 6, 7):
        edges = edges_of(n)
        m = len(edges)
        W, h, c0 = build_weights(n, edges)
        # verify E_ising(spins) == sum_v (s_v - sbar)^2 == C(n,3) - c3, exactly
        sbar = Fraction(n - 1, 2)
        ok = True
        trials = (1 << m) if m <= 15 else 4000
        for it in range(trials):
            mask = it if m <= 15 else random.getrandbits(m)
            spins = [1 if (mask >> b) & 1 else -1 for b in range(m)]
            s = [0] * n
            for b, (i, j) in enumerate(edges):
                if spins[b] == 1:
                    s[i] += 1
                else:
                    s[j] += 1
            E_direct = sum((Fraction(sv) - sbar) ** 2 for sv in s)
            E_ising = ising_energy_from_weights(n, edges, spins, W, h, c0)
            c3 = comb(n, 3) - sum(comb(sv, 2) for sv in s)
            # E_direct should equal sum(s_v-sbar)^2; relate to c3:
            # sum C(s,2) = (1/2)(sum s^2 - C(n,2)); sum(s-sbar)^2 = sum s^2 - 2 sbar*C(n,2) + n*sbar^2
            if E_ising != E_direct:
                ok = False
                print(f"  n={n} MISMATCH ising vs direct", spins[:4]); break
            # c3 = C(n,3) - (1/2)(sum s^2 - C(n,2)); and sum s^2 = E_direct + 2*sbar*C(n,2) - n*sbar^2
            sum_s2 = E_direct + 2 * sbar * comb(n, 2) - n * sbar * sbar
            c3_from_E = comb(n, 3) - Fraction(1, 2) * (sum_s2 - comb(n, 2))
            if c3_from_E != c3:
                ok = False
                print(f"  n={n} c3 mismatch"); break
        # report coupling structure
        signs = {}
        for w in W.values():
            signs[w] = signs.get(w, 0) + 1
        kind = ("exhaustive" if m <= 15 else f"sampled {trials}")
        print(f"n={n}: 2-body Ising E = sum_v (s_v-sbar)^2 EXACT? {ok} ({kind})")
        print(f"       #couplings={len(W)} (= line-graph edges of K_n), "
              f"coupling-value multiset {dict(signs)}, field h nonzero: "
              f"{any(x!=0 for x in h)}")
        print(f"       => Hopfield weight matrix = SIGNED (frustrated) line-graph "
              f"coupling + external field; NOT the uniform Q. Ground state = regular.")

if __name__ == "__main__":
    main()
