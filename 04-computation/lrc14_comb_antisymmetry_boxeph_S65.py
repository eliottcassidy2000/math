#!/usr/bin/env python3
"""
THE COMB ANTISYMMETRY LAW (boxeph-2026-07-17-S65, self-directed target)

From LEM-034(A)'s reflection x -> 1-x (sections sigma -> 6-sigma, R_s ->
R_{6-s}, boundary indices j -> 7e-j i.e. classes c -> -c, signs swap,
owner attribution PRESERVED -- the co-owner lattice at x and 1-x is the
same set, so min(own) is unchanged):

  (A) N_c^{(e)}(R_s) = - N_{(7-c) mod 7}^{(e)}(R_{6-s})   for every owner,
      class, section -- THE COMB ANTISYMMETRY LAW.
  (B) nu-hat spectra:  nu_e^{(s)}(a) = - conj(nu_e^{(6-s)}(a))
      (equivalently  = -nu_e^{(6-s)}(-a); N real).
  (C) sigma_e(s) = -sigma_e(6-s)  (owner imbalances anti-symmetric in s);
      in particular at the CENTRAL SECTION s = 3:  sigma_e(3) = 0 for EVERY
      owner of every cluster -- the owner-imbalance baseline of LEM-030
      VANISHES at s = 3.
  (D) |S^{(s)}(n)| = |S^{(6-s)}(n)| pointwise (endpoint sums reflect), so
      Q_s(w) = Q_{6-s}(w) EXACTLY for every dilation frame w -- the whole
      frame-spectrum theory needs only s <= 3.

Referee: (A)/(B)/(C) exact over all clusters x s x owners x classes;
(D) exact via the THM-880 bilinear form on w-batteries; the s = 3
vanishing checked wherever R_3 is nonempty.
"""

import sys
from math import lcm

sys.path.insert(0, '04-computation')
from lrc14_general_resonance_law_boxeph_S26 import owner_data, nu_hat
from lrc14_hyp6994_resonance_test_boxeph_S25 import Qs_bilinear

CLUSTERS = [([12, 15, 20, 21, 28, 30, 35], "balanced"),
            ([1, 2, 3, 4, 5, 36, 60], "two-owner"),
            ([1, 2, 3, 4, 5, 6, 60], "family60"),
            ([8, 9, 10, 12, 14, 15, 18], "near-AP"),
            ([1, 2, 3, 4, 5, 56, 84], "two-large"),
            ([1, 2, 3, 4, 5, 6, 70], "family70")]

if __name__ == "__main__":
    print("THE COMB ANTISYMMETRY LAW (boxeph S65)")
    print("=" * 78)
    nA = nB = nC = 0
    print("PART A/B/C -- N_c(R_s) = -N_{-c}(R_{6-s}); nu-hat; sigma_e(3) = 0")
    for E, name in CLUSTERS:
        sig3 = None
        for s in range(7):
            P, M, data = owner_data(E, s)
            P2, M2, data2 = owner_data(E, 6 - s)
            if M == 0:
                assert M2 == 0, (name, s)
                continue
            assert set(data) == set(data2), (name, s)
            for e in data:
                N, N2 = data[e]["N"], data2[e]["N"]
                for c in range(7):
                    assert N[c] == -N2[(7 - c) % 7], (name, s, e, c)
                    nA += 1
                for a in range(7):
                    z = nu_hat(data[e], a) + nu_hat(data2[e], a).conjugate()
                    assert abs(z) < 1e-9, (name, s, e, a)
                    nB += 1
                if s == 3:
                    assert sum(N) == 0, (name, e, sum(N))
                    nC += 1
            if s == 3:
                sig3 = [sum(data[e]["N"]) for e in sorted(data)]
        tag = (f"R_3 owners all sigma = 0: {sig3}" if sig3 is not None
               else "R_3 empty")
        print(f"  [{name}] antisymmetry exact; {tag}")
    print(f"  instances: (A) {nA} class checks, (B) {nB} spectrum checks, "
          f"(C) {nC} owners with sigma(3) = 0")
    print()
    print("PART D -- Q_s(w) = Q_{6-s}(w) exactly (THM-880 bilinear referee)")
    for E, name in CLUSTERS[:4]:
        P = 7 * lcm(*E)
        rows = []
        for s in range(4):
            _, M, _ = owner_data(E, s)
            _, M2, _ = owner_data(E, 6 - s)
            if M == 0:
                continue
            worst = 0.0
            for w in (1, 11, 37, 5 * max(E) + 1):
                q1 = Qs_bilinear(E, s, w)
                q2 = Qs_bilinear(E, 6 - s, w)
                worst = max(worst, abs(q1 - q2) / max(abs(q1), 1e-12))
            rows.append(f"s={s}/{6-s}: {worst:.1e}")
        print(f"  [{name}] worst rel gap by pair: " + "  ".join(rows))
    print("=" * 78)
    print("done")
