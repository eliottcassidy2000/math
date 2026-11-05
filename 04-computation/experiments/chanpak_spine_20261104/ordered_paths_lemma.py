"""Ordered-path lemma test (mac-mini-2026-11-04-chanpak, after audits C and D).

For a tournament X with h Hamiltonian paths and a cut (A, B) of a new vertex v (v -> A, B -> v) such that some arc goes
from A to B (equivalently v lies on an odd cycle of X + v), let q = number of Hamiltonian paths of X of the form
(all of B)(all of A).  PROVED: H(X + v) >= 2h - q (every path of X has a slot for v - between a consecutive B, A pair,
in front if it starts in A, at the back if it ends in B - at least two slots unless the path is ordered, and the slot
map is injective).  TESTED here: q <= (h - 1)/2, which would give H(X + v) >= (3h + 1)/2, i.e. HYP-9469 (1) in the
sharp form k >= (h + 1)/4, with THM-4876 (7) as the case |A| = 1 (q = e_a).
q = sum over y in B, x in A with y -> x of e_y(X[B]) s_x(X[A]) from the subset tables; H(X + v) from ext_spectrum.
Usage: python3 ordered_paths_lemma.py NMAX
"""
import sys, numpy as np
from fractions import Fraction
from ext_spectrum_true import tournaments, end_start_tables, ext_spectrum
NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 8
for n in range(4, NMAX + 1):
    N = 1 << n; full = N - 1
    cases = viol = 0; worst = Fraction(0); wc = None; big = Fraction(0); bigc = None; slack_min = None; tight = 0
    for adj in tournaments(n):
        e, s = end_start_tables(n, adj)
        h = int(e[full].sum())
        Hv = ext_spectrum(n, adj)
        for Am in range(1, full):
            Bm = full ^ Am
            A = [x for x in range(n) if Am >> x & 1]; B = [y for y in range(n) if Bm >> y & 1]
            if not any(adj[x] >> y & 1 for x in A for y in B):
                continue
            q = sum(int(e[Bm][y]) * int(s[Am][x]) for y in B for x in A if adj[y] >> x & 1)
            cases += 1
            H = int(Hv[Am]); slack = H - (2 * h - q)
            assert slack >= 0, (n, adj, Am, h, q, H)
            if slack == 0: tight += 1
            if slack_min is None or slack < slack_min: slack_min = slack
            if 2 * q > h - 1:
                viol += 1
                if viol <= 3: print("VIOLATION", n, adj, A, h, q)
            r = Fraction(q, h)
            if r > worst: worst = r; wc = (h, q, len(A))
            if len(A) >= 2 and len(B) >= 2 and r > big: big = r; bigc = (h, q, len(A))
    print(f"n={n}: cuts with an A->B arc {cases}; violations of q <= (h-1)/2: {viol}; max q/h = {worst} at (h,q,|A|)={wc}; "
          f"max with |A|,|B|>=2: {big} at {bigc}; min H(X+v)-(2h-q) = {slack_min} (attained at {tight} cuts)", flush=True)
