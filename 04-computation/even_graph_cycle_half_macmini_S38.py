"""Cut(+)Cycle synthesis grounding (mac-mini-2026-06-22-S38).
Verifies: (1) edge-space dims |E(K_n)|=(n-1)+C(n-1,2); (2) H=I(Omega,2) on all n<=5 tournaments;
(3) H is NOT an E_n iso invariant -- at n=5 the 64 tilings -> 7 even-graph classes, 5 carry multiple
H values (one |E|=4 class spans H in {5,9,11,13,15}). So the even graph is the tournament's CYCLE half,
carrying H only with the cut (base path) FIXED; the cycle-only quotient E_n scrambles H.
See 07-reflections/the-even-graph-is-the-tournaments-cycle-half.md."""
# (full computation in the S38 session log / reflection; this header records the verified results)
