# Proof and computation audit

Status: **VERIFIED exact replay and self-review**, not an independent-agent
audit. The session was performed by one agent. The proof conclusions are
separate from the finite universes below.

## Proof gates

1. **Graph loss:** W separates into discarded/discarded and selected/discarded
   wedges. Greedy domination pays each discarded square once. Cauchy-Schwarz
   then gives W>=2D^2/n; K_(r,r) proves the coefficient cannot improve uniformly.
2. **Shared budget:** the heavy/light inequality uses disjoint assigned edge
   sets. Its completion of a square charges eta^2 once per light group.
3. **All-size stability:** the incoming per-L tail maximum is at an endpoint
   by convexity, and both endpoints decrease. Every base inequality at
   n=276,R=138 is rational. R*eta=4 is an equality at the endpoint; the
   contradiction remains strict because aL<4D-RB. The negative-coefficient
   case is handled before multiplying by the upper bound B<=4D/R.
4. **Equality:** D=B=0 implies a matching uncertainty graph and exclusively
   transitive supported binary completions. The inherited H=1 scalar theorem,
   not the graph inequalities alone, supplies fair weights and maximal pairing.
5. **FRACTRAN:** all fractions are positive and reduced, all states positive
   integers. The observation includes original rule labels; HALT has no state
   residue. The clock stops after exactly r legal steps and does not add a
   further halt query. These conventions are essential to P2 and P4.
6. **Minimality:** any period is tested at each realizable word's explicit
   witness D_w. P4 minimizes a residue modulus. P6 separately minimizes an
   arbitrary finite observer for m=1 and the reciprocal family, under clock
   decrement. Neither asserts a finite clock-free universal carrier.
7. **Transfer boundaries:** real absolute-value stability does not reduce
   modulo an integer. Hamiltonicity of a square-sum graph does not make that
   graph a tournament. Completing orders needs the allowed-edge guard.
   Equal H does not preserve the determinant gap under fair-pair padding.

## Exact finite universes and independent implementations

`verify.py` uses only Python's standard library and active exception checks,
not assertions removed by -O. Normal and optimized outputs must agree.

* All 3^10=59049 weighted graphs on five labelled vertices with integer
  weights in {0,1,2}; all 2^15=32768 simple graphs on six vertices; 48 seeded
  sparse weighted stress cases at orders 7,11,25,64. Total **91865** graphs.
  Explicit incident-pair summation is compared with degree squares. The
  global loss bound and disjoint shared budgets are both checked.
* The exact n>=276 rational endpoints, intermediate n>=500 constants,
  scalar sign cases, and entire finite L ranges at R=138,139,150,250.
  The proof's monotonicity supplies the infinite range. A failing idealized
  n=274 endpoint is retained as a rejection control, not a counterexample
  to the mathematical conjecture.
* All 133 programs of length 0,1,2 over the eleven reduced fractions a/b
  with 1<=a,b<=4, at horizons 1..3 and terminal moduli 1..4: **1596** cases.
  The guarded compiler is compared with independent literal execution over
  two full periods (**334546** integer states). A KMP period calculation
  independently verifies minimality. There are 5408 explicit prefix-witness
  checks; the largest tested minimal modulus is 6912.
* **2688** reciprocal-family parameter cases, a,b=2..8 coprime,
  horizons 1..6, moduli 1..16, compare the general guarded compiler with P5.
  P6 is tested by literal trajectories and clock-decrement transitions on
  complete residue periods for (2,3),(3,5),(5,2), at horizons 1..6.
  One-rule 1/p controls force precision m*p^r for p=2,3,5,11.
* **729** exact triangle grid cases (p,a,b in {0,1/8,...,1}) check the joint
  reserve and its fair-edge identity. The equal-variance/different-cycle
  hostile is retained.
* Native subset Hamilton-path DP and rational Gaussian determinants check
  ordinal padding on two weighted/binary three-vertex sources through order
  seven. Independent square-sum subset DP counts all paths up to reversal
  for n=1..25, reproducing 1,1,1,3,10 at 15,16,17,23,25.

## Controls, corrections, and provenance

The first complete normal replay and its optimized replay passed. Adding the
n>=276 refinement, reciprocal minimal carrier, and square-sum controls
required fresh replays. No failed mathematical claim was promoted in this
packet. The original n>=500 derivation is retained as a valid intermediate
argument, clearly superseded in range.

The published THM-4852 prime-support-to-units-digit implication was false;
its precise repair and witness are logged in `01-canon/MISTAKES.md`, with
the MISTAKE-697 correction lineage. The historical audit itself is unchanged.
Its squarefree random programs did not cover the squareful assertion.

The initial mathematical baseline was `291b4d9757842c02aa451beb65acfade2848609b`.
During this session, incoming `400dcd6295` and `9c7eb0b13e` supplied Session139's
n>=748/n>=526 arguments, all-modulus contextual clocks, and weighted zeta
counterexamples. These were read before finalizing n>=276. The per-edge
geometric tail and coupled-mass method are explicitly credited. The new
degree cancellation and retained B budget integrate them rather than
claiming their rediscovery as novel. Source pins record the exact inputs.

The first push met concurrent mainline updates and was rejected normally.
The subsequent clean rebase incorporated `b68db2afda` (Session139 integral
and regular spectral boundaries) and `e9dd62edb6` (n>=466 and one-step
arithmetic clocks). The overlapping baseline method and one-step formula
are credited in the proof notes. The general guarded-horizon formula and
reciprocal minimal observer retain their stated distinct scopes.

Primary source reading: Conway 1987 sections 1,8,11; the authors' quasi-Riemann
TeX with the source blob named in `SYNTHESIS.md`; Gerbicz's OEIS contributor
record. The forum proof returned 403, so the all-size square-sum theorem is
cited with that limitation. No assertion about RH's resolution is made.

## Remaining boundaries

Sharp stability at 10..275 is OPEN. P4 can enumerate exponentially many
legal words. P6's smaller state space is clocked and has no terminal-residue
observable. Minimal arbitrary strong-atom carriers, realization of every
modular tournament context state, and full integer Hamilton-path spectra
remain separate open problems. A scalar or a modulus is not a replacement
for an actual complemented relation, source, joint ports, or transport phase.
