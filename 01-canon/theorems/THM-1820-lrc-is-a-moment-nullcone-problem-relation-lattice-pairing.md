# THM-1820: LRC is a moment-nullcone problem — the relation-lattice pairing

**Status:** DELIVERED (bridge identity verified; dictionary with three
experimental confirmations; wild proposals filed to the backlog). The
identity is classical (Weyl/Poisson) — the contribution is the DICTIONARY
and its measured confirmations. NO claims on LRC itself.
**Author:** boxeph-2026-07-20-S189 (HYP-8615)
**Owner:** "see if the LRC is a moment nullcone problem, see if reframes
like this give you inspiration, propose them wildly."
**Builds on:** death-star-S67 (HYP-8515: the S_m-level GMC/LRC kinship);
this goes one level deeper (the lattice level).

## 1. The answer: YES — the dictionary

**Bridge identity (Weyl/Poisson; verified: 0.073333 exact-intervals vs
0.073346 truncated-lattice, K = 40):** for integer speeds v and periodic f,
  int_0^1 prod_j f_j(v_j t) dt = SUM over the RELATION LATTICE
  {k in Z^n : k . v = 0} of prod_j fhat_j(k_j).

| GMC(2) | LRC |
|---|---|
| charge-balanced words (sum charges = 0) | relation vectors k . v = 0 |
| factorial weights 2^a a! | sinc weights sin(2 pi k delta)/(pi k) |
| moment functional E[P^m] | good/bad-set measures |G_delta| (incl-excl of lattice sums) |
| charge-imbalanced => moments vanish trivially | Q-independent speeds => only k = 0: |G| = (1-2delta)^n > 0: loneliness FREE |
| one-sided = Hilbert-Mumford unstable (S188) | trivial relation lattice = Kronecker-dense orbit |
| nullcone question: WHO achieves all-zero (answer: one-sided) | covering question: NOBODY covers below 1/(n+1) (nullcone EMPTY) |
| hyperbolic gauge (invariance) | speed scaling v -> lambda v (delta* invariant) |
| finite moment test / mod-q rungs | the repo's LRCMod ladders = the relation lattice tested mod q |
| conjugate/reality stacks (S183) | mirror pairs under t -> -t (measured: ALL component deaths pair) |
| mu_g symmetry stacks | the tight family's terminal HIGHER stack (measured: multiplicity 4 at delta*) |

Both problems: BOUND THE CONSPIRACY of a relation/charge lattice against a
product pairing. GMC classifies the conspirators (= one-sided); LRC asserts
no conspiracy fully covers below threshold.

## 2. Measured confirmations

(B2) TIGHTNESS = RELATION RICHNESS: n = 3 battery (speeds <= 9): the ONLY
delta* = 1/4 (threshold) families are (1,2,3) and scalings — with N_R = 8
small relations, vs mean 2.81 for loose families (delta* > 0.30). The
extremal LRC configurations are the relation-RICHEST — the LRC face of the
fleet's fresh "transitivity is the deepest nullcone point" (THM-1810 era).

(B3) THE delta-LADDER: sweeping delta, good-set components die at fold
events. MEASURED: every family stacks deaths in MIRROR PAIRS (the t -> -t
reality symmetry: LRC's face of the S183 reality-stacks — stacking is the
BASELINE here, not the exception); the tight family (1,2,3,4) is
distinguished by FEWER, BIGGER events — a terminal multiplicity-4 stacked
death AT delta* — the higher-symmetry (mu_g-like) stack. Prediction
partially confirmed and sharpened: tight = high-multiplicity terminal
stack, generic = mirror-pairs only.

## 3. Wild proposals (filed to INVESTIGATION-BACKLOG; one line each)

W1. q-REDEI: weight Hamiltonian paths by q^inv; is h_q(T) nonvanishing at
    roots of unity (a q-analog of the parity nullcone-emptiness)? Cheap
    census feasible.
W2. LANDAU = HILBERT-MUMFORD: score-sequence realizability (Landau's
    inequalities = dominance) as the numerical criterion for a torus
    action whose weight polytope is the score polytope. Candidate exact
    statement; if true, tournament realizability IS a GIT stability test.
W3. THE delta-RESURGENCE PROGRAM: treat delta as the resurgence parameter;
    fold events = arcs; transplant the S183-S187 machinery (odd-sector
    vectors, monodromy defects, sum rules) to the good-set ladder; the
    mirror-pairing is the reality-stack; ask what the ANALOG of the
    quantization dichotomy says about covering transitions.
W4. LRC FINITE MOMENT TEST: is |G_delta|(as a function of rational delta)
    P-recursive/piecewise-polynomial with an engine-style certificate?
    (Endpoints rational; the ladder has finitely many folds: YES-shaped;
    the content is the explicit rung bound vs the mod-q ladders.)
W5. THE {7,21} INSTABILITY GAP: read the h-monoid as a value-semigroup of
    a GIT quotient; forbidden values as non-attained stability weights.
    (Wildest; filed as tangent-grade.)

## 4. Files

04-computation/lrc_moment_nullcone_boxeph_S189.py + frozen .out (B1-B3).
