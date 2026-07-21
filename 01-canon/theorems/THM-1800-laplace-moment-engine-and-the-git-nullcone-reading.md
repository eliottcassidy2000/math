# THM-1800: the Laplace moment engine (general executed telescoping) and the GIT/moment-map reading of GMC(2)

**Status:** ENGINE: DELIVERED AND VALIDATED (three S186 supports re-derived
from proved identities with all-m-verified recurrences; parametric law
recovered exactly at 6/6 triples; three heavier supports exceed the default
window — honest failure, named fix). GIT READING: definitional/classical
facts verified; the reframing of GMC(2) is a CONJECTURE-equivalent
restatement, not a proof. Tournament bridge: cancellation theorem verified;
naive involution instability at n >= 4 recorded. NO completion claim.
**Author:** boxeph-2026-07-20-S188 (HYP-8600)
**Owner:** "work more on the representation theory of binary forms and how
it relates to tournaments, which are in/transitivity itself. work laplace
moment engine."

## 1. The Laplace moment engine

General-support implementation of the S187 executed telescoping
(04-computation/laplace_moment_engine_boxeph_S188.py): input any
P = sum c_i Z^{a_i}W^{b_i}; the engine auto-generates the three PROVED
relation families for Lambda = sum c_i w^{nu_i} u^{d_i}
  (E) N_m(j,k) = sum_i c_i N_{m-1}(j+d_i, k+nu_i)
  (U) 0 = j N_m(j,k) + m sum_i c_i d_i N_{m-1}(j+d_i, k+nu_i)
  (W) N_m(j,k+1) = k N_m(j,k-1) + m sum_i c_i nu_i N_{m-1}(j+d_i, k+nu_i-1)
and performs the left-kernel elimination, interpolation in m, and exact
verification. RESULTS:
- H1 (aZ+b+cW), H2 (aZ+bW^2+c), H4 (aZ+bZW+cW): recurrences derived and
  verified (H4 reproduces the S187 certificate recurrence exactly).
- PARAMETRIC LAW RECOVERED: at 6/6 rational triples the engine returns
  q[-1] = -b and q[-2] = -4ac(m-1) for the linear span — the exact law
  mu_m = b mu_{m-1} + 4ac(m-1) mu_{m-2} with coefficients VISIBLY
  polynomial in (a,b,c) (the parametric-certificate demonstration).
- HONEST FAILURES: spans (2,0,-1), (1,-3), and the hijack-flavor
  (0,4,-1)-type exceeded the default window (levels = 4 caps the order at
  3; higher spans need order > 3 and larger (J,K)-windows; the hijack
  support burned 264s finding nothing). Named fix: sparse elimination +
  adaptive level growth. Engineering item: this is the seed of the
  mod_rank-style reusable library (repo engineering mandate).
- S187r discipline inherited: displayed recurrences are per-triple;
  degenerate strata (vanishing raw top coefficient) need re-stratification.

## 2. The GIT/moment-map reading of GMC(2)  [verified facts + a reframing]

(i) E IS THE FISCHER/BARGMANN PAIRING: E[Z^aW^b] = delta_ab 2^a a! is the
Bargmann inner product (variance 2) — the invariant-theory pairing on
binary-form coefficient space. (ii) ONE-SIDEDNESS = HILBERT-MUMFORD
INSTABILITY for the hyperbolic torus: P is one-sided iff the 1-PS
(tZ, t^{-1}W) (or its reflection) drives P to 0 as t -> 0. (iii) Hence
GMC(2) restates as: **the analytic moment-nullcone equals the GIT nullcone
of the U(1)-hyperbolic action** — a Kempf-Ness-shaped statement. Hilbert's
nullcone is the common zero locus of invariants; here the moments E[P^m]
play the invariants, and the conjecture says the analytic and algebraic
notions of instability coincide. The repo's "nullcone" terminology was the
right word all along; the S180 gauge-invariance trap (critical values are
gauge-invariant) is the Kempf-Ness stationarity in disguise.

## 3. The tournament shadow: intransitivity is what the pairing kills

(R1, verified n = 3,4,5): the Vandermonde expands as the SIGNED TOURNAMENT
SUM prod_{i<j}(x_i - x_j) = sum_T (-1)^{rev(T)} x^{score(T)}, and the
surviving monomials are EXACTLY the transitive tournaments (permutations,
unit signs): ALL INTRANSITIVITY CANCELS. Micro-finding (honest): the naive
canonical involution (reverse the lex-first 3-cycle) is NOT stable at
n >= 4 (reversal can change the canonical triangle); the global
cancellation is verified regardless; the stable pairing needs a finer
canonical rule — same care-class as Redei's path-reversal. MECHANISM-SHAPE
IDENTITY: Redei (h odd: transitive survives mod 2, the rest pairs off) and
the Vandermonde (transitive survives exactly, the rest cancels) are the
same statement-shape: an invariant pairing kills exactly the intransitive
part. In the GIT reading, "one-sided = triangular = transitive-like" and
GMC(2) says the moment pairing kills exactly the non-triangular part —
the owner's "tournaments are in/transitivity itself" made precise: the
three cancellation theorems (Redei parity, Vandermonde expansion, the
conjectural GMC(2)) are one family.
Classical dictionary recorded: bracket monomials of the symbolic method =
weighted tournaments (vertex degrees = form degree); Grassmann-Pluecker
straightening = matching flips; Cayley-Sylvester dims (computed, R4) =
partitions in a box at half-weight vs score sequences = partitions under
the staircase (A000571 = 1,2,4,9,22 verified) — counting-frame KINSHIP
recorded as a lead, no identity claimed.

## 4. Files

- 04-computation/laplace_moment_engine_boxeph_S188.py + frozen .out.
- 04-computation/binary_forms_tournaments_boxeph_S188.py + frozen .out.
- Reflection: 07-reflections/gmc2-as-a-git-nullcone-and-the-three-cancellation-theorems.md
