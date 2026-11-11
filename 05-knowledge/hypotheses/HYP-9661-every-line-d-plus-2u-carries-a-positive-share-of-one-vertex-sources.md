---
id: HYP-9661
title: "Sources on lines. (A) Line density: for every odd N >= 129 the line D + 2U = N carries at least N/64 one-vertex sources (D = H(B), U = U_B(y)), and the share of its N/6 lattice points that are sources is bounded below and increases with N; near the diagonal U = D the share tends to 1. (B) On 23 <= h <= 73 the covered share of the window [h, 3h/2] at offset h has correlation 0.68 with the number of strong tournament classes with H = h, and the four offsets of that range with two strong classes (35, 39, 47, 63) are among its five poorest; on all attained 3 <= h <= 73 the correlation is 0.41. (C) Remark on the head sheet. (D) Sources as a local-global problem: for F1 the by-value local densities are positive in every class at each fixed odd prime and squarefree odd modulus and in every odd class modulo every power of 2 (THM-5128); that no congruence excludes any odd number, 7 and 21 included, was KNOWN; the thin-orbit analogy (kind-pasteur S11) lacks a congruence gap uniform in the modulus and a minor-arc bound, and has a counter-indication (reciprocity obstructions)."
status: >
  OPEN (opus-2026-11-10-S91, 2026-11-11). Independent audit 2026-11-11 (seven auditors; record
  05-knowledge/results/s91_f1_local_densities_20261110/audit/AUDIT_REPORT.md): PASS WITH CORRECTIONS, all applied;
  MISTAKE-761. Not recomputed by the audit, because added while it ran: the residue-class statistics of (A)
  (line_density_classes.out) and the head census of (C) (head_values_check.out).
  (A) is FINITE-EXACT as a lower bound for odd N in [129, 16384) (THM-5129 (F): a count on a finite range at rooted
  order <= 10; four pre-registered tests, all CONFIRMED; the audit notes that L1, L2 and L4 are one-sided,
  that L3 already held in sample at order 9, and that all were passed with wide margins) and a conjecture beyond; it implies Chan-Pak completeness with the finite base below 129 and
  is a quantitative form of HYP-9453 (1), (2) and of HYP-9454 (1) (klein: "the expected number of pairs on a line
  of height v grows linearly in v"), which are first; klein also has the line h + 2B = v over all rooted pairs
  (HYP-9454 (6); THM-4865 (e), (f)). Added here: the count r(N), its share and its saturation near the diagonal.
  No proof route is known. No window [h, 3h/2] is full at the offsets computed (h <= 73 except 1 and 5; 129; 257:
  THM-5129 (D), (B')); whether windows fill up for large h is OPEN.
  (B) NUMERICAL: Pearson correlation 0.68 on 23 <= h <= 73 (26 attained offsets), 0.66 on 15 <= h <= 73, 0.41 on
  all attained offsets 3 <= h <= 73 (source_supply.out, section (6)); the multiplicities are those of THM-4884 (2).
  (C) is a remark: "head value" alone is vacuous (PROVED); the question at non-sink roots (are the heads with
  D <= N/3 the tail values?) is answered for odd N < 382,688 (FINITE-EXACT, valid at every order of the base;
  added by the verification of the corrections) and OPEN beyond; "qualifying path" is used in a restricted sense
  (the word was defined first in odd_c_seed_symmetry_20261110/EXTENSIONS.md).
  (D) is a reading of THM-5128 for the F1 family only. That no congruence excludes any odd number, 7 and 21
  included, is KNOWN (Session 25; Session193 modular_sync Theorem 1; Session206 exact-modular-image); THM-5128 adds
  the size of the local densities of F1 at fixed odd primes, squarefree odd moduli and powers of 2 (odd prime-power
  moduli are NUMERICAL only; uniformity in the modulus is OPEN: it needs at least HYP-9660 (B)). The ANALOGY with thin orbits is
  kind-pasteur S11's (seam_gadget_thickness_spectrum_20261104.md; S82 note section 7); added here: the comparison
  of ingredients and the Apollonian counter-example. The theorems on thin orbits are CITED from the literature: the
  session did not re-read them; the audit (citations) read the Zaremba paper of Bourgain and Kontorovich, quotes
  Haag-Kertzer-Rickards-Stange by theorem number, and confirmed the Apollonian paper of Bourgain and Kontorovich
  and the paper of Rickards and Stange at abstract level only; the verification of the corrections checked the
  journal references of Haag-Kertzer-Rickards-Stange (arXiv:2307.02749) and of Rickards-Stange (arXiv:2401.01860)
  on their arXiv listings, and Theorem 1.6, Conjecture 1.5 and Corollary 2.6 of the former in the arXiv text.
session: opus-2026-11-10-S91
date: 2026-11-11
source: 05-knowledge/results/s91_f1_local_densities_20261110/README.md
scripts:
  - 04-computation/experiments/s91_f1_local_densities_20261110/rooted_spectrum/line_density_order10.py (+ line_density_order9.out, line_density_order10.out), line_density_classes.py (+ .out); 04-computation/experiments/s91_f1_local_densities_20261110/PREREGISTRATION.md (second registration; the file is in the experiment's root directory)
  - 04-computation/experiments/s91_f1_local_densities_20261110/rooted_spectrum/source_supply.py (section (6) of its output has the strong-class counts and the correlations of (B); it also writes slope_lists_le73.txt), source_supply_levels78.py, head_parity_check.py, head_values_check.py, census10.py (+ .out)
  - 05-knowledge/results/s91_f1_local_densities_20261110/audit/AUDIT_REPORT.md (independent audit); 05-knowledge/results/s91_f1_local_densities_20261110/audit/audit-5129/s11_supply257.out (the single-offset shares of (A), recomputed with the auditor's own code)
  - 05-knowledge/results/s91_f1_local_densities_20261110/audit/verify-corrections/HYP-9661/ (verification of the corrections; v1_supply_checks, v2_class_law_and_f1, v3_heads_all_orders, each .py with its .out; v3 is the receipt of the answer in (C))
related:
  - 01-canon/theorems/THM-5129-* (finite supply, permanent omissions, line counts), THM-5128-* (local densities of F1)
  - 05-knowledge/hypotheses/HYP-9453-*, HYP-9454-* (klein; near-join sweep, lines 3H + 2 sigma), HYP-9658-* (seeds of the 2-adic lines)
  - 01-canon/theorems/THM-4865-*, THM-5126-*, THM-5130-*, THM-1370-h-spectrum-omits-7-21-all-n.md; 05-knowledge/results/session261_creative_chanpak_20261111/, session262_interface_supply_20261111/
  - 05-knowledge/hypotheses/HYP-9659-* (tilted local limit law; the bivariate version is what (A) asks for)
  - 01-canon/theorems/THM-4884-* (2) (strong multiplicities m(h)), THM-4744-* (inventories), THM-4810-* (atlas), THM-4886-* (A) (modulus-uniform bound for the one-apex count)
  - 05-knowledge/hypotheses/HYP-9660-* (B) (uniform gap, OPEN), HYP-9480-* (D) (misses of one root change), HYP-9657-* (B)
  - 05-knowledge/results/odd_c_seed_symmetry_20261110/EXTENSIONS.md, kps_s22_sign_flip_towers_20261110/README.md section 4, session257_permutohedron_squares_20261111/history/sources/INCOMING_EXTENSION_CONTENT.md (earlier uses of qualifying and k-vertex extension)
  - 05-knowledge/results/seam_gadget_thickness_spectrum_20261104.md (kind-pasteur S11, thin-orbit analogy), one_apex_join_representation_mass_pi2_over_12_20261106.md sections 6, 7, 9 (S82 note)
  - 05-knowledge/results/session25_spectrum_wiring_20261021/session25-spectrum/PROOF.md, session193_triangulated_proof_20261105/modular_sync/PROOF.md (Theorem 1), session206/exact-modular-image/DYADIC_SURJECTIVITY.md at the repository root (no congruence obstruction, KNOWN)
---

# HYP-9661 - sources on lines

Notation of THM-5129. A one-vertex source is `(D, U) = (H(B), U_B(y))`; it gives the head value `D - E + U` and the
closed values `D + 2^t U`, `t >= 1`. Sinks appended to `B` double `U`, so the closed values of all sources are the
numbers `D + 2U`. `U >= D`, hence `D <= N/3` on the line `D + 2U = N`, which has about `N/6` lattice points with
`D` odd. `r(N)` is the number of sources on the line.

## (A) Line density (CONJECTURE; implies Chan-Pak)

1. `r(N) >= N/64` for every odd `N >= 129`.
2. The share `r(N) / (N/6)`, averaged over a dyadic block, is increasing in the block and bounded below by `0.27`.
3. For every `eps > 0` the share of the lattice points with `D <= U <= (1 + eps) D` that are sources tends to 1.

Data (THM-5129 (F); line_density_order10.out; counts on a finite range, not a proof that every line carries a
source): shares `0.296, 0.334, 0.381, 0.464, 0.518, 0.588, 0.616` on `[2^k, 2^(k+1))`, `k = 7..13`,
from rooted tournaments of order `<= 10`; minimum of `r` over a block `5, 7, 24, 63, 146, 352, 794`; top decile of
`3D/N` on `[4096, 8192)`: `0.94` (odd `U`), `0.84` (even `U`). The two-step family misses `7, 11, 13, 21, 47, 67`
and no other odd `N >= 3` below 382,688 (THM-5129 (C), THM-4865 (e)). At a single offset the supply is known exactly
at every order (THM-5129 (B), (B')): the share of the integers of `[h, 9h/8]` that are slopes is `0.67, 0.80,
0.78, 0.82, 0.85` at `h = 17, 33, 65, 129, 257`, of `[h, 2h)` it is `0.41, 0.67, 0.57, 0.66, 0.65`, and of
`[2h, 4h)` `0.29, 0.44, 0.38, 0.50, 0.48` (no output of the session's scripts prints these shares; receipt: the
audit's recomputation with its own code, audit/audit-5129/s11_supply257.out).

Residue classes (exploratory, not pre-registered, not recomputed by the audit; `line_density_classes.out`). On
`[2^11, 2^13)` the share of lattice points that are sources is `0.578` for `3 | N` and `0.557` otherwise (ratio
`1.038`; under a Poisson reading of the share, `-ln(1 - share)`, the intensities are in the ratio `1.06`); at
`5, 7, 11, 13` the ratios are `1.004, 1.001, 0.999, 1.001`, within `0.004` of 1. For comparison, in F1 the class
factor of the distinguished class of a sheet (`3 | N` on the tail sheets) is `(1 - 1/9 + 1/3)/(1 - 1/9) = 1.375`
times that of the other classes in the limit (THM-5128 (D3)). Among the sources of a line with `p | N`, the share
with `p | D` is `0.357, 0.201, 0.143, 0.092, 0.078` at `p = 3, 5, 7, 11, 13` (`1/p = 0.333, 0.200, 0.143, 0.091,
0.077`): in this statistic the closure of the source set under `(D, U) -> (m D, m U)` for attained `m` leaves
little trace. This is a saturating 0/1 statistic at order `<= 10` (a lattice point is a source or it is not; mean
share `0.564`), which compresses class contrasts; it is not a representation count, so `1.038` and `1.375` are not
directly comparable. Whether the multiplicity-weighted supply is equidistributed in residue classes is not tested.

Why a count; what is known about windows. For a fixed offset `h` the slopes are finitely many odd parts times
powers of two (THM-5129 (A2)): at most `|Delta(h)|` per octave (THM-5129 (A3)). This bound excludes a full window
at eight of the attained offsets `h <= 73` (`3, 13, 35, 39, 43, 47, 63, 67`; slope_lists_le73.txt) but not at the
others, nor at `129` and `257`: at `h = 33, 65, 129, 257` there are `34, 66, 177, 359` odd parts against
`17, 33, 65, 129` integers in `[h, 3h/2]`. What is computed (THM-5129 (D), (B')): no window `[h, 3h/2]` is full
at the offsets computed (`h <= 73` except `1` and `5`; `129`; `257`); the covered share rises from about `0.50`
(`37 <= h <= 73`) to `0.75` and `0.72` (`129`, `257`); whether windows fill up for large `h` is OPEN. (A) does
not use windows: on a line the candidate bases are the odd `D <= N/3`, about `N/6` of them, and the conjecture is
that a positive share of them have `(N - D)/2` among their slopes (between `0.30` and `0.62` in the tested range,
shares above). Linear growth of the number of sources on a line for all `N` is (A) and, in klein's form,
HYP-9454 (1). The statement to prove is a bivariate local law for `(H(B), U_B(y))`, the rooted analogue of
HYP-9659.

Refutation: an odd `N >= 129` with fewer than `N/64` sources among all rooted tournaments, or a dyadic block whose
saturated share falls below the previous one.

## (B) Supply and fibre size (NUMERICAL)

Let `m(h)` be the number of strong classes with `H = h`, over all orders (THM-4884 (2), from the inventories of
THM-4744; recounted for `h <= 73` in source_supply.out, section (6)). On `23 <= h <= 73` (26 attained offsets) the
covered share of `[h, 3h/2]` (THM-5129 (D)) and `m(h)` have Pearson correlation `0.68`; it is `0.66` on
`15 <= h <= 73` (29 offsets) and `0.41` on all attained offsets `3 <= h <= 73` (34 offsets). The value depends on
the lower cutoff, which was not registered in advance. The four offsets of `[23, 73]` with `m(h) = 2`
(`35, 39, 47, 63`: shares `0.278, 0.350, 0.208, 0.156`) are among the five poorest; the third poorest, `67`
(`9/34 = 0.265`), has `m = 6`. The other offsets with `m(h) = 2` are `9, 17, 19`, with shares `0.60, 0.44, 0.50`.
So on this range a small fibre goes with a poor window, but `m(h)` does not determine the share. By THM-5129 (A2),
(A3) a covered share bounded below at offset `h` needs `|Delta(h)| >= c h`, and `Delta(h)` is built from the
strong tournaments whose path counts divide `h`; in that sense a uniform supply at offset `h` is a statement
about multiplicities in the strong spectrum (HYP-9453 (3): at least two strong realisers for `h >= 15`;
THM-4884 (2): `m(h) >= 6` for every attained `h` with `63 < h < 625`; the atlas of THM-4810). Question: is
`|Delta(h)| >= c h` for all large attained `h`?

## (C) The head sheet (REMARK; one PROVED triviality, one question, answered below 382,688)

Every attained number is a head value for a trivial reason: if `y` is a sink of `B` then `E = U = D` and
`H(B + z) = D`. So "head value" is not a constraint; the informative heads are at roots that are not sinks, where
`H(B + z) = D + 2 U_even` with `U_even >= 1` (THM-5129 (E)). The four numbers that no tail reaches are such heads:
`11 = 5 + 2 * 3`, `13 = 5 + 2 * 4`, `47 = 17 + 2 * 15`, `67 = 33 + 2 * 17`. The head is a different sheet from
the tails: it counts only the qualifying paths with an even number of vertices, and in the join family F1 it is
the sheet `s = 0`, whose favoured residue class is `-E = -1` instead of `0` (THM-5128 (C)).

Wording. "Qualifying path" is used here in this session's restricted sense: a directed path `R` of `B` starting
at `y`, weighted by `H(B - V(R))`; the total weight is `U`. The words qualifying and k-vertex extension were
defined first in odd_c_seed_symmetry_20261110/EXTENSIONS.md (qualifying permutations of the `n + k` vertices of
an extension by `k` new vertices with arbitrary attachments; an exact compiler of degree at most `2k`; its mean
is THM-5126 (C)) and refined in Session257 (INCOMING_EXTENSION_CONTENT.md); kind-pasteur S22
(kps_s22_sign_flip_towers_20261110/README.md, section 4) first read the qualifying count as the cover count (the
slope) of a one-vertex extension. The usage here is the restriction to the one-vertex interface with a
transitive tail: an interpretation, not the general object.

FINITE-EXACT (rooted census of order `<= 9`, odd `N < 9153`; head_values_check.out; added during the audit and
not recomputed by it): every odd `N >= 3` other than 7 and 21 is a head at a root that is not a sink; with the
restriction `D <= N/3` the numbers missed are again exactly `7, 11, 13, 21, 47, 67`. So the four attained
exceptions need a base carrying more than a third of their paths (`5/11, 5/13, 17/47, 33/67` for the heads
above), and below 9153 no head with `D <= N/3` reaches a number that the tails do not. Question: is that an
identity of the two value sets, or a coincidence of the small range?

Answer for odd `3 <= N < 382,688` (FINITE-EXACT, valid at every order of the base; added by the verification of
the corrections, receipt audit/verify-corrections/HYP-9661/v3_heads_all_orders.out): an identity. (i) A tail value
is by definition a head: `B => TT_t => z` is `(B => TT_t) + z`; the base `B => TT_t` has the same `D`, its root
`y` is not a sink, and `D <= N/3` because `U >= D`. (ii) The four attained exceptions have `D > N/3` in every
head at a root that is not a sink, at every order of the base: there `B + z = B' => S`, with `S` the strong
component of `z` (`z` and the vertices of `B` reachable from `y`) and `B'` the rest, so
`H(B + z) = H(B') H(S)` and `D = H(B') H(S - z)`; the four numbers are prime and `H(S) >= 3`, hence
`H(B') = 1` and `S` is a strong tournament with that many paths, of order at most 8 (`f(9) = 75`); over the
strong classes of order `<= 8` and their vertices `z` of out-degree 1 the bases `S - z` have `D = 5` (11 and
13), `17, 23` (47) and `29, 31, 33` (67), and the census of order `<= 9` has no other `D` for these four heads.
(iii) 7 and 21 are not path counts (THM-1370), and every other odd number of that range is a tail value
(THM-5129 (C), THM-4865 (e)). So on that range the heads with `D <= N/3` at a root that is not a sink and the
tail values are the same numbers. For all `N` the equality is OPEN; it would follow from (A).

## (D) Sources as a local-global problem

**PROVED (THM-5128 (A)).** Let `p` be a prime and let two independent states in `(Z/p)^r` be glued by a
nondegenerate bilinear form (for F1 and odd `p`: `M = x_1 x_2 + 2^s y_1 y_2`, `r = 2`, and `N = M - [s = 0]`).
Then the classes `c != 0` of the form have equal probabilities to first order, and the excess of the class `0`
over their mean is `pi_1 + pi_2 - p^(-r)`: the masses of the states that vanish modulo `p`, less the generic
`p^(-r)`. The remainder is of second order in the distances of the two state laws from the uniform law. If both
states are equidistributed the law is `1/p + p^(-r) ([c = 0] - 1/p)`, positive in every class; this is immediate
and not a finding.

**KNOWN (no congruence obstruction).** That no congruence excludes any odd number, 7 and 21 included, is KNOWN:
the one-apex family alone attains every residue modulo every odd modulus and every odd residue modulo every even
modulus (Session 25; Session193 modular_sync Theorem 1; Session206 exact-modular-image).

**What THM-5128 adds** is the size of the local densities of F1 at fixed odd primes, squarefree odd moduli and
powers of 2. By value, on each sheet, the class factors at an odd prime tend to `1 - 1/p^2 + [c = c_sigma]/p`
(`c_sigma = -1` on the sheet `s = 0`, `0` on the sheets `s >= 1`), with a correction of order `1/ln X`; for
squarefree odd `m` the limiting density of a class is the product of these; the odd classes modulo every power of
2 are equidistributed (THM-5128 (C), (D3), (E); PROVED). All these densities are positive. So the residue-class
side is determined at each fixed odd prime and squarefree odd modulus (and every power of 2), with a power saving
uniform in the region but not in the modulus; uniformity in the modulus is OPEN (it needs at least HYP-9660 (B));
prime-power moduli
(`p^k`, `p` odd, `k >= 2`) are numerical only (THM-5128 (D3)).

**Consequence, for the F1 source family only.** F1 is a source family in the sense of README section 4.3 (a
one-apex source `(x_1, y_1)` closed by a one-apex module `(x_2, y_2)`): offset `D = x_1 x_2`, slope `U = y_1 y_2`,
head correction `E = 1`; the sheet `s = 0` is the head `D - E + U`, the sheets `s >= 1` are the tails
`D + 2^s U`. For this family, at a fixed odd prime, the classes other than the distinguished one are equal up to
a power-saving error, the only first-order local structure is the mass of the triples whose smaller one-apex
state is `p`-divisible (`p | x` and `p | y`), and it distinguishes the head (class `-E = -1`) from the tails
(class `0`). For a general source family with readout `D + 2^t U` nothing of this kind is proved; the bivariate
law of (A) would be the corresponding statement. The missing values `7 = N(1, 1, 1) = 3 + 2^2` and
`21 = N(1, 2, 2) = 5 + 2^4` (here `N(d, r, s) = 1 + 2^s + d 2^(r+s)` are the grid coordinates of THM-5129 (B),
not the F1 triples of THM-5128) are invisible to every such test, as they are to every congruence (KNOWN,
above). In the one-vertex interface with offset `2^s + 1` they sit where the row `d = 1` is not supplied: the
offset 3 has no dyadic slope at all, and the cone `5 + 2^m` is supplied only from `m = 6`, so that its head is
`13, 21, 37` (THM-5129 (B)). These unsupplied parts also contain attained numbers (`11, 19` at offset 3;
`13, 37` at offset 5), reached through other interfaces, so an unsupplied head is where the supply of one
interface stops, not an obstruction.

**ANALOGY (orbit and bilinear form match; the congruence gap is proved at each fixed modulus only; the minor-arc
bound is missing).** The thin-orbit / Zaremba analogy for the path spectrum is kind-pasteur S11's
(seam_gadget_thickness_spectrum_20261104.md), taken up in the S82 note
(one_apex_join_representation_mass_pi2_over_12_20261106.md, section 7); added here: the comparison of ingredients
and the Apollonian counter-example. The one-apex states are the orbit of `(psi, h, floor(h/2), 1)` under two
integer matrices (one per binary digit), F1 is a bilinear form on two orbits, and the twisted operators `T_k` have
a spectral gap (`rho(T_k) < 2`) at every fixed odd modulus (THM-5128 (B1)); a gap uniform in the modulus is
HYP-9660 (B), OPEN (a modulus-uniform bound is in canon only for the one-apex count on full blocks, THM-4886 (A)).
Orbit, bilinear form and a congruence gap uniform in the modulus are inputs of the two density-one theorems of
Bourgain and Kontorovich for thin orbits (this description of the inputs is the audit's reading of the
introduction of the Zaremba paper): almost every integer is the denominator of a reduced fraction with partial
quotients bounded by an absolute constant, for which 50 suffices (Ann. of Math. 180 (2014) 137-196), and for a
primitive integral Apollonian packing all but `O(N^(1 - eta))` of the admissible integers up to `N` are
curvatures (Invent. Math. 196 (2014) 589-650). In both the mean multiplicity is a power of `N`; for F1 it is
`0.748 ln N` (THM-5128 (D1)). So here orbit and bilinear form are present, the congruence gap is proved for each
fixed modulus only, and the minor-arc bound is missing. Such a bound (a second-moment or collision bound) would
give at most a density-one statement for F1. Completeness needs a pointwise lower bound, which no bound of this
kind gives (S82 audit M3; THM-4886, NOT CLAIMED) and which F1 alone cannot satisfy: it has misses at every scale
tested (HYP-9480 (D); to `2^35` in the S82 note, section 6). Counter-indication, and the reason no density
argument can be read as evidence for completeness: the local-global conjecture is false for every primitive
Apollonian packing outside two types (Haag, Kertzer, Rickards, Stange, Ann. of Math. 200 (2024) 749-770):
quadratic and quartic reciprocity exclude infinite families of admissible curvatures (no curvature of the packing
`(-3, 5, 8, 8)` is a square), which no congruence sees, and beyond these families only finitely many sporadic
misses are conjectured. The same happens for semigroup orbits in `SL(2, Z)`, the setting closest to the digit
semigroup: Rickards and Stange (Duke Math. J. 174 (2025) 3197-3244) refute the local-global conjecture of
Bourgain and Kontorovich for continued fractions over a finite alphabet (families whose denominators are never
squares). On sub-families the repository has obstructions of congruence type (a self-join
`x^2 + 2^s y^2 - [s = 0]` is never `7 (mod 8)`) and a statistical pattern that no single congruence explains
(among the 30 primes `4q - 1 < 3000` with `q` prime, F1 misses eight, `7, 43, 67, 163, ...`: integers with
`N (N + 1)` free of small odd primes carry the least F1 mass; HYP-9657 (B), NUMERICAL, a statement about F1 that
does not single out 7 and 21). Question: is `{7, 21}` a sporadic set in this sense, or the trace of a
reciprocity-type obstruction on the whole spectrum that happens to be finite?
