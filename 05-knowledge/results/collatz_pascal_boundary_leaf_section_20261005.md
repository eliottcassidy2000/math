# Exchangeable prices are the Pascal boundary: counter-only Collatz flows, the leaf-section tower, and what 11 and 22 are

2026-10-05, opus session `opus-2026-10-05-S5` (pascal-boundary).
**Status.** PROVED (elementary, with CITED classical inputs): Theorem A
(classification of counter-only incoming flows as Hausdorff moment arrays),
Corollaries A1-A4, Lemma B (no positive supersolution is comparable to a
power law), Proposition C (Riesz decomposition of summable supersolutions),
Proposition D (leaf-section tower). VERIFIED numerical, float64 with stated
margins: the adversarial-lift (Bellman) mass bounds of section 6.
FINITE-EXACT: all listed checks. HEURISTIC, typed and quarantined: the
dictionary of section 4 beyond its proved rows. **OPEN:** positivity at every
positive odd integer, equivalently universal Collatz. Nothing here is a canon
promotion or a literature-priority claim.

**Scope correction, 2026-10-05 (Codex source-floor audit).** Theorem A
classifies exact split arrays; it includes `E,W`, but the discounted `D`
satisfies a weighted split instead. Its support equivalence does not exclude
future proofs of support. The inherited observer obstruction concerns its
specified observations, not arbitrary finite-dimensional representations.
Proposition C concerns summable potentials and harmonic mass; it neither
classifies a Martin boundary nor requires the injection to charge every
vertex. A2's endpoint integrability condition is the general criterion;
a positive-power density condition is only a beta-family specialization.
The current text below incorporates these corrections. Exact minimal
witnesses: [scope audit](collatz_pascal_scope_audit_20261005.out).

## 0. The prompt, decoded against the repository

The pasted formula `E(n)=6(M+1)!(T+1)!/(M+T+3)!` with "total mass at most
11", "incoming-flow inequality", and "a second mixture with strict payment and
mass at most 3" is B1, (5) and the double mixture `D` of today's Codex notes
[critical beta weights](collatz_critical_beta_weights_20261005.md) and
[mixed refuel weights](collatz_mixed_refuel_weights_20261005.md). The sentence
"every odd integer coprime to 3 has a predecessor less than 22 times its size"
is the three-divisible section of
[adaptive mixture flow, section 6](collatz_adaptive_mixture_flow_20261005.md),
proved as Proposition S1 of [fusion helpers, section 3](collatz_fusion_helpers_20261005.md):
for a unit `v` the least predecessor divisible by 3 is `rho(v)=(2^(a_0(v))v-1)/3`
with `a_0(v)<=6`, hence `rho(v)<=(64v-1)/3<22v`.

The owner asks four things. Each gets one answer here, and each answer is typed.

1. **Positivity everywhere.** It stays OPEN. Changing an interior price prior
   does not change the support of the canonical rooted flow (Theorem A,
   Corollary A3). The specified finite-lookahead observers and power-law
   comparability face separate obstructions. Every summable supersolution
   is the potential of a measure on rooted sources plus a cycle-constant
   harmonic part (Proposition C). These results do not exclude every possible
   source encoding or a future proof that the canonical support is full.
2. **`11*2=22` and `{2,3,11}`.** The two numbers have different mechanisms:
   `11=3!*H_3=6(1+1/2+1/3)` comes from the beta(2,2) prior and the period
   three of sibling depths modulo three (`ord_9(4)=3`); `22=ceil(64/3)` comes
   from `ord_9(2)=6`. The product relation is NUMEROLOGY: the beta(1,2) prior
   gives `16/3`, and twice that is not a section constant. The structural
   statements are the tower of Proposition D (`22` is its first rung
   `2^(2*3^j)/3`) and the fact that the 11 is a bound that is not attained
   (section 6: the beta(2,2) mass lies in a bracket well below 11).
3. **Category theory fused with representation theory, continuous-discrete,
   graphical classification.** The honest proved content is Theorem A: the
   exact two-way split `w(L,K)=w(L+1,K)+w(L,K+1)` is the harmonicity condition
   on the Pascal graph, the Pascal graph is the Bratteli diagram of the
   gauge-invariant CAR algebra and the fusion graph of `1+chi` for `U(1)`,
   and its harmonic functions (traces) are exactly the probability measures on
   the continuous boundary `[0,1]` (Hausdorff, de Finetti; Bratteli; Vershik-Kerov).
   The Codex weights `E`, `W` are such traces pulled back along the counter
   map; `D` is a separately discounted array. Section 4 types the analogy where the
   theorems stop.
4. **Langlands.** Only the abelian (`GL(1)`) layer exists here: counter-only
   weights are mixtures of positive characters of the free monoid on the
   two letters {base, sibling} that factor through its abelianization, and Theorem A is
   their Fourier/Laplace classification. The arithmetic side at that layer is
   the discrete logarithm modulo `3^(j+1)` (the Artin coordinate of
   [HYP-9174](../hypotheses/HYP-9174-x2x3-lonely-spectrum-is-the-artin-coordinate-of-2-3.md)),
   which is exactly what the section exponent `a_0(v)` is. Word-dependent
   constructions retain information lost by counters. The inherited
   finite-lookahead theorem excludes its specified observers; it does not
   exclude arbitrary finite-dimensional exact encodings. This is a typed
   analogy, not a Langlands correspondence; see section 4.

## 1. Inheritance and board

- Closest proved mechanism: C3 and C8 of
  [three bits and critical flow](collatz_three_bits_critical_flow_20261005.md)
  (universal Collatz iff a strictly positive summable `v` with `Kv<=v`; the
  canonical critical chain flow with support exactly the rooted component),
  then P2-P5 of the adaptive note (beta(1,2) array, exact fibre payment, the
  injection measure `lambda` on rooted multiples of three with mass one).
- Canonical hostile: `53` and `113` have the same counters `(1,3)` and
  different words (Codex). Theorem A forgets the word; this example shows
  non-injectivity, not impossibility of proving a predicate constant on a fibre.
- Corrected near miss: reading `22=2*11` as structure. Section 7.
- Least-used sidecar: the order of a letter sequence, i.e. the word itself.
  Every object below either forgets it (Theorem A) or is defined by it
  (Proposition C, Lemma B).
- Board: **split identity / price measure `mu` / counter map / leaf measure
  `lambda` / section exponent `a_j` / ascending shadow / Bellman lift bound.**

**Concurrent commits read before publishing (`598879da4`, `2b38c4084`,
`5e291040b`, all Codex, same afternoon).**
[Inverse predecessor sections](inverse_predecessor_sections_20261005.md)
proves the full three-class palette `kappa_c(n)` (my unit table is its rows
`c=1,2` merged), the sharp suprema `64/3-1/(3n)` along `n=18t+1` and
`16/3-1/(3n)` along `n=18t+7`, the exact size-induction boundary (a smaller
three-divisible predecessor exists iff `n=5 mod 9`), a strictly increasing
unit inverse ray `n<R(n)<32n/3`, and the obstruction that every fixed size
cap loses all relative `W`-mass on the family `S^K(1)`; it states, as this
note does, that no identification with the primes 2, 3, 11 or an automorphic
object is used. [Level-22 oldform bridge](collatz_level22_oldform_bridge_20261005.md)
takes the owner's level literally (`S_2(Gamma_0(22))` is two oldforms of the
conductor-11 newform, with the explicit `U_2` matrix) and concludes that the
predecessor constant and the modular level are two separate proved roles of
22. [Representation positivity](collatz_representation_positivity_20261005.md)
and [Fourier atom positivity](collatz_fourier_atom_positivity_20261005.md)
show that faithful finite affine representations and continuous encodings of
`lambda` do not supply individual atoms. The updated adaptive note (P6, P7)
improves the `W`-mass bound to below `41/10` with a checked cost-eight kernel
and builds a root-deleted singular mixture `Z` (prior `(1-r)/r dr`) with
nonroot mass below `23/8`. Everything below is consistent with these: the
size-cap obstruction and Lemma B are two manifestations of a lost coordinate
(bounded arithmetic ratio does not control flow ratio). `Z` needs the separate
sigma-finite-prior argument in P7, and section 6 numerically sharpens `41/10` to `3.89` with a lower
bound `3.19`.

The Codex notes already observed that their arrays are "completely monotone in
precisely this elementary sense" (adaptive note, (8)). The step taken here is
the converse: every nonnegative split array is such a mixture, so the family
they constructed is the whole family, and its support question is the same for
every member.

## 2. The snippet, exactly

Write `U(n)=oddpart(3n+1)`, `S(n)=4n+1`, bases `b` with `v_2(3b+1) in {1,2}`,
and for a rooted `n>1` with ROOT word `(a_1,...,a_tau)` the counters
`L=tau-1` and `K=sum_i floor((a_i-1)/2)` (adaptive note (2)). Then

    E(n) = 6 (K+1)!(L+1)!/(L+K+3)! = 6 B(K+2, L+2),
    W(n) = 2 K!(L+1)!/(L+K+2)!     = 2 B(K+1, L+2),

with `(M,T)=(K,L)` in the pasted notation. `E` is the probability that a
Polya urn started with two balls of each colour produces a given colour
sequence with `K` draws of one colour and `L` of the other; `W` is the same
with one and two starting balls. Both satisfy the exact split
`w(L,K)=w(L+1,K)+w(L,K+1)`. The bound 11 is `6 int_0^1 (1+r+r^2) dr =
6(1+1/2+1/3) = 3!*H_3`; the three terms are the three sibling residues modulo
three, the 6 is `1/B(2,2)`. With the sharper proved fixed-price bound
`2+2r-r^2` of the adaptive note (P1) the same integral gives `17/2`, as the
Codex JSON already records (`beta_prior_bounds`). The root ray carries `E`-mass
exactly 3 and `W`-mass exactly 2; the rooted multiples of three carry exactly
2 and 1. The "second mixture with strict payment and mass at most 3" is the
double mixture `D` of the mixed note, section 3, with incoming ratio
`(T+1)/(T+3)`. All of this is re-verified exactly in S1 of the script, together
with the controls `E(27)=1/21296186450` at counters `(40,8)`, `53,113 ->
(1,3)`, `111 -> (23,7)`, `155 -> (29,8)`.

## 3. Theorem A: counter-only flows are the Pascal boundary

**Theorem A (classification; PROVED, with Hausdorff's moment theorem CITED).**
Let `w: N x N -> [0,infinity)` satisfy `w(0,0)=1` and
`w(L,K)=w(L+1,K)+w(L,K+1)` for all `L,K>=0`. Then there is a unique Borel
probability measure `mu` on `[0,1]` with

    w(L,K) = int_0^1 r^K (1-r)^L dmu(r)   for all L,K.

Conversely every such integral satisfies the hypotheses.

*Proof.* Put `m_K=w(0,K)`. The split says `w(L+1,K)=w(L,K)-w(L,K+1)`, so by
induction `w(L,.) = (-Delta)^L m` where `Delta h(K)=h(K+1)-h(K)`. Nonnegativity
of all `w(L,K)` is exactly complete monotonicity of `m` with `m_0=1`.
Hausdorff's theorem (1921) gives a unique probability measure `mu` on `[0,1]`
with `m_K=int r^K dmu`. Since `-Delta r^K = r^K(1-r)`, iterating gives the
display. Uniqueness is moment determinacy on a compact interval. The converse
is `r^K(1-r)^L = r^K(1-r)^(L+1) + r^(K+1)(1-r)^L`. QED.

This is de Finetti's theorem in its Pascal-graph form: the extreme points of
the convex set of split arrays are the fixed prices `r in [0,1]` (Vershik and
Kerov's boundary of the Pascal graph). The Codex arrays `E,W` are the images
of `6r(1-r)dr` and `2(1-r)dr`. The discounted array is
`D(L,K)=2W(L,K)/((L+1)(L+2))`; it instead satisfies
`D(L,K)=D(L,K+1)+(L+3)D(L+1,K)/(L+1)`.
Already `D(0,0)=1` but `D(1,0)+D(0,1)=2/9+1/3=5/9`.

**Corollary A1 (fibre payment for every `mu`; PROVED).** For `L>=1`,
`sum_(j>=0) w(L+1,K+j) = w(L,K)` exactly, with remainder `w(L,K+J)` after `J`
terms. For `L=0` the sum is `w(0,K)-mu({1})`. An atom at `r=1` therefore
leaves a formal `L=0` row defect equal to the atom, and it makes the root ray
non-summable (`w(0,j)>=mu({1})` for all `j`). So for every finite-mass member
the incoming row of every unit target is paid exactly (this is P3 of the
adaptive note for all `mu` at once).

**Corollary A2 (mass; PROVED).** The root ray has mass `int dmu/(1-r)`; the
full flow has mass at most `int (2+2r-r^2)/(1-r) dmu` (P1 integrated). Hence a
finite-mass counter-only flow exists for `mu` iff `int dmu/(1-r) < infinity`.
For `mu=beta(alpha,beta)` this is `beta>1`, with the bound
`3(alpha+beta-1)/(beta-1)-beta/(alpha+beta)` (adaptive note (10); verified
against the integral for six priors in S3).
The root-deleted mixture `Z` of the updated adaptive note (P7) uses the
non-integrable prior `(1-r)/r dr`; it is finite on `V` only because every
nonroot certificate has `K>=1`, and its root-entry ray `1/(j(j+1))` sums to
1. This is outside A2's finite-probability-prior hypotheses: the full nonroot
mass bound requires the separate cost-kernel estimate in P7, not only the
root-entry calculation.

**Corollary A3 (support; PROVED).** If `mu((0,1))>0` then every `w(L,K)>0`,
so the flow `W_mu(n)=w(L(n),K(n))` on rooted sources, zero elsewhere, has
support exactly the rooted component, independently of `mu`. If `mu` is
carried by `{0,1}` the array has zeros. Consequently **among priors charging
the interior, changing the price measure leaves the support unchanged**;
positivity everywhere of any such canonical counter-only flow
is equivalent to universal Collatz. This equivalence supplies no independent
proof of positivity, and does not rule out one by another argument.

**Corollary A4 (the KT prior is not summable; PROVED).** The
Krichevsky-Trofimov prior `beta(1/2,1/2)`, a standard universal binary
mixture, gives root-ray weights `C(2j,j)/4^j` with partial sums
`(2J+1)C(2J,J)/4^J -> infinity`. The general condition is exactly
`int dmu/(1-r)<infinity`; within the beta family it is `beta>1`.
It does not force a positive-power endpoint decay. For example, a density
proportional to `1/log(1/(1-r))^2` on `1-e^(-1)<r<1` is admissible:
substitution `u=log(1/(1-r))` gives the finite integral `int_1^infinity du/u^2`.
This density vanishes more slowly than every positive power of `1-r`.
S3 tabulates the likelihood ratio `w(L,K)/sup_r r^K(1-r)^L` for five
priors. P4 supplies the particular `W` prior's quadratic loss bound;
no uniform quadratic bound for all admissible priors or optimality of
`beta(1/2,2)` is asserted.

Scope. Theorem A says nothing about which `(L,K)` are realized by integers,
nor about the word. It classifies the price layer completely and shows that
changing an interior prior alone leaves its support unchanged. Other
objectives, such as a source-specific lower bound, can still depend on the prior.

## 4. The dictionary, typed

| Source | Target | Map | Preserved | Lost | Sidecar / hostile | Status |
|---|---|---|---|---|---|---|
| Split arrays `w(L,K)` | Harmonic functions on the Pascal graph | identity | split = harmonicity at every vertex | nothing | — | PROVED (definition) |
| Pascal graph | Bratteli diagram of the gauge-invariant CAR algebra; fusion graph of `1+chi` in `Rep U(1)` | level `n=L+K`, charge `K` | multiplicities `C(n,K)` | — | the word tree is the binary (CAR, not gauge-invariant) diagram | CITED (Bratteli 1972; Vershik-Kerov) |
| Traces of that algebra | Probability measures on `[0,1]` | moment map | convexity, extreme points `r` | — | atom at 1 = non-summable root ray (A1) | CITED + Theorem A |
| Codex `E, W` | two interior traces | priors `6r(1-r)`, `2(1-r)` | exact split and fibre payment | word, source guard | 53 vs 113; `D` has a weighted split | PROVED |
| Free monoid `{b,s}*` (words) | `N^2` (counters) | abelianization | exchangeable prices | order, realizability | inherited obstruction concerns specified bounded observations | PROVED (A3 + inherited, scoped) |
| Finite subgroups of `SU(2)` / `U(1)` / `SU(2)` | affine ADE / `A_infinity^infinity` / `A_infinity` McKay graphs | fusion with the defining representation | discrete list vs continuous boundary | dynamics | no Collatz theorem uses index values | CITED analogy only |
| `GL(1)` Langlands (characters, class field theory) | counter-only flows and residue sections | mixtures of characters of `N^2`; discrete log mod `3^(j+1)` | an abelian-coordinate analogy | the word | no constructed reciprocity correspondence | HEURISTIC |
| `S_2(Gamma_0(22))` oldforms of the conductor-11 newform (Codex level-22 note) | the predecessor constant 22 | none found | an exact phase module and graph lift (their section 3) | mass and source (their F3: the atom generating function is not a modular form) | the two 22s have different origins (their section 5) | CITED + their PROVED obstruction |

Two rows deserve a sentence each. The McKay/Jones dichotomy (finite subgroups
give the finite ADE list; the continuous groups give infinite paths with a
continuum of traces) is the shape the owner asked about, and the Pascal
boundary `[0,1]` is its `U(1)` instance, but no Jones-index value or ADE
diagram enters any Collatz statement here. The Langlands row records an
abelian-coordinate analogy, not a reciprocity theorem for Collatz. Exact
finite-dimensional matrices may retain unbounded arithmetic data; the
finite-lookahead obstruction does not exclude them.

## 5. Positivity everywhere: what any certificate must be

**Lemma B (no power-law comparison; PROVED).** Let `v>0` on positive odd
integers satisfy the incoming inequality `Kv<=v` at every nonroot target.
Then for every `l>=1`, `v(2*3^l-1) >= v(2^(l+1)-1)`, because the orbit of
`x_0=2^(l+1)-1` passes through `x_i=2^(l+1-i)3^i-1` with valuation one at
every step (S4 verifies `l<=60`). Hence for no `s>0` and `C` is
`C^(-1) n^(-s) <= v(n) <= C n^(-s)`: the ratio `x_l/x_0` tends to `(3/2)^l`.
The witnesses `x_i` are simultaneously the 2-adic and the 3-adic shadows of
the fixed point `-1` (`x_i = -1 mod 2^(l+1-i)` and `mod 3^i`). The same holds
along the `(-17)` cycle, word `(1,1,1,2,1,1,4)` with `3^7/2^11=1.0679`: the
shadows `2^(11m)t-17` (`t` even) climb to `3^(7m)t-17` (S4, `m<=3`). This is
the honest place where 11 enters the positivity problem: `2^11` against `3^7`,
the only near-coincidence of the torsion census in the eleventh note of the
Collatz thread, and the negative cycle the lookahead note uses for its
aliases.

Lemma B generalizes section 6 of the mixed note (finite-alphabet obstruction
at one normalization) to every power-law normalization: a certificate must be
unboundedly non-uniform against every power of `n`, and the non-uniformity is
forced at the shadows of the negative cycles.

**Proposition C (Riesz decomposition; PROVED).** Let `v>=0` be summable with
`Kv<=v` on the nonroot odd integers `V`, and put `lambda=v-Kv>=0`. Then

    v = G lambda + h,   G lambda = sum_(t>=0) K^t lambda,   K h = h,  h>=0,

where `lambda` is supported on rooted sources, `h` is a nonnegative
combination of indicator functions of nontrivial cycles, and
`sum_V G lambda = sum_z lambda(z) tau(z)` with `tau(z)` the number of odd
steps to 1.

*Proof.* Iterating `v=lambda+Kv` gives `v=sum_(t<T)K^t lambda + K^T v`, so
the series `G lambda` converges to at most `v` and `h:=v-G lambda>=0` is the
decreasing limit of `K^T v`; monotone convergence gives `Kh=h`. The mass of
`K^t lambda` is the `lambda`-mass of sources whose `t`-th image is still in
`V`, so `sum_V G lambda = sum_z lambda(z) #{t: U^t z in V}`; finiteness forces
`lambda=0` wherever that count is infinite, i.e. off the rooted component.
Summing `Kh=h` first forces zero harmonic mass on the predecessors of the
killed root, then recursively on every rooted source. If `n` has a divergent
orbit then `h` is nondecreasing along it and summable
over infinitely many distinct terms, so `h(n)=0`. On a cycle component, sum
`Kh=h` around the cycle: the incoming tree mass vanishes, so `h=0` on the
trees feeding the cycle and `h` is constant on the cycle. QED.

This classifies the harmonic remainder among these nonnegative summable
supersolutions; a Martin boundary requires a separately specified kernel and
is not classified here. It re-proves C3: a nontrivial cycle has an external
predecessor divisible by three, where both its rooted potential and harmonic
remainder vanish, contradicting strict positivity. A divergent source also
has both terms zero. Thus strict positivity implies every source is rooted.
It does not require `lambda(n)>0` at every vertex: the potential is positive
when the vertex lies on the finite future of a source charged by `lambda`.
The finite example `v(3)=v(5)=1`, zero elsewhere, has `Kv<=v`,
`lambda(3)=1` and `lambda(5)=0`, along `3 -> 5 -> 1`.

For the Codex flows `lambda` is supported on the leaves (odd multiples of
three) and the payment at units is exact (A1); P5 of the adaptive note is the
case `mu=2(1-r)dr` with `lambda`-mass `int r/(1-r) dmu = 1`, and the
beta(2,2) case has mass 2 (B2). Positivity everywhere is the statement that
`lambda` charges every leaf.

**Proposition D (the leaf-section tower; PROVED).** For `j>=1` and every
positive odd unit `v`, let `a_j(v) in [1, 2*3^j]` be the least `a` with
`2^a v = 1 mod 3^(j+1)`; it exists and is unique in that range because 2 is a
primitive root modulo every power of three. Then
`rho_j(v) = (2^(a_j(v)) v - 1)/3` is the least predecessor of `v` divisible
by `3^j`, and

    rho_j(v) < 2^(2*3^j) v / 3,

sharp as a ratio along `v = 1 mod 3^(j+1)`. For `j=1` this is the fusion
helpers' Proposition S1 with the table `a_0(v mod 9) = {1:6, 2:5, 4:4, 5:1,
7:2, 8:3}` and the bound `64/3 = 21.33 < 22`. The tower continues
`2^18/3 = 87381.3` and `2^54/3 = 6.0e15`. The least unit predecessor needs only
`a<=4` (two admissible residues `{4,7}` modulo nine), bound `16/3`. S5 checks
all odd units below `2*10^5` for `j<=3`, the sharpness classes, and that the
leaf predecessors of `v` are exactly `S^(3t) rho_1(v)`.

*Reading of the 22 for an induction.* A leaf `z=3m` with image `v=U(z)` and
valuation `a=v_2(3z+1)` is paid by a smaller leaf exactly when `a>a_0(v)`,
since then `rho_1(v)=S^(-3t)(z)<z` for some `t>=1`. The unpaid leaves are the
section `Sigma={z: 3|z, v_2(3z+1)<=6}` itself, which `U` maps bijectively to
the units. So "every unit has a leaf predecessor below `22v`" reduces universal
Collatz to the rootedness of `Sigma`, and `Sigma` is a conjugate copy of the
units (fusion helpers: `rho o U o U`). The bound is a fixed ordinary-size cost
of a change of coordinates, not progress on support. The adversary this
coordinate does not remove is Lemma B's: inside `Sigma` the ascending shadows
are still present.

**Specified finite-observation weights.** The inherited theorem of
[finite lookahead obstruction](collatz_finite_lookahead_weight_obstruction_20261005.md)
excludes its explicitly defined lookahead/cofactor-residue observers, and
Lemma B excludes power-law comparability. A positive summable rooted
potential requires an injection whose finite futures cover every vertex.
It suffices to find a leaf-supported injection positive at every leaf with
`sum lambda(z)tau(z)<infinity`: every unit has a leaf predecessor.
These statements do not require every possible certificate to be
non-exchangeable or infinite-dimensional. The canonical Green construction
is defined using rooted futures. The remaining obligation is an independent
proof that its support is full, or a different positive summable construction.

## 6. Hostile probe: the 11 is not attained, nor the 16/3

Two-sided brackets for the actual masses of `E` (beta(2,2)) and `W`
(beta(1,2)), without any convergence assumption.

*Lower bounds (FINITE-EXACT, float64 sums of exact-counter weights).* For
every odd `n` below `2^22` the actual ROOT word gives the counters. Summing
source weights gives the plain head; adding, for every enumerated base
`b>1`, its entire sibling ray (exact mass `w(L-1,K)` by Corollary A1) and the
full root ray gives the fibre-complete head. The Codex heads below `2^15`
(`W` 2.599557, `lambda_W` 0.572152) are reproduced to nine digits as controls.

*Upper bounds (VERIFIED numerical).* For a fixed price `r`, the base-tree
mass is bounded by an adversarial-lift Bellman value: classes are residues
modulo `3^k`, a base in class `a` has three lifts modulo `3^(k+1)`, and the
lift determines the residue class of every inverse child together with the
forbidden depth class. Let `V(a) = 1 + max_(lift) sum_b ent(lift,b) V(b)`,
where `ent(lift,b)` is the total edge weight `(1-r)r^t` of the children of
that lift landing in class `b`. Any finite fixed point `V>=1` dominates every
depth-truncated adversarial value by induction, hence the actual subtree mass
of every base in class `a`; the first generation from the root is exact by
residue. Policy iteration finds the fixed point (residual and positivity are
checked). On each `r`-interval the edge weights are replaced by their suprema
(each `r^t(1-r)/(1-r^P)` is unimodal; the mode is bracketed), and the smaller
of this bound and the proved `2+2r-r^2` is integrated against the prior. The
naive max-over-lifts matrix is **invalid** at small `r` (it triple-counts the
base child, which each lift steers to a different class); the Bellman form
does not have that defect. Near `r=1` every class passes `2/3` of its mass and
the bound tends to 3, which is also the true limit, so the integrands
`6rB(r)` and `2B(r)` are dominated by the sibling-heavy region, not by any
bound slack.

| quantity | lower bound (bases below `2^22`, full fibres) | upper bound (level `3^6`, 729 classes) | proved bounds in the Codex notes |
|---|---:|---:|---|
| `W`-mass, beta(1,2) | 3.1927 | 3.8907 | `16/3`; `41/10` after the Codex cost-eight kernel (P6) |
| `E`-mass, beta(2,2) | 5.1000 | 6.8835 | `11`, and `17/2` after P1 |
| `lambda_W` (leaves) | 0.6568 (head) | exact 1 | `1` |
| `lambda_E` (leaves) | 1.1113 (head) | exact 2 | `2` |

Pointwise base-mass bounds `B(r)` at level `3^6` against the proved
`2+2r-r^2`: `r=1/16`: 1.079 vs 2.121; `r=1/4`: 1.399 vs 2.438; `r=1/2`:
1.903 vs 2.750; `r=3/4`: 2.465 vs 2.938; `r=0.9`: 2.794 vs 2.990; `r=0.99`:
2.980 vs 3.000. Levels `3^1..3^6` give `E` bounds 7.147, 6.965, 6.926,
6.913, 6.885, 6.884 and `W` bounds 4.140, 3.987, 3.948, 3.934, 3.892, 3.891:
the residue depth has nearly saturated because the dominant region `r -> 1`
is already sharp. The plain source heads below `2^22` are `W` 2.851 and `E`
4.215; the `W`-mass by odd-step layer `tau=1..11` is 1.833, 0.514, 0.213,
0.123, 0.071, 0.043, 0.022, 0.013, 0.008, 0.005, 0.002 (head-truncated).

The beta(2,2) bracket excludes 11 and 17/2; the beta(1,2) bracket excludes
16/3. The 11 is `3!H_3` of a bound, not a mass. Since the brackets are far
from both `22/2` and any other small integer, no identity of the form
"total mass = 11" survives.

## 7. Numerology verdict and the honest elevens

- `11 = 6 + 3 + 2` equals `ord_9(2)+ord_9(4)+ord_9(8)` and `sigma(6)-1`;
  both are restatements of `6(1+1/2+1/3)` with no mechanism beyond
  `ord_9(4)=3` and the choice `1/B(2,2)=6`. NUMEROLOGY.
- `22 = 2*11`: the 2 is `ceil` and the 11 is a prior normalization.
  NUMEROLOGY. Correct statement: `rho_1(v) <= (64v-1)/3`, first rung of
  `2^(2*3^j)/3`.
  The concurrent Codex notes reach the same verdict independently (inverse
  sections, section 1; level-22 bridge, section 5).
- The elevens that are structural in this thread: the eleven halvings of the
  `(-17)` cycle (`2^11` against `3^7`, Lemma B), and `11/7` as the mediant
  convergent of `log_2 3` in the eleventh Collatz note. Neither is the mass
  bound.
- `{2,3,11}` as a set: the Collatz primes are 2 and 3; 11 enters only through
  `3^7-2^11=139` and through the torsion census `(11,7)`. No theorem places a
  third prime in the dynamics.

## 8. Reproduction, scope, limits

[Script](../../04-computation/experiments/collatz_pascal_boundary_leaf_section_20261005.py),
[output](collatz_pascal_boundary_leaf_section_20261005.out),
[JSON](collatz_pascal_boundary_leaf_section_20261005.json):

```text
python3 04-computation/experiments/collatz_pascal_boundary_leaf_section_20261005.py --level-max 6 --head-bits 22 --json 05-knowledge/results/collatz_pascal_boundary_leaf_section_20261005.json
python3 -O 04-computation/experiments/collatz_pascal_boundary_leaf_section_20261005.py --level-max 3 --head-bits 16
```

The full run records 1,811,150 explicit checks in 112 s (numba for the heads
and lift tables; 12-core machine not needed). The universe: counters
`0<=L,K<=12` for the Polya and split identities; six random discrete price
measures on a `10 x 14` rectangle with fibre remainders; the KT partial sums to
`J=299`; six priors against quadrature; the chains `l<=60` and the `(-17)`
shadows `m<=3`; all odd units below `2*10^5` for the tower `j<=3`; all odd
`n<2^22` for the heads; six residue levels with 1,180 `r`-intervals each
(width `10^-3` on `[0,0.98]`, `10^-4` above). Hostiles: the non-monotone row
`(1, 9/10, 1/2)`, the atom at `r=1`, the uniform-prior root ray, the
`(-17)` shadow with an odd multiplier (even multipliers are required for an
odd landing value; the script uses even ones), and the invalid
max-over-lifts matrix, which was observed to exceed spectral radius one at
`r=1/16` and is not used.

Limits. Theorem A classifies formal arrays; realizability of `(L,K)` by
integers and the word are outside it. Lemma B and Proposition C are
unconditional but prove nothing about any particular integer. The Bellman
bounds are float64 with `1e-7` relative margins and interval suprema; they are
not interval arithmetic and not exact rationals. The lower bounds are finite
heads; the polynomial sibling tails converge slowly (the root ray alone loses
`2/(J+1)` beyond `S^J(1)`), which is why the fibre-complete form is used.
No finite number here is treated as an unbounded statement.

## 9. Obligations and next tests

1. Independent audit of Theorem A's use of Hausdorff (the only non-elementary
   input) and of Proposition C's cycle argument.
2. Tighten the brackets by enumerating the base tree with exact fibres to depth
   `G` and bounding only the omitted subtrees with the per-interval Bellman
   values (the machinery is in place; the enumeration is the cost).
3. The only admissible positive target remains the one the Codex notes name:
   a source-dependent lower bound on `lambda(z)` across a refuel boundary for
   a named unbounded family of leaves, retaining the word. Theorem A says no
   price tuning helps; Lemma B says the family must include the ascending
   shadows; Proposition D says it may be taken inside the section `Sigma`.
