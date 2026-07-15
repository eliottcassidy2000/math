---
id: THM-856
title: Hunter tree functional at the seven-comb wall — the first-moment schema is empty, while the ideal independent-density Hunter coefficient is positive for seven combs and negative for eight; exact pair overlap is a two-sided mod-13 sawtooth around 4/169, so uniform radius-seven closure still requires projective-ratio and restricted-overlap control
status: PROVED Hunter/Kounias functional, first-moment no-go, ideal-density coefficient, exact one-comb periodicity, corrected exact global pair-overlap formula, the `2c_E/g` projective pullback bound, and the exact node/edge anomaly decomposition + FINITE-EXACT radius-seven pilots. The original S312 claims that global pair overlap is at least 4/169, that its leading term is `(a+b)^2/(169ab)`, that near-equal speeds are precisely the global minimum, and that a raw-speed `C(E)/x_min` deficit lemma can hold uniformly are REFUTED by the S14 correction below. NOT a closure of the radius-seven chart.
source: opus-2026-07-15-S312; corrected by codex-2026-07-15-S14 after live-pull referee; node/edge defect decomposition added by codex-2026-07-15-S15
depends_on:
  - THM-815 Part C   # the recursion whose union bound dies at 7 combs
related: [THM-778 (mechanical words — the near-equal residual's tool), THM-855 F6 (the moment-closure lens that led here), LRC14-FRONTIER item 3]
verification: 05-knowledge/results/seven_comb_resonance_pilot_opus_S312.out, hunter_tree_wall_crossing_opus_S312.out, hunter_pair_overlap_exact_referee_codex_S14.out
---

# THM-856 — the Hunter tree bound at the seven-comb wall

## 1. The no-go (why THM-815's schema is empty at m′ ≥ 7)

Any bound of the schema |I ∩ D_x| ≤ αL + β(x) (single comb vs single interval,
β(x) → 0) must have α ≥ 2/13: the comb D_x has density exactly 2/13, and
|I ∩ D_x|/L → 2/13 as x → ∞ (equidistribution), so α < 2/13 is FALSE for
large x. At α = 2/13 and m′ = 7 remaining combs, the coverage constraint reads
14L/13 + Σβ ≥ L — satisfied identically: NO constraint on the lifts. The same
holds restricted to any fixed safe set E in place of I (the restricted density
also tends to 2/13 — see the periodicity lemma below). First-moment schemas
cannot cross the wall; this is a theorem, not an obstacle. ∎

## 2. The crossing (Hunter's inequality)

Hunter–Kounias: for any events A_i and ANY spanning tree T on the index set,
μ(∪A_i) ≤ Σμ(A_i) − Σ_{(i,j)∈T} μ(A_i ∩ A_j). Applied to A_i = D_{x_i} ∩ E:

> **uncovered ≥ μ(E) − Σ_i μ(D_i ∩ E) + max_T Σ_T μ(D_i ∩ D_j ∩ E).**

If all single masses have their independent-density value `(2/13)mu(E)` and
the selected tree overlaps have `(4/169)mu(E)`, the ideal Hunter coefficient
at `m'` combs is

> **1-2m'/13+4(m'-1)/169=(165-22m')/169.**

For integer `m'` this is positive through `m'=7` and negative from `m'=8`.
At seven, the tree sum is `24mu(E)/169` versus the `13mu(E)/169` needed to
repair the union bound, leaving `11mu(E)/169`.  At eight, `28<39` and the
ideal tree functional has its own wall.  This proves the coefficient
calculation, not that every actual radius-seven packet has sufficiently large
restricted pair overlaps. ∎

## 3. S14 correction: the exact global pair overlap is a trapezoid

Put `D_u={t:||ut||<=1/13}`.  Write `x_1=ga`, `x_2=gb` with `(a,b)=1` and,
without loss, `a<=b`.  Multiplication by `g` preserves Haar measure, so the
answer depends on the reduced projective pair `(a:b)`, not on the raw scale.
Define

```text
T(z)=sum_(m in Z) (z-|m|)_+,
psi(theta)=theta(1-theta),       0<=theta<1.
```

Then the exact formula is

```text
mu(D_x1 intersect D_x2)
 = [T((a+b)/13)-T((b-a)/13)]/(ab)                      (3.1)
 = 4/169
   +[psi({(a+b)/13})-psi({(b-a)/13})]/(ab).            (3.2)
```

Indeed a tooth centred at `j/a` and one centred at `k/b` have centre offset
`m/(ab)`, where `m=bj-ak`.  Coprimality makes the tooth pairs bijective with
`m mod ab`.  Since `2(a+b)<13ab`, every contributing residue has a unique
least integer representative and no contributing overlap wraps around the
residue circle.  The overlap is a **trapezoid**, not the triangle used in the
original S312 draft: after multiplication by `ab` its length is

```text
((a+b)/13-|m|)_+ - ((b-a)/13-|m|)_+.
```

Summing gives (3.1).  If `r=floor(z)` then
`T(z)=(2r+1)z-r(r+1)=z^2+psi({z})`; subtracting the two squares gives
`4ab/169` and proves (3.2).  Equivalently, with
`Q(c)=r(13-r)` for `r=c mod 13` in `{0,...,12}`,

```text
mu(D_x1 intersect D_x2)
 = 4/169 + [Q(a+b)-Q(b-a)]/(169ab),                    (3.3)
|mu-4/169| <= 42/(169ab).
```

The error has both signs.  Exact counterexamples to the original lower bound
are

```text
mu(D_6 intersect D_7)=2/91 < 4/169,
mu(D_5 intersect D_9)=4/195 < 4/169.
```

Thus `(a+b)^2/(169ab)` is not the leading term; it results from omitting the
containment plateau when the narrower tooth lies inside the wider one.
Near-equality is not a global minimizer theorem.  The sign is the mod-13
sawtooth `Q(a+b)-Q(b-a)`, and the restricted overlap inside a prefix safe set
has an additional endpoint-correlation term.

There is an exact scale-normal estimate for that term.  If `E` is a union of
`c_E` intervals and `B_(a,b)=D_a intersect D_b`, then

```text
|mu(E intersect D_(ga) intersect D_(gb))
   -mu(E)mu(B_(a,b))| <= 2c_E/g.                        (3.4)
```

On every full cell `[j/g,(j+1)/g)`, multiplication by `g` maps the pullback of
`B_(a,b)` linearly onto `B_(a,b)`, so its mass is exactly
`mu(B_(a,b))/g`.  Each component of `E` has at most two partial boundary
cells; bounding their actual and expected contributions by the cell length
proves (3.4).  Thus common scale `g`, reduced ratio `(a:b)`, and the mod-13
sawtooth are the natural edge coordinates.  Raw speeds alone are not.

## 3.5 The exact node-coloured Hunter defect algebra

The correction above admits a useful exact reorganization.  Let `e=mu(E)`
and let `x_1,...,x_m` be the remaining comb frequencies.  Give vertex `i` the
single-comb anomaly

```text
s_i = mu(E intersect D_(x_i)) - 2e/13.                  (3.5)
```

For every edge `ij`, write `x_i=g_ij a_ij`, `x_j=g_ij b_ij` with the reduced
pair coprime, and define

```text
h_ij = [Q(a_ij+b_ij)-Q(b_ij-a_ij)]/(169 a_ij b_ij),
eta_ij = mu(E intersect D_(x_i) intersect D_(x_j))
         - e mu(D_(a_ij) intersect D_(b_ij)),
c_ij = e h_ij + eta_ij.                                 (3.6)
```

Here `h_ij` is the projective mod-13 defect from independent pair density,
while `eta_ij` is the endpoint/pullback defect and satisfies
`|eta_ij|<=2c_E/g_ij`.  Substitution into Hunter--Kounias gives the exact
decomposition of its lower bound on the uncovered part of `E`:

```text
L_H(E;x_1,...,x_m)
 = e - sum_i mu(E intersect D_(x_i))
     + max_T sum_(ij in T) mu(E intersect D_(x_i) intersect D_(x_j))
 = (165-22m)e/169 - sum_i s_i + MST(c),                 (3.7)
```

where `MST(c)=max_T sum_(ij in T)c_ij`.  In particular, at the seven-comb
wall,

```text
L_H = 11e/169 - sum_i s_i + MST(c).                     (3.8)
```

Thus the residual is not a statistic of seven raw speeds.  It is a coloured
complete graph with vertex colours `s_i` and edge colours
`(a_ij:b_ij,g_ij,h_ij,eta_ij)`, evaluated by a tropical spanning-tree
character.  This is the precise sense in which pair data must remain joined
to its endpoint incidences.

There is also an exact recursive classification of the edge information that
the Hunter evaluator uses.  Let `lambda_1>...>lambda_q` be the distinct values
of `c_ij`, let `F_l={ij:c_ij>=lambda_l}`, and put
`r_l=m-kappa(F_l)`, where `kappa` is the number of connected components of
the threshold graph and `r_0=0`.  Kruskal's theorem gives

```text
MST(c)=sum_(l=1)^q lambda_l (r_l-r_(l-1)).               (3.9)
```

Only connectivity-increasing edge levels contribute, but which levels do so
depends on their full incidence pattern.  For example, a connected graph of
nonnegative `c_ij` edges implies `MST(c)>=0`; it does not follow from the
number or average of nonnegative edges.  Equations (3.7)--(3.9) turn the open
uniform step into a finite-type projective edge-classification problem once
the rational-prefix periodic tables for `s_i` and `eta_ij` are fixed.

## 4. The periodicity lemma (the finite table for E-restricted masses)

For a prefix safe set E with rational endpoints and a scale-one lift
x = r + 13h: **x·(|E ∩ D_x| − (2/13)μ(E)) is EXACTLY periodic in h** (period
dividing an explicit Λ(E, r); verified: prefix {1,...,5}, r = 6: period 60,
exact). The per-comb data of the infinite radius-7 chart is a finite table of
rationals. *Proof:* the anomaly depends only on the positions of E's endpoints
in the comb's tooth-coordinate, i.e. on x·(endpoint) mod 1, which is periodic
in h with period = denominator/gcd data. ∎

## 5. Pilot verification (prefix {1,2,3,4,5}, exact ℚ)

μ(E) = 7/15, 10 components. Four random radius-7 packets (lifts h ∈ [2,40]):
Hunter bound = +0.033, +0.034, +0.036, +0.038 — ALL COERCIVE (non-coverage
PROVED per packet by one inequality; actual uncovered ≈ 0.13–0.14, i.e. the
independence prediction (11/13)⁷·μE = 0.145 is accurate to 9%). Consecutive
packets {499..505}, {32..38}: Hunter −0.001, −0.017 while actual uncovered =
0.055, 0.040 — still far from tight.  These are exact finite witnesses that
the functional can succeed and fail.  In view of (3.3), they do not prove
that consecutive packets are the unique failure locus.

## 6. What remains for a radius-7 closure (named residuals)

1. **Projective-ratio stratification.**  A bound by raw speed alone is false.
   On the legal scale-one ray

   ```text
   g=1+13k,       x=6g=6+13(6k),       y=7g=7+13(7k),
   ```

   the reduced pair is always `(6,7)` and
   `mu(D_x intersect D_y)=2/91<4/169`.  For every fixed finite union `E` of
   intervals,

   ```text
   mu(E intersect D_(6g) intersect D_(7g))
      -> mu(E) 2/91,
   ```

   because `D_(6g) intersect D_(7g)` is the pullback of
   `D_6 intersect D_7` by multiplication by `g`.  Hence the originally
   proposed `C(E)/min(x,y)` deficit from the ideal `4mu(E)/169` has a
   nonzero limiting left side and a zero right side.  The next tree lemma
   must retain every reduced ratio `(a:b)` (or prove that the maximum tree
   can avoid its deficient edges), with limiting edge weight given by (3.3).
2. **Boundary discrepancy at fixed reduced ratio.**  For `x=ga,y=gb`, use the
   proved split of the restricted edge weight into

   ```text
   mu(E) mu(D_a intersect D_b)
   + [endpoint discrepancy of the g-fold pullback].
   ```

   Equation (3.4) bounds the second term by `2c_E/g`; for rational endpoints
   its scaled value is periodic.  What remains is a maximum-tree lemma that
   combines these edgewise errors without losing the positive `11/169`
   ideal margin on the admissible projective-ratio graph.
3. **Correlated/AP-window residual.**  The two consecutive pilot packets evade
   the tree certificate but remain far from covering.  Once the projective
   densities and boundary errors are separated, use an AP-window/mechanical-
   word argument (THM-778) on any genuinely consecutive residual rather than
   attributing it to a nonexistent global Bezout minimum.

**The surviving replacement potential is Hunter's weighted-tree functional.
Its ideal-density coefficient crosses the first-moment wall at seven combs,
and the exact pilots show that it can certify real packets.  Its uniform input
is not a scalar raw-speed deficit: it is the reduced-ratio, mod-13 sawtooth,
and endpoint-discrepancy packet above.**

## 7. Exploratory cluster pilots from S312 (not an asymptotic theorem)

The companion cluster script adds three useful finite observations on
`E_[5]`, but its original “complete asymptotic schema” language is not proved.
In particular, the script labels a calibration as such and estimates one
integral by midpoint sampling; neither step is an exact all-packet argument.

1. At `x=500`, differences `d=1,2,3,5` give restricted overlaps below the
   independent baseline, while the sampled `d=8,13,40,200` rows are closer to
   it.  This supports a beat-window diagnostic.  It does **not** prove a
   threshold at `d=8`: equations (3.3)--(3.4) show that raw difference alone
   cannot determine the edge weight.
2. Three displayed multi-cluster packets have positive Hunter certificates
   when a chosen tree uses inter-cluster edges.  This is a three-packet exact
   pilot, not a proof that every large raw separation gives baseline overlap;
   projective rays with a large common scale are the explicit obstruction.
3. The seven-consecutive-speed avoid function has positive sampled integral
   over `E_[5]` (about `0.0558`) and the recorded spot values expose narrow
   denominator-seven windows.  Turning this into a uniform floor still needs
   exact integration over its Farey cells and a quantitative skew-
   equidistribution error.

Thus the cluster experiment names plausible mechanical-word sublemmas but
does not close the asymptotic radius-seven chart or reduce the remaining work
to finite bookkeeping.  Its stored numbers remain reproducible finite data.
