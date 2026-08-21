# Rule 30 adaptive observers, rank 30, and p-adic selector gaps

**Session date:** 2026-08-21  
**Status:** bounded-gap recovery, scale-fixed-period no-go, and sextic sidecar
identity `PROVED`; Rule 30 separator census and elliptic local certificate
`FINITE-EXACT`; rank lower bound `PROVED`; literature statements `CITED`; all
six named p-adic zeta singleton claims `OPEN`; proposed search transfers
`OPEN`.

Exact companions:

- [`rule30_signalizer_selected_ray_bank_codex_20260821.py`](../04-computation/rule30_signalizer_selected_ray_bank_codex_20260821.py),
  with [stored output](../05-knowledge/results/rule30_signalizer_selected_ray_bank_codex_20260821.out);
- [`rule30_scale_fixed_period_scheme_no_go_codex_20260821.py`](../04-computation/rule30_scale_fixed_period_scheme_no_go_codex_20260821.py),
  with [stored output](../05-knowledge/results/rule30_scale_fixed_period_scheme_no_go_codex_20260821.out);
- [`elliptic_rank30_optimal_local_observer_certificate_codex_20260821.py`](../04-computation/elliptic_rank30_optimal_local_observer_certificate_codex_20260821.py),
  with [stored output](../05-knowledge/results/elliptic_rank30_optimal_local_observer_certificate_codex_20260821.out);
- [`elliptic_mestre_six_tuple_sidecar_rule30_scout_codex_20260821.py`](../04-computation/elliptic_mestre_six_tuple_sidecar_rule30_scout_codex_20260821.py),
  with [stored output](../05-knowledge/results/elliptic_mestre_six_tuple_sidecar_rule30_scout_codex_20260821.out).

## 1. Outcome first

The common lesson of the three frontiers is that **separation is not closure**.

1. On the 23 physical Rule 30 signalizers `s_0,...,s_22` frozen by
   [THM-3511](../01-canon/theorems/THM-3511-rule30-orbit-signalizer-gap-renormalization-and-shallow-portrait-hostile.md),
   no one or two depth-four ray evaluations separate all signalizers, while
   104 three-ray banks do.  The lexicographically first bank is `(0,1,2)`.
   This finite separator is not a quotient of the first-return operation.

2. The failure has an exact repair.  A depth-`D+B` portrait of a signalizer
   determines its first-return gap and depth-`D` successor portrait whenever
   the gap is at most `B`; otherwise it returns the honest statement
   `gap>B`.  This is a proved **adaptive precision tariff**, not a finite
   signalizer graph or a Rule 30 prize solution.

3. Exact scale-fixed period-five and period-seven projective profiles
   linearize to eigenspaces.  Period five has only projective points; period
   seven has only points and genus-zero lines.  The only physically admissible
   eigenvalue repeats gap one, contradicting the inherited no-consecutive-gap
   law.  Thus neither fixed-profile scheme contains a physical elliptic lane.

4. The displayed elliptic curve in the prompt has 30 public rational points.
   A cardinality-minimal bank of 15 good primes gives a rank-30 map into local
   mod-two quotients.  Together with a separate two-torsion check at `p=23`,
   this proves the 30 points independent and hence proves the unconditional
   lower bound `rank E(Q)>=30`.  It does not prove exact rank 30.

5. A classical six-point elliptic construction has a hidden algebraic
   sidecar.  If a centered monic sextic has elementary symmetric coefficients
   `e_i`, then its shifted product completes to a quartic precisely on the
   nondegenerate codimension-one locus

   ```text
   2 e_5 = e_2 e_3.
   ```

   Centering alone is insufficient.  A Rule 30 symmetric tuple satisfies the
   condition automatically, which is a control rather than evidence of rank
   enrichment.

6. None of the six named p-adic zeta values is presently known irrational as
   an individual number.  The strongest applicable results are disjunctions
   or asymptotic rank statements.  They do not select the advertised
   coordinate.

The useful synthesis is therefore

```text
observer bank + native-operation compatibility + lost-coordinate sidecar.
```

The first item can identify a frozen sample.  Only all three can support a
proof that survives iteration, integer relations, or a singleton quantifier.

## 2. Inheritance pass and live concept board

### Anchor: Rule 30 signalizers

- **Closest proved mechanism:** THM-3511's first-return signalizers, exact
  gaps, and depth-four portrait atlas.
- **Canonical hostile:** `s_7` and `s_17` have the same depth-three portrait
  but gaps `8` and `5`.
- **Corrected near miss:** a shallow portrait can identify the first 23
  physical states without being a congruence for first return.
- **Least-used sidecar:** the gap itself, treated as a precision budget rather
  than as an annotation after the transition.

### Niche: high-rank elliptic curves

- **Closest proved mechanism:** local maps
  `E(Q)/2E(Q) -> E(F_p)/2E(F_p)` at good primes, followed by infinite descent.
- **Canonical hostile:** an injective mod-two image of a rational two-torsion
  point does not prove that point has infinite order.
- **Corrected near miss:** point incidence and local signature rank prove only
  what their torsion and reduction sidecars permit.
- **Least-used sidecar:** choose local primes and future point searches to add
  complementary signature dimensions, rather than using local point counts
  only as a scalar score.

### Wildcard: p-adic zeta singleton irrationality

- **Closest proved mechanism:** finite disjunctions at small weights and
  asymptotic many-value irrationality ranks.
- **Canonical hostile:** `OR_i Irr(zeta_p(i))` has no active-coordinate label.
- **Corrected near miss:** an existential irrationality theorem is not a
  coordinatewise theorem.
- **Least-used sidecar:** a labelled nonzero coefficient or minor that cannot
  vanish after projecting away the other weights.

### Board

The session compared five objects after each pull:

```text
finite separators;
operation congruences;
adaptive precision/overflow;
local independence matrices;
labelled existential selectors.
```

The new results move the board from “find a smaller statistic” to “measure the
minimum extra state needed for the statistic to commute with the operation.”

## 3. A minimum selected-ray bank for the first 23 Rule 30 signalizers

Let `pi_D(s)` be the permutation of the `2^D` binary rays induced by a
signalizer `s`.  On the exact physical universe `s_0,...,s_22` from THM-3511,
the companion enumerates all coordinate subsets of sizes one, two, and three
in the depth-four portrait.

The exact census is

```text
bank size 1:   0 separators;
bank size 2:   0 separators;
bank size 3: 104 separators.
```

The first bank is `(0,1,2)`.  It assigns 23 distinct triples to the 23 frozen
signalizers.  This proves minimum size three **only in this explicit finite
universe**.

### The bank does not carry first return

Two ambient active words already collide:

```text
word    selected output   gap   selected successor output
CA      (7,10,1)           3    (13,2,3)
CCAC    (7,10,1)          10    (9,4,7).
```

This is not repaired merely by retaining every ray at a fixed shallow depth:

```text
depth 4: CBC and AAAACBC have one portrait, but gaps 6 and 8;
depth 5: AACAAACC and CAAACAAC have one portrait, but gaps 6 and 5.
```

Thus the first failed implication is

```text
injective on a frozen physical sample
    != congruence for first-return renormalization.
```

## 4. Proved bounded-gap adaptive recovery

There is nevertheless a simple exact theorem behind the failures.

Let `d(s)` be the first moved bit on the zero ray under `s^2`, and let `R(s)`
be the section of `s^2` below the fixed prefix `0^d`.  For integers `D,B>=1`,
the depth-`D+B` portrait `pi_(D+B)(s)` determines exactly one of:

```text
d(s)<=B, together with d(s) and pi_D(R(s));
d(s)>B, reported as overflow.
```

Indeed, the source portrait determines the square portrait.  Inspecting the
image of zero through bit `B` either locates its first nonzero bit `d`, or
proves there is none in the visible prefix.  In the first case, prefix
compatibility fixes `0^d`, and for `0<=x<2^D`,

```text
R(s)(x) = (s^2(2^d x) >> d) mod 2^D.
```

Every input on the right lies inside the known depth-`D+B` tree when `d<=B`.
This proves the recovery formula without extrapolation.

The three hostile pairs give exact controls:

```text
target D   pair                    split B   full B   source depth
4          CA / CCAC                 3        10          14
4          CBC / AAAACBC             6         8          12
5          AACAAACC / CAAACAAC       5         6          11.
```

At the split budget, one state is exact and the other overflows.  At the full
budget, both exact successors are recovered.

The conceptual boundary is sharp.  A fixed-depth observer can be uniformly
closed under first return only after supplying a uniform gap bound or a
self-similar representation of arbitrary overflow.  The former is entangled
with the open bounded-gap/quasisymmetry frontier; the latter is the more
promising carry/section target.  This theorem does not establish either.

### 4.1 Scale-fixed periods five and seven are not elliptic

The first coefficient-variety experiment closes negatively, but for a useful
reason.  Write

```text
g_j=G_m(t+jq_m),        q_(m+1)=2q_m.
```

THM-3512's three-local recurrence sends a period-`n` profile to

```text
g'_j=-g_(2j)g_(2j+1)(1-g_(2j+2))/(1-g_(2j)),
```

with indices modulo odd `n`.  On the saturated locus `g_j(1-g_j)!=0`,
multiplying the fixed equations gives `product_j g_j=-1`.  Hence cyclic
amplitudes exist with

```text
g_j=-A_(j+1)/A_j.
```

Put `(T_nA)_j=A_(2j)+A_(2j+1)`.  Substitution turns `g'=g` into

```text
(T_nA)_(j+1) A_j = (T_nA)_j A_(j+1),
```

so the fixed locus is an open part of the projectivized eigenspace locus of
`T_n`.  The exact characteristic polynomials are

```text
n=5: (lambda-2)(lambda-1)(lambda+1)(lambda^2+1),
n=7: (lambda-2)(lambda-1)^2(lambda^2+lambda+1)^2.
```

At period five every eigenspace is one-dimensional, so the saturated locus
is zero-dimensional.  At period seven the repeated eigenspaces projectivize
to lines; every positive-dimensional component is genus zero.  Neither
period contains an elliptic component.

The physical filter is stronger.  All amplitude ratios are odd units, so
their valuations agree.  If `T_nA=lambda A`, then

```text
nu_2(lambda)=nu_2(1-g_j)=d>=1.
```

Every eigenvalue above except `2` is a 2-adic unit.  The `lambda=2`
eigenspace is the constant amplitude line, giving `g_j=-1` and `d=1`.
Exact scale-fixedness repeats this gap at consecutive scales, contradicting
THM-3512's proved no-consecutive-gap-one law.  Therefore there is no physical
exact scale-fixed period-five or period-seven profile.

This rules out the cheapest attempt to import elliptic geometry.  A lawful
next attempt must relax exact fixedness to a multi-scale cycle or introduce a
marked/twisted holonomy before asking whether a genus-one component appears.

## 5. The rank-30 curve: a minimal local observer certificate

For the curve

```text
y^2+xy=x^3+A x+B
```

with the coefficients in the prompt, the public ICARM record supplies 30
rational points and reports rank at least 30.  The companion independently
checks every incidence and builds a certificate using the 15 good primes

```text
43, 61, 101, 211, 223, 241, 263, 271,
283, 311, 313, 457, 521, 569, 571.
```

At every selected prime, `E(F_p)/2E(F_p)` has dimension two.  The stacked
`30 x 30` binary response matrix has rank 30.  Fifteen primes are therefore
cardinality-minimal: an elliptic mod-two quotient contributes at most two
dimensions.  Increasing-prime greedy selection produces the lexicographically
least such bank among odd good primes `<=571`, the exact audited universe.

### Why the torsion sidecar is load-bearing

Let `P_1,...,P_30` be the public points.  If

```text
sum_i n_i P_i = O,
```

injectivity of the stacked local map makes every `n_i` even.  Write
`n_i=2m_i`; then `Q=sum_i m_iP_i` is rational two-torsion.  At the separate
good prime `23`, the exact group order is `33`, so good reduction gives
`E(Q)[2]=0`.  Hence `Q=O`, and infinite descent makes the integer coefficient
vector divisible by every power of two.  It is zero.

Therefore the 30 points are independent and

```text
rank E(Q) >= 30
```

unconditionally.  The argument proves no rank upper bound.  The public claim
of exact rank 30 is conditional on additional analytic hypotheses; it must
not be blended with this lower-bound certificate.

This is the precise arithmetic counterpart of the Rule 30 lesson.  The local
bank separates coefficient parity classes, while the torsion observer makes
the separation compatible with division and infinite descent.

## 6. A codimension-one filter for six-point elliptic searches

The current ICARM curve 211 commentary records a Mestre--Fermigier search from
the integer sextuple

```text
(348,-600,-216,492,876,-900)
```

and shift `t=2326/23`.  The algebra behind that construction exposes a useful
filter.

Let

```text
p(z)=z^6+c_4z^4+c_3z^3+c_2z^2+c_1z+c_0
```

be centered, put `q(x)=p(x-t)p(x+t)`, and choose the unique monic degree-six
polynomial `g` for which `deg(q-g^2)<=5`.  Triangular coefficient matching
gives

```text
[x^5](g^2-q) = -12 t^2 (2c_1-c_3c_4).
```

In elementary symmetric coordinates, centering says `e_1=0` and
`(c_4,c_3,c_1)=(e_2,-e_3,-e_5)`.  On `t!=0`, the remainder drops to degree at
most four exactly when

```text
2e_5=e_2e_3.
```

The twelve abscissae `a_i+-t` are then forced rational incidences on the
quartic.  Distinctness, nonzero discriminant, and independence remain
separate sidecars.

Three exact controls clarify the scope.

- The ICARM tuple above satisfies the equation and gives a squarefree quartic
  with twelve distinct base abscissae.
- The centered tuple `(25,71,135,149,221,-601)` fails the sidecar by
  `15889582604000`; at `t=11` the remainder is a squarefree quintic.  Thus
  zero sum alone does not even stay in the elliptic family.
- [THM-3458](../01-canon/theorems/THM-3458-rule30-right-edge-2-adic-odometer-and-moving-observer-boundary.md)'s
  first right-edge integers `(1,7,25)` give the symmetric tuple
  `+-(1,7,25)`.  At `t=2` it gives a squarefree quartic and twelve base
  abscissae.  But every centrally symmetric sextuple makes `c_3=c_1=0`, so
  this success is automatic and is not Rule 30 enrichment.

The best generation program is therefore not “turn Rule 30 numbers into a
curve.”  It is:

1. search the non-symmetric integral locus `e_1=0, 2e_5=e_2e_3`;
2. reject collisions and singular quartics exactly;
3. score shifts with a Mestre--Nagao prime sum;
4. search for extra rational points;
5. prioritize new points whose local mod-two signatures add dimensions to a
   fixed certificate bank;
6. certify independence with torsion and height/local sidecars.

Rule 30 may schedule the finite CRT/address boxes in step 1 or 3.  A bijective
reordering alone cannot improve an exhaustive search, so it must be benchmarked
against random and Gray orders at equal budget.  The measured signal should be
earlier discovery of new local-signature dimensions or better Nagao scores,
not resemblance to a Rule 30 trace.

## 7. The six p-adic zeta claims: exact quantifiers

With the standard diagonal Kubota--Leopoldt convention

```text
zeta_p(s)=L_p(s,omega^(1-s)),
```

the bounded primary-source audit found the following status as of 2026-08-21.

| Value | Status | Strongest applicable result |
|---|---|---|
| `zeta_5(3)` | `OPEN` | some odd weight in `{3,5,7,9,11,13}` is irrational |
| `zeta_7(3)` | `OPEN` | the same six-way disjunction |
| `zeta_2(7)` | `OPEN` | one of weights `7,9,11,13` is irrational |
| `zeta_2(9)` | `OPEN` | the same four-way disjunction |
| `zeta_3(5)` | `OPEN` | no singleton theorem located; `zeta_3(3)` is known |
| `zeta_3(7)` | `OPEN` | asymptotic many-value theorems select no fixed weight |

The exact literature ledger and broader context are in
[`zagier-weight5-padic-zeta-quantifier-and-measure-frontiers-codex-20260821.md`](zagier-weight5-padic-zeta-quantifier-and-measure-frontiers-codex-20260821.md).
The recent theorem proving `zeta_2(5)` irrational is an important positive
control, but it does not change either queried `2`-adic value.

Write a published disjunction as the Boolean statement

```text
OR_(i in S) b_i = 1,      b_i = 1_{Irr(zeta_p(i))}.
```

It preserves nonemptiness and destroys the owner `i`.
[THM-3516](../01-canon/theorems/THM-3516-rule30-marked-van-der-put-carry-and-power-section-bridge.md)
has an exact Rule 30 hostile of the same type: `J` and `-J` preserve valuations, prefix
fibres, and projective ratios, yet change the marked center bit at `n=3`.
The missing zeta object is therefore not another unlabelled rank lower bound;
it is a labelled coefficient or nonzero minor that isolates one weight and
survives all cancellations.  The exact Rule 30 analogue is THM-3516's marked
signed coordinate.

## 8. Typed connection ledger

| Source | Observer map | Preserved | Destroyed | Required sidecar |
|---|---|---|---|---|
| Rule 30 signalizer | selected portrait rays | identity on 23 frozen states | gap and successor | adaptive gap/overflow and section |
| elliptic points | reductions to `E(F_p)/2E(F_p)` | parity independence | integral heights and torsion | good-reduction torsion gate |
| sextuple search | `(e_1,2e_5-e_2e_3)` | elliptic quartic locus | smoothness and point independence | discriminant, heights, local matrix |
| p-adic zeta family | existential irrationality rank | some coordinate survives | active weight | labelled coefficient/minor |

The connection is methodological, not an identification of objects.  No map
from Rule 30 states to Mordell--Weil points or p-adic zeta values has been
proved.

## 9. Ranked next experiments

1. **Adaptive Rule 30 observer compiler.**  For physical `s_m`, measure the
   least source depth needed to recover the next selected-ray bank, not just
   the least depth that distinguishes `s_m`.  Seek a recursive overflow state
   coming from
   [THM-3512](../01-canon/theorems/THM-3512-rule30-van-der-put-haar-cocycle-and-profinite-automaton-boundary.md)
   and THM-3516's carry data.

2. **Multi-scale period schemes rather than cosmetic tournaments.**  The
   exact scale-fixed period-five/seven schemes are now closed by the linear
   no-go above.  Form two-scale cycles and marked/twisted holonomy schemes.
   A genus-one component with a rational point would create a lawful elliptic
   object; it would still need amplitude and Mealy-section realization before
   saying anything about physical Rule 30.

3. **Non-symmetric sextuple census.**  Enumerate bounded primitive solutions
   of `e_1=0, 2e_5=e_2e_3`, quotient affine symmetries, and compare
   Rule-30/Gray/random address schedules under the same Nagao and point-search
   budget.

4. **Certificate-guided curve search.**  In an established elliptic surface,
   use CRT specialization boxes to maximize the incremental rank of a small
   local signature bank for known and newly found sections.  Treat this only
   as triage until exact incidence and torsion/height certification.

5. **Labelled p-adic zeta selector.**  Revisit the published linear forms
   before scalar rank extraction.  Search for a character-resolved minor whose
   nonvanishing isolates one requested weight.  A mere larger disjunction is
   not progress on a singleton.

## 10. Failure boundaries

- The three-ray bank is minimal only for 23 named depth-four portraits.
- The adaptive recovery theorem charges the actual gap; it gives no uniform
  bound and no finite graph.
- The period-five/seven no-go concerns exact scale-fixed profiles; multi-scale
  cycles and twisted holonomy remain open.
- The elliptic certificate proves `rank>=30`, not exact rank 30 and not a
  construction mechanism for the displayed curve.
- Twelve forced quartic incidences need not be independent.
- A symmetric Rule 30-derived sextuple is a positive algebra control, not
  evidence of search bias.
- Every one of the six named p-adic zeta singleton irrationalities remains
  open under the audited convention.
- The three cross-frontier maps are search and proof-design transfers, not
  arithmetic identities.

## 11. Primary external sources

- [Wolfram Rule 30 prize problems](https://rule30prize.org/): the three named
  center-column problems remain the actual targets.
- [ICARM curve 273](https://elliptic-rank.icarm.cloud/curve/273), its
  [machine-readable record](https://elliptic-rank.icarm.cloud/curve/273.json),
  [verification contract](https://elliptic-rank.icarm.cloud/api), and
  [current verifier source](https://github.com/icarm/elliptic-rank/blob/main/src/verify.ts).
- [ICARM curve 211](https://elliptic-rank.icarm.cloud/curve/211), for the
  exact Mestre--Fermigier tuple and shift used as the positive control.
- [Dujella's rank-30 record page](https://web.math.pmf.unizg.hr/~duje/tors/rk30.html)
  and [rank history](https://web.math.pmf.unizg.hr/~duje/tors/rankhist.html),
  for the lower-bound/exact-rank distinction.
- [Calegari](https://arxiv.org/abs/math/0408214),
  [Lai](https://arxiv.org/abs/2304.00816),
  [Lai--Sprang](https://arxiv.org/abs/2306.10393),
  [Lai--Lupu--Sprang](https://arxiv.org/abs/2505.23088), and
  [Lai--Sprang--Zudilin](https://arxiv.org/abs/2505.05005), for the exact
  p-adic zeta quantifiers used above.

## 12. Reproduction

From the repository root, run

```bash
python3 04-computation/rule30_signalizer_selected_ray_bank_codex_20260821.py
python3 -O 04-computation/rule30_signalizer_selected_ray_bank_codex_20260821.py
python3 04-computation/rule30_scale_fixed_period_scheme_no_go_codex_20260821.py
python3 -O 04-computation/rule30_scale_fixed_period_scheme_no_go_codex_20260821.py
python3 04-computation/elliptic_rank30_optimal_local_observer_certificate_codex_20260821.py
python3 -O 04-computation/elliptic_rank30_optimal_local_observer_certificate_codex_20260821.py
python3 04-computation/elliptic_mestre_six_tuple_sidecar_rule30_scout_codex_20260821.py
python3 -O 04-computation/elliptic_mestre_six_tuple_sidecar_rule30_scout_codex_20260821.py
```

Each ordinary and optimized stream is byte-identical to its checked output.
