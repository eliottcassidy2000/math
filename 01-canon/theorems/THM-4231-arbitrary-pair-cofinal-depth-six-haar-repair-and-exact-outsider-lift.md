---
id: THM-4231
title: "Literal-patched all-ray cofinite pair reduction and sharp Haar sidecars"
status: >
  PROVED RELATIVE TO THM-4150/4156/4170/4191/4207/4211 + VERIFIED-EXACT +
  TARGETED INDEPENDENT BOUNDARY AUDIT. For the displayed thirty-label pool,
  every nine-body becomes 1/14-Haar-safe after adjoining any two distinct
  outsiders q,r with max(q,r)>=770. A complete clean-room census of
  18,255,923,400 fixed-ray/body cases first leaves exactly 181,242 pairs
  unresolved by both orientations, with certificate cutoff 825. A cumulative
  literal census then exhausts all 45 certificate edges whose larger endpoint
  is at least 770: all 643,821,750 pair/body cases are strictly safe. Together
  with three inherited full-pool literal closures, the exact remaining proof
  residual has 181,194 edges, all in max(q,r)<=769. The complete per-ray table
  has targeted, not full-table, confirmation by the independent split-zeta
  implementation.
  The older q,r>=1290 hybrid theorem and its direct cutoff 1307 remain
  independently audited sidecars; none of 770, 825, 1290, or 1307 is claimed
  as a minimal literal threshold.
  The exact fixed-one ray is independently closed for every r>=542 by
  combining its certificate census with a literal check of the unique
  threshold-543 body at r=542. The
  depth-six activation transition 17547/17548 remains a valid but superseded
  repair-hypergraph sidecar. The 181,194-pair finite remainder and LRC(14)
  remain OPEN.
source: codex-frontier-synthesis-20260826
depends_on:
  - THM-4150-safe-set-haar-measure-universal-odd-tail-lrc14-transfer
  - THM-4156-divisor-complete-anchor-pool-haar-odd-tail-transfer
  - THM-4170-triple-deletion-matching-eventual-haar-odd-tail-transfer
  - THM-4191-complete-full-pool-newcomer-haar-transfer
  - THM-4207-two-newcomer-sharp-depth-transition-base-surplus-composition-and-variable-pool-chart-number
  - THM-4211-fixed-fifty-cofinal-two-newcomer-haar-tail-and-eighteen-label-chart
related:
  - THM-4214-two-newcomer-pascal-complete-eleven-body-haar-charts
  - THM-4223-cyclic-cut-cover-boolean-mobius-hierarchy-and-two-bad-owner-obstruction
  - THM-4227-two-outsider-scale-separated-depth-eight-haar-wedge
  - THM-4228-common-gcd-two-outsider-periodic-observable-haar-ray
  - THM-4233-pair-specific-primitive-observable-oscillation-haar-charts
all_ray_primary_script: 04-computation/lrc14_all_fixed_outsider_ray_census_primary_thm4231.cpp
all_ray_primary_output: 05-knowledge/results/lrc14_all_fixed_outsider_ray_census_primary_thm4231.out
all_ray_primary_script_sha256: 5e4dba0dc8724514aa4d2864ce5fa2b2449c28e39bda5d67672a1723820a39f8
all_ray_primary_output_sha256: a3dcb030aa7891921fa52c1bffc8497b095058e4a666802e1ccd2b14dee1e5d4
all_ray_independent_audit_script: 04-computation/lrc14_all_fixed_outsider_ray_census_independent_audit_thm4231.cpp
all_ray_independent_audit_output: 05-knowledge/results/lrc14_all_fixed_outsider_ray_census_independent_audit_thm4231.out
all_ray_independent_audit_script_sha256: a1ce292a085eee393cbc067f1a493ed8f63324b2fb3d34d460c34c38963dec74
all_ray_independent_audit_output_sha256: dd276f7927a7332907c9bc6215c967817dd463c63fc032c8094a28f16b36dc41
all_ray_symmetric_postprocess_script: 04-computation/lrc14_all_fixed_outsider_ray_symmetric_postprocess_thm4231.py
all_ray_symmetric_postprocess_output: 05-knowledge/results/lrc14_all_fixed_outsider_ray_symmetric_postprocess_thm4231.out
all_ray_symmetric_postprocess_script_sha256: afcc7ada2718e755169feef1fdc20d70fae8dad57aa8ba4839e34404511bb881
all_ray_symmetric_postprocess_output_sha256: 7c4cf0be6ee2301e11aca8803ac1dd290d607a206faabc82b6d658f807d80d43
all_ray_primary_boundary_output: 05-knowledge/results/lrc14_all_fixed_outsider_ray_primary_boundary_targets_thm4231.out
all_ray_primary_boundary_output_sha256: 4b59d5a66a198b3fc1f4b942b1642d1cae742f8e14ffaa3c4441f75a6d95c195
literal_boundary_census_script: 04-computation/lrc14_literal_boundary_45_cumulative_auditor_thm4231.cpp
literal_boundary_census_output: 05-knowledge/results/lrc14_literal_boundary_45_cumulative_auditor_thm4231.out
literal_boundary_census_script_sha256: bdb8e46c1f5eed18e0fcabf48b20938be4074e236ff8af37ecd5492c6d51357d
literal_boundary_census_output_sha256: a44372cdfc6eba17d25945dc0d1a7de90e0b1bec273dd45f0c3be99552551239
literal_boundary_controls_script: 04-computation/lrc14_literal_boundary_45_cumulative_controls_thm4231.py
literal_boundary_controls_output: 05-knowledge/results/lrc14_literal_boundary_45_cumulative_controls_thm4231.out
literal_boundary_controls_script_sha256: 581194f9cb517a02d3db733989a4ba66429ec4ef4d1fbf26b4029154d51533ca
literal_boundary_controls_output_sha256: fbba123507da124bf7ae345b63aac7827735395b593851175c6825bf8975cd34
literal_boundary_postprocess_script: 04-computation/lrc14_literal_boundary_45_residual_postprocess_thm4231.py
literal_boundary_postprocess_output: 05-knowledge/results/lrc14_literal_boundary_45_residual_postprocess_thm4231.out
literal_boundary_postprocess_script_sha256: f73724902da4b79dc17f7751acfd94df21b104166113ea03165d41ac77eb0e0d
literal_boundary_postprocess_output_sha256: 88539ac6f2f889b8293e9057977c432077d0df9b2dc8500a5ba8cd3773652db5
direct_primary_script: 04-computation/lrc14_two_outsider_direct_body_cofinal_primary_thm4231.cpp
direct_primary_output: 05-knowledge/results/lrc14_two_outsider_direct_body_cofinal_primary_thm4231.out
direct_primary_script_sha256: 48ccab283874ff35f362e6d3c71b36e958dbe0dc444307dadb23b0ec62194bea
direct_primary_output_sha256: 65116c5a9e921ce6cd0a0834baf6749c24e25e922879dc28500b41c754e8f4ca
direct_independent_audit_script: 04-computation/lrc14_two_outsider_direct_body_cofinal_independent_audit_thm4231.cpp
direct_independent_audit_output: 05-knowledge/results/lrc14_two_outsider_direct_body_cofinal_independent_audit_thm4231.out
direct_independent_audit_script_sha256: 933dcb5b03a4ef1ff5a455821968ee2f9d566558c7b64e6575e0079cc64f11bb
direct_independent_audit_output_sha256: b3793a60144fd1f2b04dd7ab7ec51ab1544179c2b99236e168238b0104427adf
literal_control_script: 04-computation/lrc14_two_outsider_direct_body_cofinal_literal_controls_thm4231.py
literal_control_output: 05-knowledge/results/lrc14_two_outsider_direct_body_cofinal_literal_controls_thm4231.out
literal_control_script_sha256: 33414032ad707feaaba988724b655085d14559137361e432ae2d6497c7271172
literal_control_output_sha256: f6a2cb94306c55ca39c43ab8e6973abd9ea9af6cc2f3d08b670566f6475858be
literal_independent_audit_script: 04-computation/lrc14_two_outsider_direct_body_cofinal_literal_independent_audit_thm4231.py
literal_independent_audit_output: 05-knowledge/results/lrc14_two_outsider_direct_body_cofinal_literal_independent_audit_thm4231.out
literal_independent_audit_script_sha256: 9db0bae7ccdcae2609f8a1e237435e2ffd6f092322950d548e98aa2cd6402350
literal_independent_audit_output_sha256: 4102df02b3c35678e8968d60fdf1e306346c0d71097d1aa644eb77d544c9b12b
hybrid_q1290_primary_script: 04-computation/lrc14_hybrid_q1290_cofinal_quadrant_primary_thm4231.py
hybrid_q1290_primary_output: 05-knowledge/results/lrc14_hybrid_q1290_cofinal_quadrant_primary_thm4231.out
hybrid_q1290_primary_script_sha256: 62c41846be5e9b49a3ad67ae744de70dd6cd2f1f45b4c25eafdf731ae13853d3
hybrid_q1290_primary_output_sha256: 9698e9ffc59418fbe7b636ec3a9cf4e73c7ffe5a76f2dfaf48422993f70e3420
hybrid_q1290_independent_audit_script: 04-computation/lrc14_hybrid_q1290_cofinal_quadrant_independent_audit_thm4231.py
hybrid_q1290_independent_audit_output: 05-knowledge/results/lrc14_hybrid_q1290_cofinal_quadrant_independent_audit_thm4231.out
hybrid_q1290_independent_audit_script_sha256: af5bcb4ea1be63464a37c9e48c900b19f30996dc868072224399a42d3175abff
hybrid_q1290_independent_audit_output_sha256: 3b4d00959b58b3416bac5ef95a540ed112ba068a8f4c9d6b343646d1807cc761
fixed_one_primary_script: 04-computation/lrc14_fixed_one_outsider_cofinal_tail_primary_thm4231.cpp
fixed_one_primary_output: 05-knowledge/results/lrc14_fixed_one_outsider_cofinal_tail_primary_thm4231.out
fixed_one_primary_script_sha256: 278426831e51e7052eab980c94274d30cb21255f8a10cc0a54eab18042ce6bf2
fixed_one_primary_output_sha256: dde24be6a541365ab7b0be6a0399aafa544a241d387bcdfb0e31eff0ae08b0c9
fixed_one_literal_script: 04-computation/lrc14_fixed_one_outsider_cofinal_tail_literal_controls_thm4231.py
fixed_one_literal_output: 05-knowledge/results/lrc14_fixed_one_outsider_cofinal_tail_literal_controls_thm4231.out
fixed_one_literal_script_sha256: c973a9141c88d3404c9734c93e2d0f87ba39986e164c044a5de5d0b065cb5dbc
fixed_one_literal_output_sha256: b6ea6eb6bb9dcd21ce5cfe25bc160eaf95b7a6932a7ca12415a4ff8a8864688a
fixed_one_independent_audit_script: 04-computation/lrc14_fixed_one_outsider_cofinal_tail_independent_audit_thm4231.cpp
fixed_one_independent_audit_output: 05-knowledge/results/lrc14_fixed_one_outsider_cofinal_tail_independent_audit_thm4231.out
fixed_one_independent_audit_script_sha256: 1c3508225ae9ee988d39a796886b888921fc37fd5ac0d1153b24a9bac91cbe2a
fixed_one_independent_audit_output_sha256: 83faac9ba0740295f5b9bc179ce908eb0660e47daa9085976352fa609b6e1d04
activation_primary_script: 04-computation/lrc14_arbitrary_pair_cofinal_depth6_primary_thm4231.py
activation_primary_output: 05-knowledge/results/lrc14_arbitrary_pair_cofinal_depth6_primary_thm4231.out
activation_independent_audit_script: 04-computation/lrc14_arbitrary_pair_cofinal_depth6_independent_audit_thm4231.py
activation_independent_audit_output: 05-knowledge/results/lrc14_arbitrary_pair_cofinal_depth6_independent_audit_thm4231.out
activation_primary_script_sha256: 794f0df69956e46c5c73ad6489498b1bd404b9d3643722b129f4b63b092c890a
activation_primary_output_sha256: ddaa2fdd3c822126a5c51c48b526450f34cd22537a40421bba42524ed5c51834
activation_independent_audit_script_sha256: a52740967b84dcef6be68e9ad362cdc51ea9552234c6f5cd8e725a87f0d947f9
activation_independent_audit_output_sha256: 0e7dffe54514f3d46d8076b7dda0c69418ec774c6d47ba95aa2b50bb73f76cec
hash_basis: LF-normalized bytes
audit: >
  PASS / ACCEPT WITH SCOPED INDEPENDENCE. The primary midpoint-cell split-zeta program and a separate
  endpoint-toggle flat-coefficient implementation each exhaust all
  C(30,9)=14,307,150 bodies and recover the same seven-field labelled ledger,
  unique threshold-1307 extremizer, adjacent ceiling slacks, and exact
  activation sidecar. Independent endpoint and midpoint sweeps agree on all
  289 exceptional pairs in the hybrid 1290 quadrant. Independent midpoint and
  endpoint implementations also
  agree on the complete fixed-one-ray census and its unique threshold-543
  certificate extremizer; independent literal controls close that body at
  r=542. A primary split-zeta full-range census and a clean-room truncated
  Boolean-lattice census each exhaust all 18,255,923,400 fixed-ray/body cases
  and agree on positivity, the unique oriented maximum 931, and the normalized
  minimum. The clean-room path freezes the complete K(q) table and its
  181,242-edge two-orientation residual; primary reruns independently confirm
  its load-bearing q=744,824,825 boundary. Full-table digest agreement is not
  claimed. A cumulative joint-wall census exhausts all 45 certificate edges
  with endpoint at least 770 and finds all 643,821,750 cases strictly safe;
  `-O2/-O3` builds byte-agree, and a separate direct midpoint scan confirms
  all 45 unique extremizers. The old activation artifacts independently
  retain the 17547/17548 cover transition.
---

# THM-4231 -- literal-patched all-ray cofinite pair reduction and sharp Haar sidecars

**PROVED RELATIVE TO THM-4150/4156/4170/4191/4207/4211 + VERIFIED-EXACT + TARGETED
INDEPENDENT BOUNDARY AUDIT; LRC(14) REMAINS OPEN.**

## 1. Statement

Retain the THM-4156 pool

```text
P={8,10,15,16,20,30,40,42,60,63,80,84,85,88,95,
   120,126,132,143,145,168,170,176,190,193,240,252,
   264,286,290}.                                        (1)
```

For a finite positive label set `S`, put

```text
G_S={y in R/Z:min_(s in S)||sy||>=1/14},
alpha=4/63.                                              (2)
```

> **Cofinal arbitrary-pair theorem.** For every two distinct integers
> `q,r>=1290` and every `B in binom(P,9)`,
>
> ```text
> mu(G_(B union {q,r}))>=alpha.                          (3)
> ```

More strongly:

> **Cofinite arbitrary-pair theorem.** For every two distinct outsiders
> `q,r notin P` with `max(q,r)>=770` and every `B in binom(P,9)`,
>
> ```text
> mu(G_(B union {q,r}))>=alpha.                          (3a)
> ```

Equation `(3)` is retained as an independently literal-audited hybrid
subtheorem. Equation `(3a)` reduces the combined certificate-plus-literal
proof remainder to exactly `181,194` finite outsider pairs with maximum at
most `769`; it does not assert that any residual pair is unsafe, prove a
minimal literal threshold, give physical entry of an arbitrary LRC row into
this pool, or prove LRC(14).

## 2. Inheritance and corrected object

The closest analytic mechanism is THM-4170's discrepancy estimate for a
finite union of circle intervals. The canonical hostile is THM-4207's failure
of marginal repair-deck composition. The corrected near miss is to keep
sharpening deletion certificates before testing the target body itself.
MISTAKE-520 exposed the same error one level earlier: ambient full-pool safety
already implied all fixed-pool subsets by heredity.

The least-used sidecar is the exact component count of the literal body safe
set `G_B`. Once `(M_B,c_B)` is retained for every labelled body, no deletion
repair is needed at the proved direct/hybrid cutoff. The old depth-six repair deck
remains an independently useful incidence and activation sidecar.

## 3. Direct-body Bonferroni theorem

For `B in binom(P,9)`, set

```text
U_B=G_B,          M_B=mu(U_B),          c_B=#components(U_B). (4)
```

THM-4170 equation `(9)` applies to every union of `c_B` circle intervals:

```text
mu(U_B intersect G_s)>=(6/7)M_B-6c_B/(49s).             (5)
```

Bonferroni inside the fixed set `U_B` gives

```text
mu(U_B intersect G_q intersect G_r)
 >=(5/7)M_B-(6c_B/49)(1/q+1/r).                         (6)
```

Define

```text
s_B=45M_B-4.                                            (7)
```

The exact census below proves `s_B>0` for every labelled nine-body. Put

```text
kappa_body(B)=ceil(108c_B/(7s_B)).                      (8)
```

If `q,r>=kappa_body(B)`, then

```text
(6c_B/49)(1/q+1/r)<=12c_B/(49kappa_body(B))<=s_B/63.
```

Since `(5/7)M_B-alpha=s_B/63`, equation `(6)` proves the `q,r>=1307`
subcase as soon as

```text
max_(B in binom(P,9)) kappa_body(B)<=1307.              (9)
```

No independence or equidistribution of `G_q` and `G_r` is assumed. This direct
argument first proves the cofinal theorem with `1307`; Section 6 supplies the
finite exact exceptional ledger needed to lower the final statement to `1290`.

## 4. Complete exact body census

The fixed-pool arrangement has common denominator

```text
D=18,241,159,416,480,
7,134 walls and 7,133 open cells.                       (10)
```

The primary exhausts every one of

```text
binom(30,9)=14,307,150                                  (11)
```

labelled bodies. All `s_B` are strictly positive. The activation range is
`101` through `1307`; the maximum is unique:

```text
B*={170,176,190,193,240,252,264,286,290},
M_B* ticks=4,579,301,272,924,
c_B*=618,
s_B* ticks=133,103,919,615,660,
kappa_body(B*)=1307.                                    (12)
```

The second-largest threshold is `1290`, also attained once. The body-mass and
component ranges are

```text
4,131,630,720,966 <= M_B ticks <= 7,498,291,720,920,
80 <= c_B <= 618.                                       (13)
```

For each body in increasing thirty-bit mask order, the primary freezes

```text
(body,kappa_body,mass_ticks,components,s_B ticks,
 first_disjoint_repair_activation,composite_threshold). (14)
```

The bytewise FNV-1a-64 fingerprint of the complete seven-field ledger is

```text
5b8c08ad9d02622a.                                       (15)
```

### Split-zeta construction

Let `F_i` be the failed-label mask on cell `i`, `ell_i` its exact length, and
index cells cyclically. Then

```text
M_B=sum_i ell_i 1_(F_i intersect B=empty),
c_B=sum_i [1_(F_i intersect B=empty)
           -1_((F_(i-1) union F_i) intersect B=empty)]. (16)
```

The primary midpoint classifier aggregates `(16)` to `2,939` mass and `1,457`
signed component-boundary coefficients. It splits the thirty labels into two
fifteen-label halves, performs a complete fifteen-bit subset-zeta transform on
one half for every filtered choice on the other, and visits each body exactly
once.

## 5. Independent exact audit

The referee imports neither the midpoint classifier nor the primary split
schedule. It constructs the arrangement by endpoint enter/leave toggles,
builds a flat disjoint-coefficient census, and independently enumerates all
`14,307,150` bodies. It recovers

```text
cells=7,133,
effective mass/component coefficients=2,939/1,457,
all-positive surpluses=14,307,150,
unique maximum=1307 at B*,
unique second maximum=1290,
unique minimum=101,
FNV ledger=5b8c08ad9d02622a.                              (17)
```

It also reconstructs all `140,082` strict depth-six repair rows and agrees on
the first disjoint repair activation of every body. Two additional commutative
digest controls and exact left/right ceiling inequalities are frozen in its
output.

## 6. Sharp direct boundary and hybrid literal quadrant

At `Q=1306`, exactly `14,307,149` bodies pass the direct test and `B*` is the
unique failure. At `Q=1307`, all bodies pass. The exact symmetric gaps have
opposite signs, so `1307` is the minimum of the uniform direct certificate.
The unique second-largest body threshold is `1290`; hence every `B!=B*` is
already direct-safe for distinct `q,r>=1290`.

It remains to patch `B*`. Order the outsiders as `q<r`. If `q>=1307`, the
direct theorem applies. For every `1290<=q<=1306`, exact arithmetic proves
that the first direct-safe larger label is

```text
r_0(q)=2614-q.                                          (18)
```

Both the strict failure at `r_0(q)-1` and success at `r_0(q)` are checked by
integer cross-multiplication. The unresolved literal region therefore has

```text
sum_(q=1290)^1306 (2614-2q-1)=33+31+...+1=289 pairs.   (18a)
```

An endpoint-toggle sweep and a separately written wall/midpoint classifier
exhaust all 289 pairs. Both report the same semantic ledger hash
`ab05deaad2b90c981c15abd2420a74ae74ba85e303c98f825bd3548d7db01860`,
zero nonpositive rows, and the unique smallest literal margin

```text
(q,r)=(1300,1305),
63mu-4=18886235531/2585198330>0.                        (18b)
```

Thus `B*` is safe throughout the missing triangle, and `(3)` holds for every
distinct `q,r>=1290`. As a representative hostile control,
`B* union {1306,1307}` is strictly safe even though its direct bound is
negative. The number `1290` is a proved hybrid cutoff, not a claimed minimal
literal threshold; `1307` remains the exact direct-certificate cutoff.

## 7. Superseded depth-six repair headline

For comparison, retain the exact strict depth-six repair activation

```text
rho(B)=min_{R intersect B=empty} kappa_2(R),             (19)
```

with `rho(B)=infinity` if there is no strict repair. The previous version of
this theorem proved

```text
|A_6^(2)(17547)|=54,563 with one nine-cover W,
|A_6^(2)(17548)|=54,566 with no nine-cover.              (20)
```

The unique cover was

```text
W={85,88,143,168,193,240,252,264,290}.                  (21)
```

Both old primary/referee artifacts remain canonical and independently audited.
The direct census now gives

```text
M_W ticks=4,802,564,195,362,
c_W=506,
kappa_body(W)=995,
rho(W)=17548.                                            (22)
```

Thus the cover was itself uniformly direct-safe long before its last repair
activated. If

```text
theta(B)=min(kappa_body(B),rho(B)),                      (23)
```

then `max_B theta(B)=1307`. In fact the first strict depth-six repair activates
only at `3077`, so the repair branch is dormant at the optimum. Equation `(20)`
is still an exact sharp activation-filter transition, but its `17548` safety
headline is superseded by `(3)`.

## 8. Exact outsider lift and odd-tail transfer

Fix distinct `q,r>=1290`. Every eleven-subset of

```text
P union {q,r}                                            (24)
```

is safe:

- with no outsider, THM-4156 proves the stronger full-pool bound and
  monotonicity handles every subset;
- with one labelled outsider, THM-4191 handles every ten-subset of `P`; and
- with both outsiders, `(3)` handles every nine-subset of `P`.

The exact labelled partition is

```text
binom(30,11)+2binom(30,10)+binom(30,9)
 =54,627,300+2(30,045,015)+14,307,150
 =129,024,480=binom(32,11).                              (25)
```

For every two-outsider body, every positive integer `c`, and every two
distinct positive odd integers `a,b`, THM-4150 supplies some `x in R/Z` with

```text
min_(v in 2c(B union {q,r}) union {a,b})||vx||>=1/14.   (26)
```

These are genuine infinite LRC(14) families of thirteen relative speeds, not
a proof for an arbitrary instance.

## 9. Independently audited fixed-one ray

The same retained-coordinate method is much stronger once one outsider is
fixed literally. For `B in binom(P,9)`, put

```text
V_B=G_(B union {1}),  M_B^(1)=mu(V_B),  c_B^(1)=#components(V_B). (26a)
```

THM-4170 gives

```text
mu(V_B intersect G_r)>=(6/7)M_B^(1)-6c_B^(1)/(49r).    (26b)
```

The complete census proves `54M_B^(1)-4>0` for all bodies and, with

```text
kappa_1(B)=ceil(54c_B^(1)/(7(54M_B^(1)-4))),           (26c)
```

finds the unique maximum

```text
max_B kappa_1(B)=543,
B*={170,176,190,193,240,252,264,286,290},
M_B*^(1) ticks=3,885,436,686,322,
c_B*^(1)=528,
(54M_B*^(1)-4) ticks=136,848,943,395,468.              (26d)
```

Therefore every `B union {1,r}` is direct-safe for every `r>=543`. The primary
midpoint/split-zeta census and an independent endpoint-event/down-set census
each exhaust all `14,307,150` bodies and recover the same extremizer, mass,
component count, surplus, cutoff, and unique second cutoff `530`. Their
coefficient aggregations and digests are independently frozen.

At `r=542`, every body but `B*` passes the analytic certificate. The literal
endpoint controls nevertheless give

```text
mu(G_(B* union {1,542}))=689745548341/3874102039080,
63mu-4=443770815701/61493683160>0.                      (26e)
```

The unique second certificate cutoff is `530`, so every body other than `B*`
is already direct-safe at `r=542`; equation `(26e)` supplies the missing
literal `B*` case. Thus the whole fixed-one ray is safe for every `r>=542`.
The value `543` is certificate-sharp, while `542` is not claimed
literal-minimal. With the zero/one layers, this also fills all `C(32,11)`
faces of every chart `P union {1,r}` for `r>=542`, and THM-4150 again supplies
the odd-tail LRC(14) families.

## 10. Complete fixed-ray atlas and exact certificate residual

For every

```text
q in {1,...,1306}\P,        B in binom(P,9),             (27)
```

put

```text
V_(q,B)=G_(B union {q}),
M_(q,B)=m_(q,B)/D_q,
c_(q,B)=#components(V_(q,B)),
s_(q,B)=54m_(q,B)-4D_q,
K_q(B)=ceil(54c_(q,B)D_q/(7s_(q,B))).                   (28)
```

THM-4170 gives

```text
mu(V_(q,B) intersect G_r)
 >=(6/7)M_(q,B)-6c_(q,B)/(49r),                         (29)
```

so `r>=K_q(B)` proves the target. Define `K(q)=max_B K_q(B)`. A primary
split-zeta census and a clean-room endpoint-event/truncated-Mobius census each
exhaust

```text
1276 binom(30,9)=18,255,923,400                         (30)
```

ray/body cases. Every `s_(q,B)` is positive. Both paths agree on the unique
oriented maximum

```text
max_q K(q)=931,
q=1305,
B={20,170,190,193,240,252,264,286,290},                 (31)
```

and on the normalized surplus minimum, attained at `q=50` and
`B={8,15,80,84,88,95,120,145,170}`. The clean-room ledger freezes

```text
full XOR/SUM=e2b35ef11bbc30a7/da92ed0c54582151,
FNV(q,K(q),least maximizing body)=e5b533495b4d0d6f.    (32)
```

Full-table digest agreement with the primary is not claimed. The primary
transcript independently agrees on the full case universe, positivity, `(31)`,
and the normalized minimum, and fresh primary replays check the load-bearing
values

```text
K(744)=825,             K(824)=781,             K(825)=711. (33)
```

For outsiders `a<b`, the two orientations certify the pair whenever

```text
b>=K(a)  or  a>=K(b).                                  (34)
```

The exact residual graph of this certificate therefore has an edge precisely
when both strict reverse inequalities hold. Exact postprocessing of the frozen
`K(q)` table gives

```text
|E_res|=181,242,
FNV(E_res)=8a4e1370fb023907,
max endpoint=824,
unique edge at endpoint 824={744,824}.                  (35)
```

It follows that for every two distinct outsiders `a,b notin P`, every
`B in binom(P,9)`, and

```text
max(a,b)>=825,                                          (36)
```

one has `mu(G_(B union {a,b}))>=4/63`. Indeed, if both labels are at most
`1306`, `(34)--(35)` applies. If `a<b`, `a<=1306<b`, then
`b>931>=K(a)`. If `1307<=a<b`, the hybrid quadrant `(3)` applies. Section 8
therefore fills every `C(32,11)` face in each chart `P union {a,b}` satisfying
`(36)`, and THM-4150 supplies the corresponding odd-tail families.

The value `825` is the exact cutoff of this two-orientation fixed-ray
certificate because `{744,824}` survives both inequalities. It is not claimed
literal-minimal, and no residual pair is asserted unsafe.

## 11. Literal boundary peel and combined residual

Inside the certificate graph `(35)`, put

```text
L={ (a,b) in E_res : b>=770 }.                          (36a)
```

Exact postprocessing gives `|L|=45`. A fresh cumulative joint-wall census
builds the literal geometry for every pair in `L` and exhausts

```text
45 binom(30,9)=643,821,750                              (36b)
```

pair/body cases. Every case has `63mu-4>0`; there are no zero or negative
cases, and every pair has a unique minimizing body. The complete labelled
ledger has

```text
XOR/SUM=fb8c4e194710163a/1990d28492857d56,
ordered FNV=029d6427584abd92.                           (36c)
```

Independent `-O2/-O3` builds byte-agree. A separate direct midpoint-wall
scanner, which uses neither the Boolean coefficient table nor its Mobius
transform, reproduces the body, component count, and reduced margin of all 45
minimizers.

Three further full-pool closures already lie in proved canon:

```text
H={{1,542},{49,50},{50,51}},                            (36d)
```

from Section 9, THM-4211 equation `(7)`, and THM-4207 equation `(5)`,
respectively. All 48 pairs in `L union H` are distinct members of `E_res`.
Removing them gives the exact combined proof residual

```text
|E_rem|=181,194,
FNV(E_rem)=3874fecac4ecbd8a,
max endpoint=769,
endpoint-769 layer={{616,769},{721,769}}.               (36e)
```

Therefore `(3a)` holds. For `a<b<=1306`, either `(34)` certifies the pair or
it belongs to `E_res`; if also `b>=770`, `(36a)--(36c)` closes it literally.
If `a<=1306<b`, then `b>931>=K(a)`. If `1307<=a<b`, the hybrid quadrant `(3)`
applies. Section 8 then fills every `C(32,11)` face in each such chart, and
THM-4150 supplies the corresponding odd-tail families.

The value `770` is exact only for this certificate-plus-literal proof graph:
the two pairs in the last line of `(36e)` remain unclosed by it. Neither pair
is asserted unsafe, and `770` is not claimed literal-minimal.

## 12. Pair-plane and method boundaries

The unbounded two-outsider problem is now reduced to the exact finite label
graph `(36e)`, not merely to a family of rays. Literal or stronger
pair-correlated treatment of its `181,194` edges remains open.

THM-4227's scale wedge, THM-4228's common-gcd region, and THM-4233's displayed
primitive families are coverage-subsumed by `(3a)`. Their scale order, common
divisor, primitive shape, and cyclotomic-zero methods remain useful inside the
finite residual. THM-4234 supplies complementary restricted positive controls
on the residual `q=50` chart; it is not a full-pool closure.

## 13. Boolean sidecar and connection contract

Equation `(16)` is a lower Boolean zeta transform plus a cyclic transition
count. This is the lawful transfer from THM-4223's owner-refined Boolean
hierarchy: retain the exact label tensor and the adjacency required by the
next operation. There is no intrinsic pairwise orientation on the cells, so no
tournament is imposed.

```text
source:       fixed-q cyclic failure masks and adjacent mask unions
target:       every labelled (q,B), its K(q), and the patched residual graph
map:          disjointness zeta -> (M_(q,B),c_(q,B)) -> K(q) -> E_res -> E_rem
preserved:    q/body labels, exact mass, components, both pair orientations
destroyed:    component addresses and literal q/r phase alignment
sidecar:      adjacent mask unions and literal endpoint sweep
hostile:      {616,769},{721,769}, unclosed but not asserted unsafe
decisive test: max endpoint and exact digest of E_rem.     (37)
```

## 14. Reproduction

From the repository root, compile both C++ paths with C++20, optimization, and
OpenMP, then compare stdout byte-for-byte with their frozen outputs:

```bash
g++ -O3 -fopenmp -std=c++20 \
  04-computation/lrc14_two_outsider_direct_body_cofinal_primary_thm4231.cpp \
  -o /tmp/thm4231-primary
OMP_NUM_THREADS=12 /tmp/thm4231-primary

g++ -O3 -fopenmp -std=c++20 \
  04-computation/lrc14_two_outsider_direct_body_cofinal_independent_audit_thm4231.cpp \
  -o /tmp/thm4231-referee
OMP_NUM_THREADS=12 /tmp/thm4231-referee 1306

python3 -B \
  04-computation/lrc14_two_outsider_direct_body_cofinal_literal_controls_thm4231.py

python3 -B \
  04-computation/lrc14_two_outsider_direct_body_cofinal_literal_independent_audit_thm4231.py

python3 -B \
  04-computation/lrc14_hybrid_q1290_cofinal_quadrant_primary_thm4231.py

python3 -B \
  04-computation/lrc14_hybrid_q1290_cofinal_quadrant_independent_audit_thm4231.py

g++ -O3 -fopenmp -std=c++20 \
  04-computation/lrc14_fixed_one_outsider_cofinal_tail_primary_thm4231.cpp \
  -o /tmp/thm4231-q1-primary
OMP_NUM_THREADS=12 /tmp/thm4231-q1-primary 1

g++ -O3 -fopenmp -std=c++20 \
  04-computation/lrc14_fixed_one_outsider_cofinal_tail_independent_audit_thm4231.cpp \
  -o /tmp/thm4231-q1-referee
OMP_NUM_THREADS=12 /tmp/thm4231-q1-referee

python3 -B \
  04-computation/lrc14_fixed_one_outsider_cofinal_tail_literal_controls_thm4231.py

g++ -O3 -fopenmp -std=c++20 \
  04-computation/lrc14_all_fixed_outsider_ray_census_primary_thm4231.cpp \
  -o /tmp/thm4231-all-q-primary
OMP_NUM_THREADS=12 /tmp/thm4231-all-q-primary 1 1306

g++ -O3 -fopenmp -std=c++20 -Wall -Wextra -Werror \
  04-computation/lrc14_all_fixed_outsider_ray_census_independent_audit_thm4231.cpp \
  -o /tmp/thm4231-all-q-referee
OMP_NUM_THREADS=12 /tmp/thm4231-all-q-referee

python3 -B \
  04-computation/lrc14_all_fixed_outsider_ray_symmetric_postprocess_thm4231.py

OMP_NUM_THREADS=12 /tmp/thm4231-q1-primary 744
OMP_NUM_THREADS=12 /tmp/thm4231-q1-primary 824
OMP_NUM_THREADS=12 /tmp/thm4231-q1-primary 825

g++ -O3 -fopenmp -std=c++20 -Wall -Wextra -Werror \
  04-computation/lrc14_literal_boundary_45_cumulative_auditor_thm4231.cpp \
  -o /tmp/thm4231-literal-boundary
OMP_NUM_THREADS=12 /tmp/thm4231-literal-boundary

python3 -B \
  04-computation/lrc14_literal_boundary_45_cumulative_controls_thm4231.py

python3 -B \
  04-computation/lrc14_literal_boundary_45_residual_postprocess_thm4231.py
```

The frozen referee transcript is the `1306` invocation above; the no-argument
mode performs a binary-search diagnostic and emits additional lines. The
primary output is stable under `-O2/-O3` and one/twelve threads. The all-ray
referee uses a structurally different endpoint geometry and truncated
Boolean-lattice census; its full run is intentionally expensive. The
postprocessor freezes the complete `K(q)` table and exact two-orientation
graph. The three final primary invocations reproduce the sharp boundary
controls in their concatenated frozen output. The cumulative literal source
imports the clean-room geometry engine but rebuilds each two-outsider joint
wall arrangement; its output covers exactly the 45 edges in `(36a)`. The
direct-control script independently midpoint-scans all 45 extremal eleven-
speed sets, and the final postprocessor verifies `(36d)--(36e)`. The older
Python activation paths reproduce the exact sidecar `(20)`.

**QED.**
