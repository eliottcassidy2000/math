---
id: THM-3910
title: "LRC(14) auxiliary-center erosion and t-sheet variance response"
status: >
  PROVED RELATIVE TO CITED LRCUpTo13 + VERIFIED-EXACT + INDEPENDENTLY
  HOSTILE-AUDITED.  In THM-3878's conditional t>=U 11+2 slice, a fixed-radius
  auxiliary-center response closes 41 of the 57 scale-one certificate
  survivors, leaving 16 scale-one types and the scale-two (1,9) type: 17
  total.  A compact-to-open body-component response and an exact integer
  t-sheet/Fourier carrier give further sufficient filters; the former closes
  all 58 AP11 controls.  Pairwise tournament/Gram data are proved insufficient
  by a native AP11 cubic-response collision.  The t<U slice and LRC(14) remain
  open.
source: root / THM-3878 response audit, 2026-08-23
audit: >
  INDEPENDENT HOSTILE AUDIT PASS.  Two independent rational open-component
  constructions agree on the strict radius-1/182 erosion, all 398 finite
  auxiliary cells, the 41/16 split, the unique a=p killer law, both
  erosion-boundary contact controls, and the scale-two nonclosure.  A
  standalone wall/topology audit confirms the compact-to-open law, AP11's
  14 positive arcs plus four isolated walls, the owner-labelled 1/77 arc,
  the complete beta staircase, and both q<=25 hostiles.  Exact sheet-count
  companions verify 4,214 AP11 scale-one rows, the CV and integer tariffs,
  boundary scales and scale-two controls.  An independent wall-cell engine
  reproduces the native pairwise-moment collision.  Every stored transcript
  byte-matches normal and optimized Python.  A separate proof audit checked
  types, strict inequalities, duplicate-speed handling, scale-two weights,
  and every equality condition.
depends_on:
  - THM-3878-lrc14-eleven-plus-two-harmonic-absorption-seam-collapse
  - LRCUpTo13
related:
  - THM-579-lrc14-covering-floor-sheetcount-CV-criterion
  - THM-1042-component-length-obstruction-to-additive-certificates
  - THM-2221-tournament-context-cut-metric-and-pinned-transport-response
  - THM-3377-path-colour-deletion-compiler-and-skew-current
  - THM-3380-hamiltonian-deletion-layer-monoid-semiring-and-small-order-boundaries
  - THM-3729-rooted-pfaffian-response-and-sign-root-deletion-average
  - THM-3961-lrc14-signed-endpoint-cross-phase-and-fixed-scale-two-family
script: 04-computation/lrc14_auxiliary_center_erosion_response_thm3910.py
output: 05-knowledge/results/lrc14_auxiliary_center_erosion_response_thm3910.out
script_sha256: d4c6e9598cc30ca1155d8e407e3065aa2034556a8ac04349890d1baa6796f261
output_sha256: 9e0eea75a04e1ad4963a91c1d19c4eb73bb1d69b26cdaecc4230b6e91713ffb8
semantic_sha256: 10ce5bf6777272db163be6954f9f6fc54119eb8ae38710b261113789c14a0e1e
independent_audit_script: 04-computation/lrc14_auxiliary_center_erosion_independent_audit_thm3910.py
independent_audit_output: 05-knowledge/results/lrc14_auxiliary_center_erosion_independent_audit_thm3910.out
independent_audit_script_sha256: f88974b156f26001b93378160aa200eae4118d1c26e7c10859698e0f74d1acf6
independent_audit_output_sha256: c12658316a45ae13f1d1df0bb185d7535395aaecc0f20211fd66587f438de36d
independent_audit_semantic_sha256: fa8e359289f0a2468b63df89d8335c4717bebad43cc5afe1545ba982c93d6f83
component_script: 04-computation/lrc14_safe_component_response_thm3910.py
component_output: 05-knowledge/results/lrc14_safe_component_response_thm3910.out
component_script_sha256: 2a2c9da20e79e21f471ecf4d7d8fd443798923d52ed3f781844a618a48d01d3b
component_output_sha256: fa87eb80629da11220efe3345736be58bcf6cebba7c979d4cb4068b643488222
component_semantic_sha256: 01fcb7b28f890a56648314611dd79ee740f9f47e9426d6f616c126e4a352e914
component_independent_script: 04-computation/lrc14_safe_component_response_independent_audit_thm3910.py
component_independent_output: 05-knowledge/results/lrc14_safe_component_response_independent_audit_thm3910.out
component_independent_script_sha256: f4e359a27ab31d9dd8bf31477415ac02c0678428c8c2a709da649c48b5763a19
component_independent_output_sha256: 689a04e17b162c8b62b3d4d62e6a4b046b8cb517f93bdd4e3073f92c2213655c
component_independent_semantic_sha256: f6f5bcb852da0efa5072a53cd211e1cc57974f623e0cabfae9edaea5f35853e9
lift_count_script: 04-computation/lrc14_t_sheet_lift_count_variance_thm3910.py
lift_count_output: 05-knowledge/results/lrc14_t_sheet_lift_count_variance_thm3910.out
lift_count_script_sha256: a8663fee9f6294560102df3ba8bcdef67f07d6aafeddcb3ec5c57faab08307a8
lift_count_output_sha256: b1d2901f5e1d43142f1cbb4fb97648403ef52073b92683c93a9ff8da06b3f279
lift_count_semantic_sha256: 770436c3cac664b0c5b954e92f64c153d9739a099ad63511d608117055cdd821
lift_boundary_script: 04-computation/lrc14_t_sheet_lift_count_ap11_boundary_thm3910.py
lift_boundary_output: 05-knowledge/results/lrc14_t_sheet_lift_count_ap11_boundary_thm3910.out
lift_boundary_script_sha256: e6c48b6ac9307278a86d44e65f136d8921652dbd27e5006ce7967097fcda092e
lift_boundary_output_sha256: d468a2c29c83c0812f88b3cce60b49f9594ee9a772c53ee71c5e11ebd7263924
lift_boundary_semantic_sha256: e95e34b5206949c14b6757c5e8b36e2e27284cca74f983ce6126b548f00e0e68
scale2_script: 04-computation/lrc14_t_sheet_lift_count_scale2_thm3910.py
scale2_output: 05-knowledge/results/lrc14_t_sheet_lift_count_scale2_thm3910.out
scale2_script_sha256: 544b4e0a023d0b17ed244aa1d33b1086cb49626fd0258c84516dcfe710283d85
scale2_output_sha256: 909d4cc2be5479c78e730a9cc3f26c744612870a6be95305506ece24b7fe562c
scale2_semantic_sha256: 749fc3c20341ce4c199d4e07f34f5c477b55194ebc049af1fc82f008e414f251
tournament_script: 04-computation/lrc14_tournament_response_exact_scout_thm3910.py
tournament_output: 05-knowledge/results/lrc14_tournament_response_exact_scout_thm3910.out
tournament_script_sha256: ec3faea956d5d6c9f87b84e68ac1bead08d821fd66a1452b22592185a829b2cd
tournament_output_sha256: d503ecab39631d03ef33aadafdcb0e805ca0e1dc02053730e5d781fb9338038e
tournament_semantic_sha256: 0ea2b13f5bee6318ca5872189a867e19c760111c58c33d0d047a357fe5e004ac
tournament_independent_script: 04-computation/lrc14_tournament_response_collision_independent_audit_thm3910.py
tournament_independent_output: 05-knowledge/results/lrc14_tournament_response_collision_independent_audit_thm3910.out
tournament_independent_script_sha256: a038be527cd4d435c102e4b8c09f998767b86310e9a7c686e7ecd033dd38d040
tournament_independent_output_sha256: ec8bcaf8457f70e6d5985f5f09bfdaa8ab51a5e1af4db6f4bbbda1eec8554acb
hash_basis: raw LF bytes
---

# THM-3910 -- auxiliary-center erosion and the t-sheet response

**PROVED RELATIVE TO CITED LRCUpTo13 + VERIFIED-EXACT + INDEPENDENTLY
HOSTILE-AUDITED.**  This theorem acts only inside THM-3878's exact `11+2`
branch.  Its strongest census is conditional on `t>=U`; it says nothing new
about `t<U` or other THM-3818 component shapes.  In particular,
**LRC(14) remains OPEN.**

## 1. Setup

Work on `T=R/Z` with normalized Haar measure.  Put

```text
delta=1/14,
D_v={x in T: ||v x||<delta},
G_u=intersection_i D_(u_i)^c,
U=max_i u_i,
mu=meas(G_u).
```

Here `u` is an eleven-speed component.  In THM-3878's scale-one branch the
two-speed obstruction is

```text
A_pq=D_p union D_q,
beta_1(p,q)=largest open component length of A_pq.       (1)
```

For odd `p,q` in scale two, let `C_pq` be THM-3878's `w=2z` image of
`A_pq intersect (A_pq-1/2)`, and let `beta_2(p,q)` be its largest open
component length.

## 2. Fixed-radius auxiliary-center response

Assume `t>=U`.  For every positive integer `a`, cited `LRCUpTo13` applied to
the at-most-twelve distinct speeds in `u union {a t}` gives a time `y_0` at
which all retained speeds have clearance at least `1/13`.  Repeated speeds
are deduplicated; the resulting smaller pack has an at least as strong cited
clearance.  Put `w_0=t y_0`.  Then

```text
||a w_0||>=1/13.                                        (2)
```

If `|eta|<=1/182`, the time `y_0+eta/t` changes every body phase by at most

```text
U|eta|/t<=1/182=1/13-1/14.
```

Consequently the body remains `1/14`-safe on the **closed** `w`-arc

```text
w_0+[-1/182,1/182].                                    (3)
```

The auxiliary speed is used only to locate the center; it need not remain
safe on all of (3).  If the full scale-one row were a counterexample, every
point of (3) would belong to the **open** set `A_pq`.  Thus a necessary
condition is

```text
there exists w_0 with ||a w_0||>=1/13 and
w_0+[-1/182,1/182] subset A_pq.                         (4)
```

The same implication holds in scale two with `C_pq`: if the original row is
bad, both physical lifts are pair-bad at every body-safe quotient time.

If the largest obstruction component has length `beta`, its largest eroded
core has length `gamma=beta-1/91`.  Every component of
`{w:||a w||<1/13}` has length `2/(13a)`.  Therefore

```text
a gamma>2/13                                             (5)
```

forces the largest erosion core to contain an `a`-deep point and makes the
multiplier search finite.

Two independent exact component engines give

```text
scale-one input                                      57
closed by fixed-radius auxiliary response            41
remaining                                             16
finite multiplier cells                              398
per-row terminal cutoffs                            2..12
```

Every closed pair has the unique killing multiplier `a=p`.  Grouped by `p`,
the 41 closures are

```text
p=2: q=15,21,23
p=3: q=14,17,19,20,22,26,31,38
p=4: q=13,19,21,25,37,43,49,51
p=5: q=17,18,24,29,36,39,41,42,48,53,54,63
p=6: q=17,19,23,41,47,53,65
p=7: q=13,15,22.
```

The exact remaining scale-one list is

```text
(1,3),(1,4),(1,9),(1,10),(2,3),(2,9),(3,7),(3,8),
(4,7),(5,6),(5,12),(6,11),(7,10),(8,9),(8,21),(9,11).  (6)
```

For the scale-two `(1,9)` obstruction, no multiplier kills the necessary
center condition.  Hence THM-3878's conditional `58`-type ledger becomes

```text
16 scale-one types plus the scale-two (1,9) type = 17. (7)
```

Strict open/closed semantics are load-bearing.  At `(p,q,a)=(4,13,4)` and
`(7,13,7)`, the deep set has only erosion-boundary contacts; the closed source
arc would touch the complement of the open obstruction, so those contacts do
not survive (4).

## 3. Compact-to-open body-component response

Let `lambda_+(u)` be the largest positive length of a closed connected
component of `G_u`; isolated safe walls contribute zero.  A counterexample in
the conditional branch necessarily satisfies

```text
scale one:  t lambda_+(u)<beta_1(p,q),
scale two:  t lambda_+(u)<beta_2(p,q).                  (8)
```

Map a longest closed body-safe arc by `y -> t y`.  If its image wraps the
circle, containment in the proper open obstruction is impossible.  Otherwise
its image is a compact connected arc of exact length `t lambda_+(u)`.  A
compact arc contained in an open component has strictly smaller length than
that component.  This proves (8), including equality.

For AP11, `u=(1,...,11)`, the corrected THM-1042 topology contains

```text
J=[1/14,13/154],       |J|=1/77.                        (9)
```

The left endpoint has unique equality owner `1`, the right endpoint unique
owner `11`, and the interior is strict.  Since `U=11`, `U|J|=1/7`.
On the 57-pair bank, `beta_1<=1/7`, with equality exactly at

```text
(1,3),(1,4),(1,9),(1,10).
```

Thus (8), including its strict equality boundary, closes all 57 AP11
scale-one controls for every `t>=11`.  It also closes AP11 scale two because
`beta_2(1,9)=2/63<1/7`.

The exact cumulative response staircase is

```text
U lambda_+(u)>=1/35  closes 15 of 57,
U lambda_+(u)>=1/28  closes 34 of 57,
U lambda_+(u)>=1/21  closes 46 of 57,
U lambda_+(u)>=1/14  closes 52 of 57,
U lambda_+(u)>=1/7   closes 57 of 57.                  (10)
```

The generic body witness used in THM-3878 supplies only
`U lambda>=1/42`, which closes none of the 57.  Formula (8) is a body-shape
sidecar rather than a restatement of the old cyclic bound.

## 4. The exact integer t-sheet carrier

Let `f=1_(G_u)` and define

```text
N_t(w)=sum_(a=0)^(t-1) f((w+a)/t),
m=integral N_t=t mu.                                  (11)
```

For scale one put `H=A_pq^c`, `h=meas(H)`, `b=1-h`.  Changing variables on
the `t` inverse branches gives

```text
meas(G_u \ (D_(tp) union D_(tq)))
  =(1/t) integral_H N_t(w) dw.                         (12)
```

Thus THM-3878's missing conditioned one- and two-danger moments are labelled
moments of one integer carrier.

With `fhat(n)=integral f(x) exp(-2 pi i n x) dx`, one has

```text
Nhat_t(k)=t fhat(k t),
Var(N_t)=integral (N_t-m)^2
        =t^2 sum_(k!=0)|fhat(k t)|^2.                  (13)
```

Hence `Var(N_t)/m^2 ->0` by the `l2` Fourier tail.  Every fixed
positive-measure body and fixed pair type is lonely for all sufficiently
large `t`; this is nonuniform in the body.

## 5. Variance and the sharp integer failure tariff

If (12) vanishes, `N_t=0` almost everywhere on `H`.  Cauchy--Schwarz gives

```text
Var(N_t)>=m^2 h/b,
CV(N_t)^2<h/b  implies positive safe mass.             (14)
```

Put `k=floor(m/b)` and `theta=m/b-k`.  Convexity on the nonnegative integers
sharpens the necessary failure tariff to

```text
Var(N_t)>=m^2 h/b+b theta(1-theta).                    (15)
```

Equality in the abstract integer-response problem holds exactly when
`N_t=0` on `H`, while on the obstruction it takes values `k,k+1` in relative
proportions `1-theta,theta`.  Realizable lift counts have extra geometric
constraints, so this is not a converse construction.

If `G_u` has `r` positive-length components, its indicator has total
variation `2r`, and

```text
|fhat(n)|<=r/(pi |n|),
Var(N_t)<=r^2/3.                                       (16)
```

Thus the explicit sufficient condition is strictly

```text
t>(r/mu) sqrt(b/(3h)).                                 (17)
```

On the **finite-exact 57-pair bank only**, the unique minimum is

```text
min h/b=51/19 at (p,q)=(1,10),
b=19/70, h=51/70.                                     (18)
```

It follows that all 57 scale-one types close whenever

```text
t mu/r>sqrt(19/153).                                  (19)
```

For scale two, let `ell(w)` count the pair-safe physical lifts of the two
solutions of `2x=y`.  Then

```text
meas(full safe)=(1/(2t)) integral N_t(w) ell(w) dw.    (20)
```

The useful quotient set `H_2={ell>0}` is the complement of `C_pq`.
Positivity in (20) is equivalent to `integral_(H_2)N_t>0`, although the two
quantities are not generally equal.  Therefore (14)--(17) apply as positivity
gates with `H_2`.  For the exact `(1,9)` quotient,

```text
b=4/63, h=59/63, h/b=59/4,
t mu/r>2/sqrt(177)                                    (21)
```

is sufficient.

## 6. Why pairwise tournament data cannot carry the target

Write `g=1_(G_u)`, `a=1_(D_(tp))`, `c=1_(D_(tq))` and

```text
M0=integral g,       Mp=integral g a,
Mq=integral g c,     Mpq=integral g a c.
```

The exact response is

```text
F(z,w)=M0+z Mp+w Mq+z w Mpq,
safe mass=F(-1,-1)=M0-Mp-Mq+Mpq.                      (22)
```

Pairwise Gram/cut data discard `Mpq`.  On AP11, the fixed survivor
`(p,q)=(6,41)` at `t=12` and `t=36` has identical complete one/two data

```text
Mp=11/1176,       Mq=1/123,
meas(D_(tp) intersect D_(tq))=5/246,                  (23)
```

but

```text
Mpq(t=12)=1/861,       Mpq(t=36)=89/61992.             (24)
```

The safe masses differ by `17/61992`.  An independent wall-cell engine
reproduces (23)--(24).  Hence an ordinary tournament, pairwise Gram matrix,
or cut semimetric cannot compile this target.  The minimal positive-mass
sidecar is the cubic response, equivalently the rooted integer carrier.
Endpoint equality still needs an atomic owner sidecar.

## 7. Exact controls and limitations

```text
AP11 scale-one distinct rows, 57 pairs, 11<=t<=84       4,214
positive safe mass                                      4,214
raw CV gate closures                                    4,113
integer-tariff closures                                 4,145
raw CV hostile scales                              t=13,26 only
all-pair raw CV closure                  t>=27 (exact to 87, BV from 88)

AP11 scale-two odd scales 11<=t<=99                         45
positive / exact-CV closures                            45/45

pairwise AP11 fibres in t=12..42                         1,739
fibres split by cubic response                              16
```

These are controls, not an arbitrary-body census.  The theorem gives no
uniform bound on `r/mu`, does not close the 17 types in (7) for an arbitrary
body, does not enter `t<U`, and does not prove LRC(14).

## 8. Reproduction

```text
python -B 04-computation/lrc14_auxiliary_center_erosion_response_thm3910.py
python -B 04-computation/lrc14_auxiliary_center_erosion_independent_audit_thm3910.py
python -B 04-computation/lrc14_safe_component_response_thm3910.py
python -B 04-computation/lrc14_safe_component_response_independent_audit_thm3910.py
python -B 04-computation/lrc14_t_sheet_lift_count_variance_thm3910.py
python -B 04-computation/lrc14_t_sheet_lift_count_ap11_boundary_thm3910.py
python -B 04-computation/lrc14_t_sheet_lift_count_scale2_thm3910.py
python -B 04-computation/lrc14_tournament_response_exact_scout_thm3910.py
python -B 04-computation/lrc14_tournament_response_collision_independent_audit_thm3910.py
```

Repeat with `python -B -O`; every stored transcript is byte-identical.
