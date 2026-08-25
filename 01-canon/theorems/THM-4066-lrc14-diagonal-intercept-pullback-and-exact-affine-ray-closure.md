---
id: THM-4066
title: "LRC(14) diagonal-intercept pullback and exact affine-ray closure"
status: >
  PROVED RELATIVE TO CITED LRCUpTo13 + VERIFIED-EXACT + INDEPENDENTLY
  AUDITED. Common unit dilation of every exception conjugates its labelled
  lift masks through one diagonal affine-intercept coordinate. The primitive
  d=3 shape (1,4,5) and d=4 shape (2,7,9) have exact two-arc spoiled-phase
  windows. Consequently 3H union h(1,4,5) closes when 24h>=11 max(H), and
  4H union h(2,7,9) closes when 21h>=22 max(H), strictly improving the
  THM-4052 cones on these certificate-positive rays. This closes two exact
  affine families, not LRC(14).
source: codex-frontier-synthesis-creative-20260825b / LRC affine-intercept lane
audit: >
  PASS. The primary engine performs literal rational wall-cell subdivision,
  checks 63,200 labelled phase conjugacies, every strict endpoint, the
  nonunit hostile, runner typing, and explicit gains. The independent path
  starts from the THM-4030/4032 affine-centre boxes, projects their real
  interval intersections to the circle, and verifies conjugacy algebraically.
  Normal and optimized outputs are byte-identical; both scripts have zero
  assert nodes and zero float literals.
depends_on:
  - LRCUpTo13
  - THM-4030-lrc14-d4-affine-defect-lattice-boundary
  - THM-4032-lrc14-d3-affine-defect-lattice-boundary
  - THM-4052-lrc14-affine-component-width-escape-cones
  - THM-4062-lrc-divisor-star-affine-intercept-obstruction
related:
  - THM-4024-lrc14-complete-divisor-incidence-envelope
  - THM-4041-lrc14-d2-affine-defect-edge-boundary
  - THM-859-hamming-six-dilation-conjugacy-and-order-one-gate
script: 04-computation/lrc14_diagonal_intercept_affine_ray_thm4066.py
output: 05-knowledge/results/lrc14_diagonal_intercept_affine_ray_thm4066.out
independent_audit_script: 04-computation/lrc14_diagonal_intercept_affine_ray_thm4066_independent_audit.py
independent_audit_output: 05-knowledge/results/lrc14_diagonal_intercept_affine_ray_thm4066_independent_audit.out
script_sha256: 0d38a5d320b47f1815202e3732e889e9b377ae485b748d25ee965cbe8d923306
output_sha256: 5445d332c437c4b0178fcfac1cc292ad9a098d4e7fb9a2ed5b5bca8c4d6380e7
independent_audit_script_sha256: fcb1a36a062c3869d871095b64a1d7817230589736697b579bc88eb822f92f3d
independent_audit_output_sha256: 0f223fc284b21ee697ba1a6062511ac6957c9afa69d152a23d5e4a3e27ebb54d
hash_basis: raw LF bytes
---

# THM-4066 -- diagonal-intercept pullback and exact affine-ray closure

**PROVED RELATIVE TO CITED LRCUpTo13 + VERIFIED-EXACT + INDEPENDENTLY
AUDITED.** THM-4062 identifies the missing selected-lift coordinate as an
affine intercept, while THM-4030/4032 characterize the remaining `d=4,3`
spoiled phases one exception triple at a time. This theorem finds a lawful
one-dimensional orbit inside that intercept field and uses its exact
component widths to close two infinite equality-shape families. LRC(14)
remains open.

Write `||x||` for distance to the nearest integer.

## 1. Inheritance and the diagonal intercept

For a divisor `d`, an exception `delta`, and a chosen real representative
`y` of a divided-pack phase, retain THM-4062's labelled mask

```text
D_delta^(d)(y)={j in C_d: ||delta(y+j)/d||<1/14}.      (1)
```

For a finite exception set `E`, put

```text
Sigma_d(E)={y in R/Z:
              union_(delta in E)D_delta^(d)(y)=C_d}.  (2)
```

Thus `Sigma_d(E)` is the strict fully-spoiled phase set. Its representative
gauge is harmless only because `(2)` forgets a simultaneous label rotation;
the individual masks in `(1)` retain that rotation.

The closest proved mechanism is THM-4062's `(q,a,tau)` factorization. The
canonical hostile is its `d=4` same-static-packet pair `(2,9,11)` and
`(2,1,3)`. The corrected near miss is static divisor/depth completion. The
least-used sidecar is the full intercept orbit as `y` moves. On a common
exception ray `hE_0`, all intercepts become functions of the one coordinate

```text
s=hy in R/Z.                                           (3)
```

This remains a dynamic coordinate; it is not another static packet.

## 2. Exact diagonal-intercept pullback

Let `h>=1` and `gcd(h,d)=1`. For `S subset C_d`, write

```text
hS={hj mod d:j in S}.                                  (4)
```

Then for every exception `delta` and real representative `y`,

```text
h D_(h delta)^(d)(y)=D_delta^(d)(hy).                  (5)
```

Indeed, for `k=hj mod d`, write `hj=k+nd`. Then

```text
h delta(y+j)/d=delta(hy+k)/d+n delta,                 (6)
```

and the last term is integral. Multiplication by `h` permutes `C_d`, so
`(5)` gives the exact set identity

```text
Sigma_d(hE)=[h]^(-1)Sigma_d(E),                       (7)
```

where `[h](y)=hy` is the degree-`h` circle cover. Reducing the real number
`hy` modulo one rotates the right-hand labels by the integer part of `hy`,
exactly as THM-4062's covariance requires; coverage in `(7)` is unchanged.

Haar measure is therefore preserved. If `Sigma_d(E)` has proper open
components of lengths `L_1,...,L_c`, then `Sigma_d(hE)` has exactly `hc`
components: each `L_i/h` occurs `h` times. This is a labelled-lift analogue
of THM-859's safe-set dilation conjugacy, with the affine phase gauge kept.

The condition `gcd(h,d)=1` is load-bearing. Section 7 gives a hostile when it
fails.

## 3. The primitive `d=3` window

For the primitive exception shape `(1,4,5)`, direct strict-mask arithmetic
gives

```text
Sigma_3({1,4,5})
 =(11/56,3/14) union (11/14,45/56).                   (8)
```

Here is a proof with no phase scan. On `0<=y<1`, speed one can kill only
label zero when `y<3/14`, or label two when `y>11/14`. In the first range,
the only possible complementary assignment is speed five on label one and
speed four on label two. The three strict windows are

```text
y<3/14,
11/70<y<17/70,
11/56<y<17/56.                                       (9)
```

Their intersection is `(11/56,3/14)`. The reverse assignment is empty.
Reflection gives `(11/14,45/56)`. On the two components, in speed order
`(1,4,5)`, the mask words are respectively

```text
({0},{2},{1}),             ({2},{0},{1}).             (10)
```

Every endpoint in `(8)` is excluded. The literal strict endpoint masks are

```text
y=11/56: ({0},{},{1}),       y=3/14:  ({},{2},{1}),
y=11/14: ({},{0},{1}),       y=45/56: ({2},{},{1}).   (11)
```

Each leaks the displayed complementary label; changing `<1/14` to
`<=1/14` covers all three labels. Thus openness is exact, not notation.

By `(7)`, whenever `3` does not divide `h`,

```text
y in Sigma_3(h{1,4,5})
 iff hy mod 1 lies in the two arcs in (8).             (12)
```

The scaled spoiled set has `2h` components, each of length `1/(56h)`, and
total measure `1/28`.

## 4. The primitive `d=4` window

For the primitive typed shape `(2,7,9)`,

```text
Sigma_4({2,7,9})
 =(5/49,1/7) union (6/7,44/49).                       (13)
```

Speed two kills the parity pair `{0,2}` for `y<1/7`, the pair `{1,3}` for
`y>6/7`, and no full parity pair in between. In the first range, the only
possible complementary assignment is speed seven on label one and speed
nine on label three. Their windows give

```text
(5/49,9/49) intersect (5/63,1/7) intersect [0,1/7)
 =(5/49,1/7).                                         (14)
```

The reverse assignment is empty, and reflection gives the second arc.
The mask words, in speed order `(2,7,9)`, are

```text
({0,2},{1},{3}),           ({1,3},{2},{0}).            (15)
```

Again every endpoint is genuinely strict:

```text
y=5/49: ({0,2},{},{3}),      y=1/7:  ({},{1},{}),
y=6/7:  ({},{2},{}),         y=44/49:({1,3},{},{0}).  (16)
```

The closed masks cover all four labels. For every odd `h`, `(7)` now gives

```text
y in Sigma_4(h{2,7,9})
 iff hy mod 1 lies in the two arcs in (13).            (17)
```

There are `2h` components, each of length `2/(49h)`, and the total measure
is `4/49`.

The independent audit obtains `(8)` and `(13)` through a genuinely different
route. It enumerates the exact affine-centre boxes of THM-4032/4030, takes
their real open-interval intersections, and projects them to the phase
circle. For `d=3`, two real intervals project to the two arcs in `(8)`; for
`d=4`, four real lift intervals project pairwise to the two arcs in `(13)`.

## 5. Exact affine-ray closure

Let `H` be any ten-element set of distinct positive integers and put

```text
M=max H.                                               (18)
```

The following two LRC(14) families close.

1. If `3` does not divide `h` and

   ```text
   24h>=11M,                                           (19)
   ```

   then

   ```text
   V_3=3H union {h,4h,5h}                             (20)
   ```

   is `1/14`-lonely.

2. If `h` is odd and

   ```text
   21h>=22M,                                           (21)
   ```

   then

   ```text
   V_4=4H union {2h,7h,9h}                            (22)
   ```

   is `1/14`-lonely.

The runner count and disjointness are exact. In `(20)`, none of
`h,4h,5h` is divisible by three. In `(22)`, `2h` is two modulo four while
`7h,9h` are odd. Hence each exception triple is disjoint from `dH`, its
three entries are distinct, and `(20)/(22)` contain exactly thirteen
nonzero speeds. With the stationary runner, these are LRC(14) rows. The
ten nonzero speeds in `H` receive the cited lower-dimensional clearance
`1/11`.

Choose a cited pack phase `y_0` with

```text
||r y_0||>=1/11              for every r in H.         (23)
```

The one-Lipschitz argument of THM-4052 puts the closed circle arc

```text
I={y:dist(y,y_0)<=3/(154M)}                           (24)
```

inside the `1/14`-safe set `G(H)`. Its length is `3/(77M)`. Every lift
`x_j=(y+j)/d` preserves all pack speeds `dH`. Therefore a non-lonely full
row would force

```text
I subset Sigma_d(hE_0).                               (25)
```

For `d=3`, condition `(19)` is exactly

```text
3/(77M)>=1/(56h).                                    (26)
```

For `d=4`, condition `(21)` is exactly

```text
3/(77M)>=2/(49h).                                    (27)
```

A connected closed arc cannot lie inside a strictly shorter open component,
or inside an open component of equal length. Since a connected subset of
the disjoint open set in `(25)` lies in one component, `(26)/(27)` contradict
`(25)`. This proves both closures, including equality.

## 6. Strict gains and physical controls

For the `d=3` ray, THM-4052 uses the largest-tooth bound

```text
3/(7E)=3/(35h),
```

whereas the exact width is `1/(56h)`. The gain factor is

```text
(3/35)/(1/56)=24/5.                                  (28)
```

Equivalently, the new threshold is `h/M>=11/24`; the old coarse cone needs
`h/M>=11/5`.

For `d=4`, `4/(7E)=4/(63h)` becomes `2/(49h)`, a gain of

```text
(4/63)/(2/49)=14/9.                                  (29)
```

The new threshold is `h/M>=22/21`, versus `44/27` in the old cone.

Both gains occur on literal typed rows below the old cones. Take
`H={1,...,10}`.

For `d=3,h=5`,

```text
V_3=(3,6,9,12,15,18,21,24,27,30,5,20,25),
24h=120>=110,             5h=25<110.                 (30)
```

It has the physical `11+2` decomposition

```text
body=(6,9,12,15,18,21,24,27,5,20,25),
pair=(3,30),                                             (31)
```

with primitive body, `c_3=8`, and all speeds distinct. The phase `x=1/33`
has clearance `1/11`.

For `d=4,h=11`,

```text
V_4=(4,8,12,16,20,24,28,32,36,40,22,77,99),
21h=231>=220,             27h=297<440.               (32)
```

Its typed decomposition is

```text
body=(8,12,16,20,24,28,32,36,22,77,99),
pair=(4,40),                                             (33)
```

with primitive body and `c_4=8`; `x=1/44` has clearance `1/11`.

The primitive shapes `(1,4,5)` and `(2,7,9)` are certificate-positive in
THM-4032 and THM-4030, respectively. Thus these are not certificate-negative
closures in disguise. The new operation is exact component width along the
diagonal intercept orbit.

## 7. Controls, general ray functional, and loss ledger

The `d=2` shape `(3,5)` is a consistency control. THM-4041 gives

```text
Sigma_2({3,5})
 =(13/35,8/21) union (13/21,22/35).                   (34)
```

Its common odd dilation has `2h` components of width `1/(105h)` and measure
`2/105`. Comparing with the eleven-pack safe arc `1/(42M)` gives `5h>=2M`,
exactly the `(3,5)` specialization already present in THM-4052. The genuinely
new gains are at `d=3,4`.

The unit condition in `(7)` cannot be dropped. At `d=3,h=3`, the set
`(3,12,15)` no longer has the equality type: every exception is divisible by
three, so each mask is either all of `C_3` or empty. Its fully-spoiled set is
the union of the ordinary danger sets for speeds `(1,4,5)`, with exact
measure

```text
51/140,                                                (35)
```

whereas a false degree-three pullback of `(8)` would retain measure `1/28`.
The independent interval union and the literal wall engine agree on `(35)`.

More generally, for a primitive equality-type exception shape `E_0`, define

```text
W_d(E_0)=maximum component length of Sigma_d(E_0),    (36)
```

with `W_d=0` when the spoiled set is empty. THM-4032/4030 make `(36)` a
finite exact rational computation from the affine-centre boxes. The same
proof yields the reusable ray certificate

```text
gcd(h,d)=1 and W_d(E_0)/h<=3/(77M)
 ==> dH union hE_0 is 1/14-lonely,          d in {3,4}. (37)
```

Thus, for fixed `H`, every fixed primitive exception ray has a finite initial
residual; all later common dilations close. THM-4052 replaces `(36)` by a
coarser largest-tooth bound.

The connection ledger is:

| entry | exact content |
|---|---|
| source | labelled lifts over `y`, with exception ray `hE_0` |
| target | diagonal intercept `s=hy` and primitive cover window `Sigma_d(E_0)` |
| map | `y->hy`, `j->hj mod d` |
| preserved | strict full spoilage, Haar measure, component spectrum up to `1/h` |
| destroyed | inverse sheet in `C_h`, absolute placement against `G(H)`, owner labels if the permutation is dropped |
| sidecar | inverse component/sheet address, `h mod d`, and the pack-safe set or its closed arc |
| hostile | nonunit `d=h=3`, where measure changes from `1/28` to `51/140` |

## 8. Scope and replay

This theorem closes two common-dilation affine shape families and supplies
the general exact ray certificate `(37)`. It does not classify arbitrary
primitive exception shapes, decide arbitrary intersections with `G(H)`,
close the remaining finite THM-3818 producer, handle other `11+2` types, or
prove LRC(14).

Reproduce both exact paths from the repository root:

```text
python3 -B 04-computation/lrc14_diagonal_intercept_affine_ray_thm4066.py
python3 -B -O 04-computation/lrc14_diagonal_intercept_affine_ray_thm4066.py
python3 -B 04-computation/lrc14_diagonal_intercept_affine_ray_thm4066_independent_audit.py
python3 -B -O 04-computation/lrc14_diagonal_intercept_affine_ray_thm4066_independent_audit.py
```

Both normal/optimized pairs reproduce their frozen raw-LF outputs. **QED.**
