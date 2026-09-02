# Half-turn minimum depth: an exact lower-dual and its cubic frontier

**Status:** DRAFT RESEARCH REFLECTION. The pointwise paired-minimum lemma,
the abstract marginal obstruction, and the displayed finite computations are
proved/exact as marked below. Universal positivity of the cubic functional is
**OPEN**. This file is not canon and does not claim LRC(14).

## 1. Inheritance pass

The closest proved scalar mechanism is
`THM-4335-lrc14-owner-permutation-component-address-and-minority-renewal.md`.
For a minority row

```text
S_h={2h} union W,             |W|=12,             every w in W odd,
```

it gives, with `E=G_{2h}` and physical safe mass `F`,

```text
F=A_0+Omega-6/7.                                      (1)
```

That identity is exact but forgets where the multiplicity lies.

The closest half-turn mechanism is
`THM-3385-odd-fibre-doubling-projection-and-half-even-complement-clocks.md`,
refined by
`THM-3387-exact-cyclic-sheet-cover-atlas-and-q2-gcd-graph.md`. Under the
degree-two quotient `y=2t`, an odd tail blocks at most one of the two
physical sheets. THM-3387 retains the Boolean full-cover locus and proves
the pairwise gcd graph, but not the depth of a full cover on the anchor-safe
core. `MISTAKE-382` is the endpoint-grid warning: Boolean open-cell words
must not silently be identified with the underlying point set.

The least-used sidecar here is therefore not another owner graph. It is the
two-sheet **minimum occupation depth**. This also distinguishes the present
construction from
`THM-4092-parity-weighted-antipodal-k-comb-density-compiler.md`, whose odd
antipodal comb merges the two sheets into one weight-two union.

The polynomial technology is inherited rather than new:

- `THM-406-covering-depth-master-object-factorial-moments-and-spectral-identity.md`
  identifies binomial moments with overlap volumes and records the
  Bonferroni direction;
- `THM-534-lrc-sector-moment-lp-dual-certificate.md` uses integer-state
  polynomial **majorants** of a zero-depth indicator;
- the present application uses a polynomial **minorant** of a paired free-
  sheet count on the new state `q=min(a,b)`.

`THM-638-signed-pair-mass-law-rational-thresholds.md` is the most immediate
candidate for controlling the next layer, because the cubic penalty begins
with simultaneous depth four and necessarily mixes the two sheet signs.

The canonical tight AP `{1,...,13}` is a scope guard, not an admissible
test: after anchoring at speed `2`, tail `4` is even and can hit both members
of a half-turn pair. Its usual witness `t=1/14` remains equality-safe.

## 2. The paired-minimum object

Put

```text
E_0=G_{2h} intersect [0,1/2).
```

For `t in E_0`, define

```text
a(t)=#{w in W: ||wt||<1/14},
b(t)=#{w in W: ||w(t+1/2)||<1/14},
q(t)=min(a(t),b(t)).                                  (2)
```

Because every `w` is odd,

```text
w(t+1/2)=wt+1/2 mod 1.
```

The two danger conditions for one label are disjoint. Consequently

```text
a+b<=12,                         0<=q<=6.             (3)
```

Let

```text
H_r=mu({t in E_0:q(t)=r}),                0<=r<=6.    (4)
```

The free members of the physical pair are counted by

```text
1_(a=0)+1_(b=0).
```

It follows exactly that

```text
F=H_0+B_0,     B_0=mu({t in E_0:a=b=0}).              (5)
```

In particular, `F>0` iff `H_0>0`. The extra `B_0` in `(5)` is important:
`H_0` counts pairs with at least one free member, whereas `F` counts free
physical points.

Under `y=2t`, condition `q>=1` is precisely THM-3387's two-sheet full-cover
predicate, restricted to `y in G_h`. Thus `(4)` refines the old Boolean
support by six multiplicity levels without changing the target predicate.

## 3. A finite exact lower-dual lemma

**PROVED.** Let

```text
p(q)=sum_(j=0)^d c_j binom(q,j)                       (6)
```

satisfy

```text
p(0)<=1,                         p(r)<=0 for 1<=r<=6. (7)
```

Then, pointwise on every half-turn pair,

```text
1_(a=0)+1_(b=0) >= p(min(a,b)),                       (8)
```

and hence

```text
F >= integral_(E_0) p(q(t))dt
  = sum_j c_j M_j,
M_j=sum_r binom(r,j)H_r.                              (9)
```

The proof is the seven-value check encoded in `(7)`: at `q=0` the left side
is at least one, and at `q>=1` it is zero. This is a seven-constraint rational
LP certificate format. Degree six recovers `H_0` exactly.

Three useful duals are

```text
p_2(q)=1-q+(1/3)binom(q,2)=(q-1)(q-6)/6,             (10)
p_3(q)=1-q+binom(q,2)-binom(q,3)
      =-(q-1)(q-2)(q-3)/6,                           (11)
p_6(q)=sum_(j=0)^6(-1)^j binom(q,j)=1_(q=0).         (12)
```

Their value vectors on `q=0,...,6` are

```text
p_2:  1,0,-2/3,-1,-1,-2/3,0,
p_3:  1,0,0,0,-1,-4,-10,
p_6:  1,0,0,0,0,0,0.                                (13)
```

Under the normalization `p_2(0)=1,p_2(1)=0`, the coefficient `1/3` in
`(10)` is sharp because the `q=6` constraint forces it. For the cubic,

```text
integral p_3(q)dt=H_0-H_4-4H_5-10H_6.               (14)
```

Therefore cubic positivity is exactly the new arithmetic target

```text
H_0>H_4+4H_5+10H_6.                                  (15)
```

Writing `T_j=mu(q>=j)`, the same condition is

```text
H_0>T_4+3T_5+6T_6.                                   (16)
```

## 4. First marginals do not imply the cubic target

This failure is exact, not heuristic and not stochastic independence.
Normalize the full base `X=[0,1/2)` to probability one. For each of twelve
labels, let `(L_i,R_i)` be its two disjoint sheet indicators. The exact odd-
comb marginals are

```text
E L_i=E R_i=1/7.                                      (17)
```

With probability `5/7`, make all labels absent. With probability `2/7`,
choose a uniformly labelled six-subset `A`, put the labels of `A` on the
left sheet and its complement on the right. Then every marginal in `(17)`
holds, but

```text
P(q=0)=5/7,              P(q=6)=2/7,
E p_3(q)=5/7-10(2/7)=-15/7.                           (18)
```

The value `-15/7` is sharp from these data. On the seven integer states,

```text
p_3(q)>=1-(11/6)q,                                   (19)
```

the chord through `q=0,6`. Also `Eq<=Ea=12/7`, so `(19)` gives the lower
bound `-15/7`, attained in `(18)`.

In fact the abstract admissible `q`-laws are exactly the probability laws
with `Eq<=12/7`. Necessity is immediate. For sufficiency, at state `q=r`
start with disjoint uniformly labelled `r/r` sets. With probability

```text
lambda=2(12/7-Eq)/(12-2Eq),                           (20)
```

put every unused label on a uniformly chosen side, and otherwise leave the
unused labels absent. This preserves `q=r`, and label and side randomization
give `(17)`. Here randomization is only a finite convex construction proving
an LP obstruction; it is not an assumption about runner events.

There is an equally sharp anchor-core abstract obstruction. Use the
subprobability measure `nu=2dt` on the core, of mass `6/7`; put `q=0` on
mass `4/7`, `q=6` on mass `2/7`, and make the removed anchor strip of mass
`1/7` empty. All full-base marginals remain `(17)`, while the core cubic is
`-16/7`.

Thus the missing coordinate is genuine shared-phase/interval geometry, not
another first-marginal calculation.

## 5. A current-sensitive deterministic survivor

Put `J=integral |a-b|`. Since

```text
q=(a+b-|a-b|)/2,                                     (21)
```

the chord `(19)` gives the pointwise inequality

```text
p_3(min(a,b))
 >=1-(11/12)(a+b-|a-b|).                             (22)
```

This is the strongest affine bound in `q`, with equality at `q=0,6`.

The normalization matters:

- On the full base `X` with probability-one normalization,
  `E(a+b)=24/7`. Thus

  ```text
  E p_3 >= -15/7+(11/12)E|a-b|,
  E|a-b|>180/77  ==>  E p_3>0.                        (23)
  ```

- On the anchor core, retain `nu=2dt`, so `nu(E_0)=6/7`. Deleting the
  anchor-danger strip only gives
  `L_E=integral_(E_0)(a+b)dnu<=24/7`. The row-adaptive condition is

  ```text
  J_E>L_E-72/77.                                      (24)
  ```

  Without another load sidecar, the uniform sufficient conditions are

  ```text
  J_E>192/77                  under nu=2dt,
  integral_(E_0)|a-b|dt>96/77 under physical dt,
  E(|a-b| | E_0)>32/11        under core probability. (25)
  ```

The often tempting `180/77` threshold belongs to the full base, not directly
to the anchor-restricted certificate. Transporting it across the deleted
anchor strip requires controlling that strip. The exact profiles below all
leave the chord bound negative, so `(22)` identifies a clean target but does
not close even the known controls.

## 6. Exact controls and hostiles

The companion performs a rational wall sweep with strict open danger teeth
and closed anchor components. It records all `H_r`, `F`, `B_0`, total load,
absolute current, and the three duals.

### 6.1 Quadratic universality is refuted

At `h=420`, take

```text
W_0=(9,841,1681,2521,3361,4201,5041,5881,
     6721,7561,8401,9241).                            (26)
```

This row is primitive, passes every denominator `2,...,14`, and the fixed
half-turn clocks are blocked by `5041,5881`. Exact evaluation gives

```text
integral p_2(q)dt
=-164696925438924148369133272304940284525104
 /4325895791028726433155124208050350034149135
~=-0.03807232846.                                    (27)
```

Thus `(10)` cannot be a universal positivity certificate. The first failed
implication is that a positive `H_0` need dominate the negative `q=2,...,5`
charges while the `q=6` constraint fixes the quadratic curvature. The cubic
survives this row:

```text
integral p_3(q)dt
=13755395733596017612242204554383223168195
 /288393052735248428877008280536690002276609
~=0.04769669589.                                     (28)
```

### 6.2 The cubic has not failed, but its margin can be small

No arithmetic row with nonpositive cubic was found. The following are exact,
not evidence of a theorem:

- Complete census of all `12`-subsets of the first `15` odd speeds and
  `1<=h<=30`: `13,650` rows, no nonpositive cubic. The minimum is

  ```text
  1446157/38057019
  ```

  at `h=2` and
  `(1,3,5,7,11,13,17,19,23,25,27,29)`.

- A modular near-zero scope hostile

  ```text
  h=4,
  W=(1,7,9,17,25,33,41,49,57,65,73,81)
  ```

  has `q_max=6` and

  ```text
  integral p_3(q)dt=373566525247/19204120390455
                    ~=0.01945241530.                  (29)
  ```

  It is outside the denominator-complete residual: denominators
  `6,10,12,14` are missing. A numerical complete scan of the `293,930`
  twelve-subsets of
  `{1} union {8j-1,8j+1:1<=j<=10}` selected exactly this row; the displayed
  value itself was then recomputed by the exact wall sweep.

- A large-height modular scope hostile

  ```text
  h=420,
  W=(1,839,841,1681,2521,3361,4201,5041,5881,
     6721,7561,8401)
  ```

  has `q_max=6`, misses only denominator `9`, and gives

  ```text
  integral p_3(q)dt
  =28959147561368969682651461711510842816
   /844629198923707571893175739849000114881
  ~=0.03428622595.                                   (30)
  ```

  Perturbing and optimizing the affine families at large height did not
  produce a negative row; that exploratory statement is not exhaustive.

### 6.3 Inherited `h=420` controls

For THM-4335's two joint controls, the quadratic and cubic values are both
positive:

```text
P=1287:
 p_2=551054571711629821461569179
     /23938321074154037836892907300 ~=0.02301976693,
 p_3=14783540775565989336572164
     /153450776116372037415980175   ~=0.09634060609;

P=9009:
 p_2=37049456465357700644747
     /1579774373005611947264100     ~=0.02345237212,
 p_3=16362496058066898986258
     /169261539964886994349725      ~=0.09666989950. (31)
```

Both have `q_max=4`. Their exact `F` agrees with the occupation computation
in THM-4335, and `(5)` is audited independently as `F=H_0+B_0`.

The additional height-sharp THM-366-complete control

```text
h=420,
W=(1,735,1155,1365,1679,1681,3255,3359,3361,
   3609,5039,5041)                                   (32)
```

also gives positive values

```text
p_2=1877339701087410921940449985333
    /133247162596101650207195390175870 ~=0.01408915330,
p_3=17339729627705789770092924238
    /180551710834826084291592669615    ~=0.09603747064. (33)
```

Here `q_max=5`. This is informative because the row blocks the old grid/unit
banks yet remains comfortably cubic-positive; the cubic is not merely
restating those certificates.

## 7. Source-to-target connection contract

| Field | Exact content |
|---|---|
| source | THM-4335's scalar occupation identity and THM-3387's degree-two sheet deck |
| target | the core-relative measure vector `(H_0,...,H_6)` and its finite lower-dual cone |
| map | pair `t,t+1/2`, or quotient `y=2t`; labelled sheet counts `(a,b)` map to `q=min(a,b)` |
| preserved | strict danger status, half-turn location, anchor-core membership, minimum cover depth, Haar mass, and the predicate `F>0 iff H_0>0` |
| lost by `q` | tail labels and tooth addresses, imbalance sign, the separate counts `(a,b)`, endpoint order, renewal transitions, and cross-component reuse |
| lost by `(A_0,Omega)` | all half-turn pairing and depth placement |
| needed sidecars | signed current `a-b`; labelled mixed-sign intersections; THM-4335 addresses `(w,n,q_i)`; component index and deleted anchor-strip load |

This answers the novelty audit narrowly: polynomial moment duals and the
two-sheet Boolean quotient already existed. The new object is their
core-relative **minimum-depth refinement**, together with the lower-dual
direction and exact occupation link `(5)`.

## 8. Cheapest continuation

The cubic lane is worth one focused continuation, not a claim of closure.
The cheapest decisive targets are:

1. Prove or refute `(15)` on the denominator-complete `420|h` residual.
   Search affine residue families and gate-complete perturbations, but
   exactify only candidates near zero.
2. Bound the balanced-depth tail

   ```text
   H_4+4H_5+10H_6
   ```

   by `H_0`. THM-3387's gcd graph controls the Boolean support but needs
   labelled four-fold extensions; THM-638 supplies exact signed pair masses
   for the first layer.
3. Use the row-adaptive current criterion `(24)`, not the coarser uniform
   threshold. The missing quantity is the joint core load/current pair
   `(L_E,J_E)` plus the anchor-strip sidecar.
4. If cubic positivity fails, solve the seven-state rational LP after adding
   whatever exact arithmetic upper bounds are available for `H_4,H_5,H_6`.
   Degree six always returns `H_0` exactly, so this is a finite certificate
   format even when low degree fails.

No stochastic independence enters any statement above.

## Reproduction

```bash
python3 -B 04-computation/lrc14_halfturn_min_depth_dual_probe_occupation_dual_20260902.py \
  | diff -u 05-knowledge/results/lrc14_halfturn_min_depth_dual_probe_occupation_dual_20260902.out -
```
