# LRC14 disconnected-low affine limit and tail certificate

**Status.** PROVED analytic reduction + exact rational limit atlas; the
finite residual is certified by the companion compiler and promoted as
THM-3360. This note concerns the disconnected-low physical pair floor only.
Together with the connected-low theorem it independently recovers the current
matched-residue reflected closure, but it does not close projected, rung,
physical-entry, or arbitrary-`k<=1` branches of LRC(14).

## 1. Inherited boundary

The closest proved mechanism is the arbitrary-ratio Dirichlet reduction in
`lrc14-generalized-dirichlet-resonance-reduction-20260812.md`. After the
large-ruler, far-ratio, `3:5`, and `g>=4` theorems and the later head-263
certificate, its sole residual is:

* one of the 2,530 body-safe small-ruler contexts `(L,j,e,f)`;
* a nonzero-resonance affine ray

  `p=p0+dn, q=q0+(d+a)n, dq-(d+a)p=c`,

  with `1<=d<=8`, `0<=a<=7d`, `1<=|c|<=46`, and `1<=p0<=d`;
* raw `p>=264`, `p<q<8p`, and `gcd(p,q)<=3`;
* the inherited Dirichlet admissibility condition `9|c|<=p`.

There are 22,890 canonical ray rows. This is a cover, not a disjoint
partition. MISTAKE-376 forbids routing that finite bank as closed without an
actual raywise certificate.

## 2. Exact homogenized limit

Let `D=d+a`, `k=gcd(d,D)`, `P=d/k`, and `Q=D/k`. Put

```text
alpha=p0-e/L,       beta=q0-f/L,
u=(ej mod L)/L,     v=(fj mod L)/L,
A=Lc+De-df.
```

Body safety means that the physical overlap is literally

```text
I_n=int_0^1 chi((knP+alpha)x-u) chi((knQ+beta)x-v) dx,       (1)
```

where `chi` is the circular radius-`1/14` indicator. There are no clipped or
extra endpoint teeth.

For the periodized tent `T_R(theta)=sum_m [R-|theta+m|]_+`, define

```text
Phi_(P,Q)(theta)
 = [T_((P+Q)/14)(theta)-T_(|P-Q|/14)(theta)]/(PQ).          (2)
```

The exact ray limit is

```text
J = int_0^1 Phi_(P,Q)(theta0-A t/(kL)) dt,
theta0=(-D(ej mod L)+d(fj mod L))/(kL).                     (3)
```

Equation (3) is rational: a primitive of each tent is piecewise quadratic
with rational breakpoints. It depends on `(d,a,c,L,j,e,f)` but not on the
residue representative `(p0,q0)`.

The exact atlas has 17,206 distinct `(d,a,c)` triples and 79 context lanes
not already handled by the monotone midpoint envelope, hence 1,359,274
rational limit evaluations. Its minimum is

```text
J_min = 709/48048
```

at `(d,a,c;L,j,e,f)=(3,8,-1;168,90,12,1)`, with strict gap

```text
J_min-1/294 = 1273/112112 > 0.                              (4)
```

Thus no affine ray approaches the target wall. What remained was an
effective, uniform convergence estimate.

## 3. The blockwise averaging lemma

Set

```text
N=kn,              epsilon=alpha/P,
B=P beta-Q alpha=A/(kL),
delta=B/[P(N+epsilon)].
```

The change of variable `y=(N+epsilon)x` turns (1) into

```text
I_n = 1/(N+epsilon) int_0^(N+epsilon)
      chi(Py-u) chi(Qy+delta y-v) dy.                      (5)
```

Let

```text
h(s)=int_0^1 chi(Py-u) chi(Qy+s-v)dy.
```

Translation of one circular interval gives the elementary sharp-enough
Lipschitz estimate

```text
|h(s)-h(t)| <= 2|s-t|.                                    (6)
```

The same constant controls the within-block frequency change, and this point
needs proof because the shift is `delta y`, not constant. Interpolate

```text
b_lambda(y)=chi((Q+lambda delta)y+s),   0<=lambda<=1.
```

Both endpoint frequencies are at least one: `Q>=1`, while

```text
Q+delta=P(Lq-f)/(Lp-e)>1.
```

For a frequency `R>=1`, one boundary family has roots
`y=(m+theta)/R` in `[0,1]`. If `R=M+r`, `M>=1`, summing the arithmetic
progression separately for `theta<=r` and `theta>r` gives

```text
sum_(0<=m+theta<=R) (m+theta) <= R^2.
```

Distributional differentiation in `lambda` therefore costs at most
`|delta|` per boundary family: the root Jacobian contributes `1/R`, the
motion contributes `|delta|y`, and their sum is
`|delta| sum(m+theta)/R^2`. There are two boundaries, so each full block
differs from its frozen-shift version by at most `2|delta|`. Endpoint roots
follow by one-sided approximation and do not add mass.

Split (5) into its first `N` unit blocks and a final segment of length
`epsilon`.

1. The full-block average is at most `1/7`. The final segment is bounded by
   the first clock, whose mean is `1/7` and whose mean-zero primitive has
   oscillation `6/(49P)`. Normalization plus that segment therefore costs at
   most `[2epsilon/7+6/(49P)]/(N+epsilon)`.
2. Replacing `delta(j+y)` by `delta j` on each full block costs at most
   `2|delta|` after normalization.
3. Integrating (6) inside each Riemann cell shows that the left sum for
   `h(delta N t)` costs at most `|delta|`.
4. Replacing its endpoint `delta N` by the limiting endpoint `B/P` costs at
   most `|B|epsilon/[P(N+epsilon)]`.

Therefore

```text
|I_n-J| <=
 [2epsilon/7+6/(49P)+|B|(3+epsilon)/P]/(N+epsilon).       (7)
```

This estimate keeps the finite first-clock offset; dropping it would be an
invalid asymptotic substitution. If `z0=Lp0-e` and `z=Lp-e`, (7) becomes

```text
|I_n-J| <=
 [2z0/7 + 6L/49 + 3|A|/k + |A|z0/(dL)] / z.              (8)
```

Every quantity in (8) is rational and affine in the ray parameter. Combining
(8) with the exact positive gap `J-1/294` gives an explicit finite head on
every ray/context row.

## 4. Orthogonal exits and exact residue

The compiler uses four proof-preserving reductions.

* The monotone midpoint envelope leaves only 79 context lanes, with active
  counts `(79,10,4)` for raw gcd `g=1,2,3`.
* Since `gcd(p,q)` divides `c`, it is constant on each class
  `n mod |c|`. Classes outside `g<=3` are removed without iteration.
* The many-turn theorem fires whenever
  `floor(p/d)|A|/z>=5` **and `p>=273`**; these are exact affine exit
  conditions.  Its weaker `p>=264` conclusion beats `D_max/5` and would
  already suffice after summing a five-edge tree, but it does not by itself
  prove the advertised edgewise `1/294` floor.  The compiler therefore sends
  the entire strip `264<=p<273` through the exact residual path.
* Before literal evaluation, the contextwise THM-3350 midpoint bound is
  compared exactly with `1/294`.

The compiler starts each formal ray only when `9|c|<=p`; points before that
threshold are not incidences supplied by the Dirichlet reduction.  The
`p>=273` scope repair also sends the entire admissible strip
`264<=p<273` through exact certification.  This is why the conclusion below
is an edgewise `1/294` floor on every Dirichlet-admissible incidence rather
than only the weaker five-edge-average closure available from the many-turn
theorem in that strip.

Every remaining row is evaluated by an integer `__int128` port of the
canonical Euclidean floor-moment mass engine. Deterministic controls compare
the port against both the canonical Python engine and an independent exact
implementation of (3). The ordinary and optimized wrapper runs reproduce
the same semantic summary.

## 5. Consequence

Together with the head-263 certificate, (3)--(8) discharge the 22,890-ray
residual in the generalized Dirichlet reduction. Hence every high
cross-component edge in every disconnected-low physical context has overlap
at least `1/294`. The explicit five-edge complete-multipartite tree then has
credit at least `5/294`, exceeding the exact singleton debt by

```text
570672686921/494472037809750 > 0.
```

Thus the disconnected-low branch closes. This is a substantial LRC(14)
subbranch theorem, not a proof of LRC(14).

## 6. Relation to the concurrent horn-tree proof

THM-3355 reached the same branch closure independently. Its centered-grid
lower bound proves the edgewise `1/294` floor outside four oriented
scale-one lanes `(L,j,e,f)=(168,90,12,f)`, `f in {1,2,3,4}`, and then pays
those analytic horns by a debt-sensitive weighted tree. The affine atlas
retains all four lanes. Its exact limit minimum occurs at the first of them,
yet is still `1273/112112` above the wall, and the complete admissible finite
heads have no failure. Therefore the horns are defects of that particular
centered-grid estimate, not exceptions to the physical floor.

The two proofs preserve different information. THM-3355 keeps the six-level
component graph and can trade edge quality against the actual singleton debt;
THM-3360 forgets the surrounding assignment but proves every individual high
pair safe. The former is the shorter branch closure, while the latter is the
reusable pairwise strengthening.

## Reproduction

```text
python 04-computation/lrc14_disconnected_low_affine_tail_kps_s171.py --threads 12
python -O 04-computation/lrc14_disconnected_low_affine_tail_kps_s171.py --threads 12
```

Companion artifacts:

* `04-computation/lrc14_disconnected_low_affine_tail_kps_s171.cpp`
* `04-computation/lrc14_disconnected_low_affine_tail_kps_s171.py`
* `05-knowledge/results/lrc14_disconnected_low_affine_tail_kps_s171.out`
