# Half-turn cubic depth: nonlinear current, Brownian energy, and an AP-step hostile

**Status:** FINAL RESEARCH REFLECTION.  The pointwise identities, the
Brownian covariance specialization, and the displayed rational wall sweeps
are exact.  Universal positivity of the `q=min(a,b)` cubic remains **OPEN**.
The bivariate cubic considered below is a valid pointwise lower dual, but its
uniform positivity is **REFUTED**.  This file does not claim LRC(14).

## 1. Inheritance and board

The source object is the paired occupation profile from
`lrc14-halfturn-min-depth-dual-occupation-dual-20260902.md`.  For twelve odd
tails over an anchor-safe half-base, write

```text
a(t)=low-sheet danger depth,
b(t)=high-sheet danger depth,
q(t)=min(a(t),b(t)),             a+b<=12.
```

The inherited open target is

```text
Q_3 := integral [1-q+C(q,2)-C(q,3)]
     = H_0-H_4-4H_5-10H_6 > 0.                       (1)
```

The canonical hostile is the modular `h=420` row with six-fold depth from
the previous reflection.  The corrected near miss is the quadratic dual,
which is already negative on exact denominator-complete rows.  The
least-used sidecar is the signed sheet current `C=a-b`, not another scalar
load moment.

The live concepts in this pull were:

1. the minimum-depth law `(H_0,...,H_6)`;
2. one-sheet Bonferroni depth;
3. signed current and its gcd-reduced mod-14 covariance;
4. the anchor-strip loss between full and core energy;
5. affine/AP residue families as hostile generators.

## 2. The minimum cubic is a nonlinear current certificate

Put

```text
r(x)=1-x+C(x,2)-C(x,3)=-(x-1)(x-2)(x-3)/6.           (2)
```

On the integer states `0,...,12`, `r` is nonincreasing and has the exact
plateau

```text
r(1)=r(2)=r(3)=0.                                    (3)
```

Consequently the nonpolynomial minimum operation has the exact form

```text
r(min(a,b))=max(r(a),r(b)).                           (4)
```

Let `E` be a sheet-symmetric half-base region (in particular the anchor core
`G_(2h) intersect [0,1/2)`).  Reflection swaps `a` and `b`, so if

```text
B_3(E)=integral_E r(a)=integral_E r(b),
V_E=integral_E (a-b)^2,
```

then `(4)` gives the exact identity

```text
integral_E r(min(a,b))
 =B_3(E)+(1/2) integral_E |r(a)-r(b)|.                (5)
```

Writing `S=a+b` and `C=a-b`, direct factorization gives

```text
|r(a)-r(b)|
 =|C|[C^2+3(S-4)^2-4]/24.                            (6)
```

The bracket is nonnegative on every unequal admissible integer state.  It
vanishes precisely between different points of the plateau `{1,2,3}`.
Thus the old `q`-cubic is an ordinary one-sheet third Bonferroni truncation
plus an exact **nonlinear cubic-current rebate**.  This is why it can retain
signal invisible to polynomial moments in `(a,b)` of the same nominal degree.

The exact state audit also proves

```text
|r(a)-r(b)| >= max(C^2-4,0)/6.                        (7)
```

Equality holds when `a=b`, when both depths lie in `{1,2,3}`, and at
`(a,b)=(0,4),(4,0)`.  Hence

```text
Q_3(E) >= B_3(E)+(1/12) integral_E max(C^2-4,0)
        >= B_3(E)+(V_E-4|E|)/12.                      (8)
```

For the anchor core `|E|=3/7`, the last term is

```text
B_3(E)+(V_E-12/7)/12.                                (9)
```

This is a valid smaller arithmetic certificate, but not a universal closure:
the exact AP-step rows below keep `(9)` negative even though `(1)` is
positive.  The positive-part version in the first line of `(8)` is sharper;
on the denominator-complete modular hostile below it is positive, whereas
the plain variance version is negative.

## 3. THM-638 becomes a three-step Brownian current kernel

For an odd speed `w`, put

```text
sigma_w(t)
 =1_(||wt||<1/14)-1_(||w(t+1/2)||<1/14).             (10)
```

Then `C=sum_w sigma_w`.  For two odd speeds let

```text
g=gcd(u,v),             U=u/g,             V=v/g,
d_14(n)=distance from n to 14Z.
```

The exact signed-pair specialization is

```text
K(u,v):=integral_0^1 sigma_u sigma_v
       =[d_14(U+V)-d_14(U-V)]/(7UV).                  (11)
```

For the odd reduced residues, define

```text
residue        1  5  3 | 13  9 11 | 7
epsilon        +  +  + |  -  -  - | 0
ell            1  2  3 |  1  2  3 | 0.
```

The numerator in `(11)` factors exactly as

```text
M(r,s)=2 epsilon(r)epsilon(s) min(ell(r),ell(s)).     (12)
```

Thus the residue numerator is the rank-three Gram matrix of nested signed
features, or twice the covariance of Brownian motion at times `1,2,3`, with
opposite residues reversing orientation and residue seven null.  This is a
literal factorization, not a stochastic-independence assumption.

Because `C^2` is half-turn invariant,

```text
V_[0,1/2)=(1/2) sum_(u,v in W) K(u,v).                (13)
```

But the target lives on the anchor core.  If

```text
A_h=D_(2h) intersect [0,1/2),
```

then exactly

```text
V_E=V_[0,1/2)-integral_(A_h) C^2.                    (14)
```

The Brownian kernel computes the first term and says nothing by itself about
the second.  The **anchor-strip current energy** is therefore the precise
position sidecar lost by full-circle pair covariance.  Even after it is
retained, second-order energy alone is diagnostic rather than sufficient:
the AP-step controls have negative variance certificate `(9)`, while their
exact nonlinear rebate `(6)` makes `(1)` positive.

## 4. A symmetric bivariate cubic and its exact failure

The proposed mixed cubic is

```text
P(a,b)=4/5-(a+b)/5
       +(7/15)[C(a,2)+C(b,2)]-(2/5)ab
       -(3/5)[C(a,3)+C(b,3)]
       +(2/15)[C(a,2)b+aC(b,2)].                      (15)
```

It is a valid pointwise lower dual:

```text
P(a,b)<=1_(a=0)+1_(b=0),             a+b<=12.        (16)
```

There is a short finite-lattice proof.  With `n=a+b,d=a-b`,

```text
P(a,b)=-[(n-8)(n-6)(n-2)+(11n-48)d^2]/120.           (17)
```

For `a,b>=1`, check `n=2,3,4` directly; at `n=5` the two terms in the
bracket are positive; at `n=6,8` the first term vanishes and the second is
nonnegative; at `n=7` parity gives `d^2>=1`; and for `n>=9` both terms are
positive.  On the boundary `b=0`,

```text
1-P(a,0)=(a-3)(3a^2-7a-2)/30>=0                     (18)
```

for integer `a>=1`; the other boundary is symmetric, and `(0,0)` is
immediate.  The exact equality states are frozen in the companion output.

Uniform positivity of its integral is nevertheless false.  Define

```text
B={1,3,5,7,9,11,13,15,17,19,21,45}.                 (19)
```

For every odd multiplier `m`, the paired state of `mB` at `t` is exactly the
state of `B` at `mt`; multiplication by an odd number preserves which member
of the half-turn pair is low/high.  The exact full-half sweep proves

```text
q_B(t)<=3 for every t.                                (20)
```

Hence the original `q`-cubic on every pullback `mB` is simply `H_0`, with no
negative depth-four charge.  At `h=420`, any `m` coprime to `840` makes
`{840} union mB` primitive.  It is denominator-complete: `840` supplies the
even denominators and `9m,11m,13m` supply `9,11,13`.

For the especially clean choice `m=127`, tail `45m=5715` lies in the failed
fixed-clock band `(5040,6720)`.  Thus the row

```text
W=(127,381,635,889,1143,1397,1651,1905,2159,
   2413,2667,5715)                                    (21)
```

is primitive, denominator-complete, and not removed by the fixed half-turn
clock.  Exact integration gives

```text
q_max=3,
Q_3=H_0=2337561827/25869073230
         ~=0.09036125130,
integral P=-304329184657/517381464600
         ~=-0.58821045105.                            (22)
```

The bivariate dual fails because it assigns large negative values to states
such as `(12,0)` and `(11,0)`, even though those states have a free sheet.
It discards precisely the imbalance/current coordinate restored by `(5)`.
The exact audits at multipliers `1` and `19` show the same mechanism; `(21)`
adds the sharper fixed-clock scope control.

## 5. Hostile search for the original cubic

A deterministic coordinate search over modular perturbations was used only
as a scout.  Its most attractive coarse-mesh rows were recomputed by the
independent rational wall sweep; mesh aliasing at high speeds changed their
ranking, so none is presented as finite evidence.

The strongest exact denominator-complete modular row retained from that
search is

```text
h=420,
W=(1,837,841,1681,2521,3361,4201,5041,5881,
   6721,7561,8401).                                   (23)
```

It is primitive, covers every denominator `2,...,14`, and has fixed-clock
blockers `5041,5881`.  Its exact value is

```text
Q_3=
36616988743388624099778353018234375602
/842615780094330438229544808407166979923
~=0.04345632922>0.                                    (24)
```

Here the exact decomposition is

```text
B_3(E) ~= 0.00514505287,
nonlinear current rebate ~= 0.03831127635.             (25)
```

The truncated-current lower bound (first line of `(8)`) is already positive,
about `0.02121679040`, but the plain variance certificate `(9)` is negative,
about `-0.07697858009`.  This sharply separates three levels of information:

In order of increasing retained information, the three levels are

```text
full Brownian pair energy
    -> add the anchor-strip sidecar to recover core energy
    -> retain the truncated/nonlinear current profile.  (26)
```

The arrows are information refinements, not numerical inequalities; in this
row the full-half energy is larger than the core energy.

No exact row with `Q_3<=0` was found.  This is sampled search evidence only;
the inequality `(1)` remains open.

## 6. Connection contract and next test

| field | content |
|---|---|
| source | paired tail occupations `(a,b)` on the minority anchor core |
| target | `Q_3=H_0-H_4-4H_5-10H_6` |
| map | one-sheet Bonferroni score plus absolute nonlinear current, equation `(5)` |
| preserved | exact paired danger status, all depths, sheet imbalance, physical Haar mass |
| lost by pair covariance | anchor-strip location and all current moments beyond order two |
| sidecars | anchor-strip energy; preferably the joint `(S,|C|)` profile or its truncated current |
| decisive next test | bound `integral_E max(C^2-4,0)` or the exact rebate `(6)` after conditioning on the mod-14 Brownian residue classes |

The most promising proof split is now explicit.  Low-current balanced rows
must be controlled by the depth-four tail itself; high-current/AP rows should
pay the nonlinear rebate.  A scalar variance theorem cannot bridge both,
and another unrestricted cubic polynomial in `(a,b)` cannot restore the
minimum operation.

## Reproduction

```bash
python3 -B 04-computation/lrc14_halfturn_cubic_depth_exact_audit_cubic_depth_20260902.py \
  | diff -u 05-knowledge/results/lrc14_halfturn_cubic_depth_exact_audit_cubic_depth_20260902.out -
```

The LF-byte SHA-256 values are

```text
script  c35afe7fdf2ff40c5a8ccbb84e1e34370472e98a1b536d15f01171966d8572e1
output  c4acd92d448b8667fe5004cc145edd51408bb297e88cda714cf1817ebac71d9e
```

The sampled scout, which is not load-bearing, is

```text
04-computation/lrc14_halfturn_cubic_depth_search_cubic_depth_20260902.py
```
