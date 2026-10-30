# Point-defect phases: complete affine-flat classification and a coordinate-cost gap

Session 117. Elementary exact deductions; no literature-priority claim. This is
**affine Ramsey geometry**, not a statement about graph Ramsey clique numbers.

Let V=F_2^n, n>=3, and

    f(x) = epsilon + ell_r(x) + 1[x=v],  ell_r(x)=r.x.

All values and sums are binary. Write c=epsilon+ell_r(v), the value of the
*undefected* affine function at the defect. A monochromatic d-flat is an affine
subspace W=a+U, dim U=d, on which f is constant. Denote the Gaussian binomial
coefficient by G(m,d)=[m choose d]_2, with G(m,d)=0 outside 0<=d<=m.

## 1. Complete geometric classification

* Every singleton is monochromatic.
* A flat avoiding v is monochromatic exactly when U is contained in ker ell_r.
* A line containing v is monochromatic exactly when its nonzero direction u
  satisfies ell_r(u)=1. Its color is c+1.
* No flat of dimension at least two containing v is monochromatic.

Proof. Off v the restriction is affine, hence constant precisely when its
linear part vanishes. On a flat containing v the restriction differs at exactly
one point from an affine function. If that affine function is constant, the
one-point change destroys constancy whenever d>=1. Otherwise it has 2^(d-1)
points of each color; one change makes it constant if and only if d=1. In that
case its two values differ exactly when ell_r(u)=1. This proves every case.

## 2. Uniform counts, including colors

Let M_d be the number of monochromatic d-flats. Then M_0=2^n. For d>=1:

    r=0:   M_d = (2^(n-d)-1) G(n,d).
    r!=0:  M_d = (2^(n-d)-1) G(n-1,d) + 1[d=1] 2^(n-1).

Thus for r!=0 the line count is

    M_1 = (2^(n-1)-1)^2 + 2^(n-1)
        = 2^(2n-2)-2^(n-1)+1.

Proof. Each d-dimensional direction has exactly 2^(n-d) affine cosets, exactly
one containing v. There are G(n,d) directions when r=0, and G(n-1,d) directions
in ker ell_r when r!=0. The through-v lines contribute the 2^(n-1) vectors
u with ell_r(u)=1. These directions are all distinct over F_2.

For r=0 every positive-dimensional monochromatic flat has color epsilon. For
r!=0 and 1<=d<=n-1 the color-resolved counts are

    M_d(color c)   = G(n-1,d) (2^(n-d-1)-1),
    M_d(color c+1) = G(n-1,d) 2^(n-d-1) + 1[d=1] 2^(n-1).

For d=n both counts vanish. Singleton colors have counts 2^(n-1)-1 for c and
2^(n-1)+1 for c+1 when r!=0; for r=0 they have counts 2^n-1 for epsilon and one
for epsilon+1. The color formulas follow by dividing the cosets of each kernel
direction equally between the two affine colors, then removing the v coset.

In particular the largest monochromatic dimension is exactly n-1. There are
2^n-1 monochromatic hyperplanes for r=0, but exactly **one** for r!=0:

    {x : ell_r(x)=ell_r(v)+1}.

The n>=3 hypothesis matters here: in dimension two the exceptional through-v
lines are themselves hyperplanes.

## 3. Exact affine orbits; do not erase the output-color distinction

For fixed n>=3, there are two input-affine orbits among r!=0 defects, indexed
by c. To see sufficiency, translate v to zero, then choose an invertible linear
map carrying one nonzero functional to another. Necessity follows from the
truth-table weight: c=0 has weight 2^(n-1)+1 and c=1 has weight 2^(n-1)-1.
Alternatively the unique affine phase plus one-point defect is preserved under
input-affine transformations. Uniqueness holds because two distinct affine
functions disagree on at least 2^(n-1)>2 points, whereas two single-point
defects disagree on at most two points.

If output complementation is also allowed, those two orbits merge. In
particular the subfamily epsilon=0,v=0,r!=0 is one GL(n,2) orbit without
requiring output complementation. This is the subfamily used below.

Its Walsh coefficients, using W_f(s)=sum_x (-1)^(f(x)+s.x), satisfy

    W_f(s)=2^n 1[s=r]-2.

More generally W_f(s)=(-1)^epsilon 2^n 1[s=r]
                         -2(-1)^(epsilon+(r+s).v).

Thus every nonzero r in the normalized subfamily has the same Walsh multiset,
one coefficient 2^n-2 and all others -2. Its full colored affine-flat incidence
structure is isomorphic under an invertible input map, not merely equal in
counts or in largest-flat dimension.

## 4. The precise circuit bridge, and its limit

The companion `../circuits/POINT_DEFECT_CIRCUITS.md` proves, in the full
16-binary-gate B2 model with unrestricted sharing and one exact output,

    C(f)=n-1                  if wt(r)=0 or 1,
    C(f)=n+wt(r)-2            if wt(r)>=2,

for n>=3. Its all-dimension proof and independent finite base certificate are
separate from the affine-flat theorem above. The present section imports that
companion result; it is not a new independent proof of the gate theorem.

Within the single GL orbit epsilon=v=0,r!=0, exact gate costs run through all
integers from n-1 to 2n-2, by choosing support weights 1 through n. Consequently:

1. Full affine-flat geometry and the Walsh multiset do not determine B2 cost.
2. Any lower-bound function depending **only on the input-affine orbit** and
   valid for every member of this orbit is at most n-1 on the orbit.
3. Such a bound misses the hardest coordinate presentation by n-1 gates, or
   a factor two. This is an exact obstruction, not an asymptotic conjecture.
4. If invertible affine input reparametrization is declared free, the minimum
   ordinary B2 cost over presentations of this orbit is n-1. In the actual B2
   model that reparametrization is not free: its linear combinations must be
   computed if used as preprocessing. Input permutation/translation is much
   weaker than arbitrary invertible mixing of input bits.

This identifies a concrete missing ingredient in an affine-symmetry approach:
a useful sharp gate lower bound must retain coordinate-sensitive information,
or prove an additional lower bound for realizing the coordinate change. It
does not imply that arbitrary preprocessing always costs the displayed gap,
since alternate circuits may bypass that particular preprocessing strategy.

No NP lower-bound improvement follows: the whole family has circuits of size
at most 2n-2. These functions are exact examples of why an impressive symmetry
classification or affine Ramsey statistic is not by itself a circuit-hardness
amplifier. This does not preclude harder constructions that use the family as
one component; those constructions require their own proof.

## 5. Independent finite verification

Run `python check_affine_flats.py` and `python -O check_affine_flats.py`.
The checker builds actual F_2 subspaces and their affine cosets, evaluates f
literally on every point, and checks the classification, total counts, and
positive-dimensional color counts. It covers n=3,4,5; every r; defect points
0,1,2^n-1; and both epsilon values. Checks raise explicit exceptions, so Python
optimization cannot disable them. The all-n proofs above are not inferred
from these finite checks. Receipt: `check_affine_flats.json`.
