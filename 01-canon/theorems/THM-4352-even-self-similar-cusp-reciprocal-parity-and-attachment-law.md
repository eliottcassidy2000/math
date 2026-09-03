---
id: THM-4352
title: "Even self-similar cusp reciprocal parity and attachment law"
status: >
  PROVED FORMAL-LOCAL + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.
  For m=2g and 1<=r<m, the normalized self-similar tail is connected,
  has two marked infinity points, genus floor((m-r-1)/2), and persistent
  delta floor(r/2), whose intrinsic sum is g-1. Reciprocity inverts the
  Newton slope but has an exact parity defect. A graph contribution +1 is
  proved only under the explicit hypothesis that the two distinct marks
  attach nodally into the same connected complement. Oriented types form
  odd-length square blocks and reciprocal orbits form triangular blocks.
  No global attachment, planar-Jacobian, or seam-entry claim is made.
source: root + next-wildcard/even-cusp-referee agents, 2026-09-02
depends_on: []
related:
  - THM-4341-odd-self-similar-cusp-reciprocal-tail-duality
  - THM-4344-clean-cubic-infinity-exit-planar-jacobian-extinction
  - THM-4351-double-normal-owner-zero-cubic-infinity-exit-planar-jacobian-extinction
mistake_firewall:
  - MISTAKE-540
primary_script: 04-computation/even_self_similar_cusp_reciprocal_parity_attachment_thm4352.py
primary_output: 05-knowledge/results/even_self_similar_cusp_reciprocal_parity_attachment_thm4352.out
primary_script_sha256: 24616ce61de3f8684366f16dd306bfc795d1507dfc641048edc03c9e336500e0
primary_output_sha256: 0e7298c0d2b2233328a0d79aa4cadade61bb9fa972b002e8f78cd1a84dcacbb0
referee_script: 04-computation/even_self_similar_cusp_reciprocal_parity_attachment_independent_referee_thm4352.py
referee_output: 05-knowledge/results/even_self_similar_cusp_reciprocal_parity_attachment_independent_referee_thm4352.out
referee_script_sha256: 55666a09c494be9cdbd7dc76576a9b55c232b8041d69e19f228ca6e92f27d305
referee_output_sha256: 6de1de07e244d356be77697774f7de4654e9cdc7e335708220559b823f77c37b
hash_basis: raw LF bytes
audit: >
  PASS AFTER INCIDENCE AND VALUATION-SCOPE REPAIRS. The primary checks
  250,000 oriented types through g=500 and 3,751,004 exact identities. The
  independent symbolic referee checks 397,940 identities through g=90 over
  Q(c), including squarefreeness, both infinity charts, Riemann--Hurwitz,
  the conductor, primitive displayed valuations, graph hostiles, reciprocal
  excess, and closed index inverses. The graph unit is conditional, B is
  integral in the threshold statement, and d*r/(m-r) is called a possibly
  fractional pre-clearing scale. Normal and optimized streams byte-match
  both frozen outputs.
---

# THM-4352 -- Even self-similar cusp reciprocal parity and attachment law

**PROVED FORMAL-LOCAL + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.
THE INTRINSIC DEFICIT, PARITY LAW, AND INDEXINGS ARE UNCONDITIONAL UNDER THE
DISPLAYED LOCAL HYPOTHESES. THE GRAPH `+1` IS CONDITIONAL ON THE TWO MARKED
ATTACHMENTS. NO PLANAR-JACOBIAN EXTINCTION OR GLOBAL SEAM CLAIM IS MADE.**

## 1. Statement and inheritance

Let `k` be algebraically closed of characteristic zero. Put `m=2g>=2`, take
`1<=r<m`, and let `u(rho)` be a unit of `k[[rho]]` with `c=u(0)!=0`.
Consider the completed one-parameter family

```text
C_(m,r): x^2=z^m-(tz)^r u(tz)
        =z^r[z^(m-r)-t^r u(tz)].                        (1)
```

After the honest base change and normalized weighted chart in Section 2,
the reduced exceptional function field is

```text
E_(m,r): Y0^2=Z^epsilon[Z^(m-r)-c],
epsilon=r mod 2.                                        (2)
```

Its smooth projective model is connected, has exactly two marked points
above `Z=infinity`, and has genus

```text
gamma_r=floor((m-r-1)/2).                               (3)
```

The normalized-away horizontal singularity is `A_(r-1)` with

```text
delta_r=floor(r/2),              gamma_r+delta_r=g-1.    (4)
```

The final identity is intrinsic. If, in addition, a semistable nodal global
completion attaches the two distinct infinity marks by two distinct nodes
to vertices in the same connected complementary subcurve, then the tail
adds one graph cycle and

```text
gamma_r+delta_r+1=g=delta(A_(m-1)).                      (5)
```

If the marks land in different connected components, or the construction
retains only a single attachment edge, the graph increment is zero; `(5)` is
then false. Identifying the two tail points to each other is a different
self-node operation and is not covered by this sentence.

The inheritance pass was:

- closest proved mechanism: THM-4341's odd self-similar cusp classifier;
- canonical hostile: THM-4344's even `A5` tail, where normalization alone is
  one genus unit short and two attachment addresses repair it;
- corrected near miss: replacing odd `m` by even `m` does not preserve the
  reciprocal genus/delta swap;
- least-used sidecar: the pair of marked infinity points, which is invisible
  in the affine tail equation but determines the graph correction.

The Anchor / Niche / Wildcard board was

```text
tail function field | parity | conductor | reciprocal ray
two infinity marks | graph incidence | square block | triangular quotient. (6)
```

## 2. Honest base change, normalization, and primitivity

Make the finite dominating base change and substitutions

```text
t=tau^[2(m-r)],       z=tau^[2r]Z,       x=tau^[mr]Y.    (7)
```

All terms in `(1)` have `tau`-weight `2mr`. Dividing by this common power
gives

```text
Y^2=Z^r[Z^(m-r)-u(tau^(2m)Z)].                          (8)
```

Write `r=2q+epsilon`, `epsilon in {0,1}`, and set `Y0=Y/Z^q`. The normalized
total chart is

```text
Y0^2=Z^epsilon[Z^(m-r)-u(tau^(2m)Z)],                   (9)
```

whose reduced exceptional curve at `tau=0` is `(2)`.

In this displayed base-changed total space, the divisorial valuation on
`(tau,z,x)` is

```text
(1,2r,mr).                                               (10)
```

Its first coordinate is one, hence the vector is primitive. Therefore `(2)`
is the exceptional function field in this chart, not the genus of a
positive-degree root-coordinate cover. This does not assert descent to an
un-base-changed weighted quotient.

## 3. Squarefreeness, the two ends, and the intrinsic deficit

The branch polynomial in `(2)` has degree

```text
d_branch=m-r+epsilon.                                   (11)
```

Because `m` is even, `(11)` is even and at least two. The nonzero roots of
`Z^(m-r)-c` are simple in characteristic zero; if `epsilon=1`, zero is an
additional simple root. Thus the polynomial is squarefree and the double
cover is connected.

At infinity put

```text
w=Z^-1,             A=Y0/Z^(d_branch/2).
```

The local equation has the form

```text
A^2=1-c*w^(m-r).                                        (12)
```

At `w=0`, the two points `A=+1,-1` are distinct and smooth. There is no
ramification there. Riemann--Hurwitz for the degree-two map, with
`d_branch` simple finite branch points, gives

```text
gamma_r=(d_branch-2)/2=floor((m-r-1)/2).                (13)
```

The square `Z^(2q)` removed from `(8)` leaves an `A_(r-1)` horizontal
singularity of delta `q=floor(r/2)`. Splitting the two parities in `(13)` now
proves the intrinsic identity `(4)`.

Adding one tail vertex and two distinct edges to one connected complementary
component changes `E-V+1` by `2-1=1`, proving `(5)` under its incidence
hypothesis. If the edges join two different components, the component count
drops by one and cancels this increment. Retaining only one attachment edge
also contributes zero. A self-identification of the two tail points is a
separate loop construction. The completed germ alone does not choose among
these global incidence types.

## 4. Reciprocal slope, parity defect, and residue excess

The numerical involution

```text
iota(r)=m-r                                                 (14)
```

sends the Newton slope

```text
lambda_r=r/(m-r)              to              1/lambda_r.  (15)
```

Unlike odd `m`, it has the fixed point `r=g`. The exact replacement for the
odd genus/delta swap is

```text
gamma_r=delta_(m-r)-1_[r even],
delta_s=floor(s/2).                                      (16)
```

The parity tax in `(16)` is intrinsic and is not the graph unit in `(5)`.

Let `d` be a positive integer and write `t=sigma^d`. For an integral `B` and
a relative residue

```text
omega=sigma^kappa z^B dz/x,                             (17)
```

the balanced ray gives the common, possibly fractional, pre-clearing
`sigma`-scale

```text
ord(z)=d*r/(m-r),                 ord(x)=g*ord(z).        (18)
```

Its valuation excess beyond `kappa` is therefore

```text
Delta_r=(B+1-g)d*r/(m-r).                               (19)
```

For fixed `(m,d,B)`, reciprocity gives the exact rational identity

```text
Delta_r*Delta_(m-r)=((B+1-g)d)^2.                       (20)
```

For the `r` row, `(m-r)|d*r` makes both scales in `(18)` integral; for the
reciprocal row the corresponding condition is `r|d*(m-r)`. Without the
appropriate divisibility they are rational pre-clearing scales. For integral
`B`, the excess in `(19)` is positive for every split exactly when `B>=g`;
at `B=g-1` it is zero. A positive prefactor `kappa` can
make the full order positive below this excess threshold, so no converse for
the full order is claimed.

## 5. Conductor and three hostile boundary rows

The threshold in Section 4 matches the normalization conductor. For

```text
R=k[[z,x]]/(x^2-z^(2g)),
N=k[[z]] direct-sum k[[z]],                              (21)
```

normalization sends `a(z)+x*b(z)` to

```text
(a+z^g b, a-z^g b).                                     (22)
```

Since two is invertible, its image is precisely the pairs congruent modulo
`z^g`. Hence the conductor is `z^gN`; `z^(g-1)(1,0)` proves sharpness.

Three small rows block the natural overstatements:

```text
m=6,r=2: E is Y0^2=Z^4-c, gamma=1, delta_r=1,
         but delta_(m-r)=2; the naive swap fails.

m=4,r=2: r is fixed by reciprocity, but gamma=0, delta_r=1;
         fixed does not mean tail and defect agree.

m=6,r=3: E is Y0^2=Z(Z^3-c), gamma=1, delta_r=1;
         only a proved same-complement two-end incidence restores 3.       (23)
```

The hypotheses `c!=0`, characteristic zero, and `r<m` are load-bearing.
At `c=0` the branch polynomial is not squarefree; in characteristic dividing
`m-r` it can be inseparable; for `r>=m` the finite-tail model changes to a
horizontal-normalization problem.

## 6. Exact square and triangular natural-number indexings

Let

```text
S_even={(g,r): g>=1, 1<=r<=2g-1}.                       (24)
```

The oriented index

```text
n(g,r)=(g-1)^2+r                                       (25)
```

is a bijection `S_even -> N_(>0)`. For fixed `g`, its block is

```text
(g-1)^2+1, ..., g^2,                                   (26)
```

of length `2g-1`, the `g`th odd number. The explicit inverse is

```text
g=ceil(sqrt(n)),                 r=n-(g-1)^2.            (27)
```

With `j=r-g`, the same index is

```text
n=2T_(g-1)+1+j,                -(g-1)<=j<=g-1.          (28)
```

Under the integer extension `T(-g)=T(g-1)`, the block center is equally
`2T(-g)+1`. Reciprocity sends `j` to `-j`; `j=0` is the fixed split.

On reciprocal orbits put `h=min(r,2g-r)` and define

```text
N(g,{r,2g-r})=T_(g-1)+h.                               (29)
```

This is a bijection from `S_even/iota` to the positive naturals. Its inverse
is

```text
g=min{q:T_q>=N},                 h=N-T_(g-1).            (30)
```

For `h<g`, one side bit recovers `r=h` or `r=2g-h`; at `h=g` the orbit is
fixed. Thus the odd square `(2h-1)^2` may be replaced by its rank `h` only
after retaining the genus block and, off the center, the side bit.

## 7. Information loss and the tournament boundary

Two exact collisions of compressed data show what the natural indices omit:

1. the reduced slopes for `(g,r)=(3,2)` and `(6,4)` both equal `1/2`, but the
   tail genera and conductor thresholds differ;
2. the unmarked curve `(2)` does not record whether its two ends attach into
   one or two complementary components, so it does not determine graph genus.

The lawful carrier is

```text
(natural block index; reciprocal side/fixed flag; coefficient c;
 two marked attachment addresses; residue buffer).                      (31)
```

Reciprocity resembles reversal of an oriented edge away from the center, but
`r=g` is a genuine fixed type and the two ends carry independent addresses.
An ordinary tournament cannot encode both features. The exact discrete
object is an involution with fixed points plus a two-ended packet, or,
equivalently, a valuation preorder with its tie polynomial and boundary
marks. This is an information-loss theorem, not a transfer of tournament
theorems to cusp geometry.

## 8. Reproduction and boundaries

Run from the repository root:

```bash
python3 -B 04-computation/even_self_similar_cusp_reciprocal_parity_attachment_thm4352.py
python3 -B -O 04-computation/even_self_similar_cusp_reciprocal_parity_attachment_thm4352.py
python3 -B 04-computation/even_self_similar_cusp_reciprocal_parity_attachment_independent_referee_thm4352.py
python3 -B -O 04-computation/even_self_similar_cusp_reciprocal_parity_attachment_independent_referee_thm4352.py
```

The primary checks all 250,000 oriented rows through `g=500`, both index
inverses, block contiguity, genus/delta, parity, reciprocal excess, and the
named hostiles in 3,751,004 exact checks. The independent symbolic referee
checks squarefreeness over `Q(c)`, both infinity points, Riemann--Hurwitz,
incidence hostiles, conductor sharpness, primitivity, and the closed inverses
through `g=90` in 397,940 checks. Normal and optimized executions byte-match
both frozen raw-LF outputs.

This theorem supplies a formal-local classifier and exact indexing scheme.
It does not assert that reciprocal germs or tails are isomorphic, identify
the global targets of the two marks, prove a Keller-map extinction, enter a
Jacobian seam, or cover `c=0`, `r>=m`, or positive characteristic.
