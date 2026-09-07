# The concrete higher-cusp curve has abelian affine complement

**Status: PROVED ANALYTICALLY + INDEPENDENTLY AUDITED. Nori's theorem
is CITED; curve and resolution identities are FINITE-EXACT. JC(2) remains
OPEN.** September 6, 2026.

## 1. Target and recovered mechanism

Let `C subset A2_C` be the image of

```text
U(t)=t^4+t^2,             V(t)=t^6+t^5+t^2.              (1)
```

We prove that `pi_1(A2_C\C)` is abelian. Consequently it has no transitive
degree-six permutation representation in which a generic meridian acts
as a single three-cycle and three fixed labels. In particular, (1)
cannot be the whole nonproperness set of a nonautomorphic planar
polynomial Keller map.

The starting mechanism is the audited
[odd-cusp divisor spectrum](planar_jc48_sep06_odd_cusp.md), whose actual
sheet transport and Euler budget force precisely that degree-six local
passport for (1). The local passport is a genuine abstract possibility,
so it cannot itself exclude (1). The named hostile is its five-sheet
nontrivial block plus one fixed label. The corrected near miss is the
initial count of four nodes, repaired to **five** by the missing pair sum
`p=1`. The underused sidecar is the full divisor at infinity. It is
essential for transporting a projective-complement theorem to the actual
affine complement.

The concept board is: actual local passport; vertical discriminant;
resolution multiplicities; retained line at infinity; self-intersection;
and global monodromy. A targeted current-repository search found no Nori
application to (1). We make no external priority claim. Instead of
extracting numerical braids, we apply the following primary result.

**CITED Nori.** Madhav V. Nori,
[*Zariski's conjecture and related problems*](https://www.numdam.org/article/ASENS_1983_4_16_2_305_0.pdf),
*Ann. Sci. École Norm. Sup.* (4) **16** (1983), 305–344,
[DOI 10.24033/asens.1450](https://doi.org/10.24033/asens.1450),
Proposition 3.27, journal p. 331, with the surface hypotheses also stated
in Introduction II, p. 306. For curves `D,E` meeting transversely in a
smooth projective surface `X`, with `D` nodal and every irreducible
component `D_i` satisfying `D_i^2>2r(D_i)`, the kernel of

```text
pi_1(X\(D union E)) -> pi_1(X\E)
```

is abelian. Here `r(D_i)` counts its nodes. The actual primary PDF
statement and hypotheses were read. We use neither a theorem limited
to ordinary cusps in the original plane nor an assertion that any
projective abelian-complement statement automatically applies affinely.

## 2. Exact normalization and complete finite singularity inventory

The map in (1) is finite onto its image: `t` satisfies the monic equation
`t^4+t^2-U=0`. Its divided differences are

```text
N(s,t)=(U(s)-U(t))/(s-t)=(s+t)(s^2+t^2+1),
M(s,t)=(V(s)-V(t))/(s-t).
```

Set `p=s+t`, `q=st`. On the branch `p=0`,
`M(s,-s)=s^4`, so there is no off-diagonal collision there. Near its
unique point `(s,t)=(0,0)`, `s^2+t^2+1` is a unit, and the pair algebra
is exactly `(s+t,s^4)`, of length four. On `p!=0`, the first equation
forces `q=(p^2+1)/2`, and the second gives

```text
4M=H(p)=-(p-1)(p^4+2p^3+4p^2+8p+1),
H(0)=1,      Disc H=-11776000,
Res_p(H,p^2+2)=123.                                     (2)
```

Thus there are five distinct nonzero pair sums. For each, `s,t` are the
two roots of `z^2-pz+(p^2+1)/2`; their discriminant `-p^2-2` is nonzero
by (2). This produces ten distinct points of the ordered pair scheme.
The independent parameter projection is

```text
Res_s(N,M)=-t^4 J(t),
J=(t^2-t+1)(2t^8+4t^7+8t^6+6t^5+12t^4+10t^3+15t^2+5t+5),
Disc J=-166571520000000000.                              (3)
```

In particular the pair scheme is zero-dimensional, so the finite map
has generic degree one: otherwise a generic image point would contribute
off-diagonal pairs in a one-dimensional family. Since `A1` is normal,
(1) is the normalization. Its intrinsic pole pair is `(4,6)`.

For completeness, no three distinct parameters have the same image.
If the common image is `(a,b)`, the remainder modulo `U-a` is

```text
V-b = -t^3+(a+2)t^2+a t-a-b  mod (U-a).                 (4)
```

A three-parameter common fibre would force this nonzero cubic to divide
the quartic `U-a`. Put `A=a+2`; comparison with
`(t+A)(t^3-A t^2-a t+a+b)` then forces

```text
f(A)=A^2+A-1=0,       j(A)=(A-2)(A^2+1)=0.
(A^2-A-1) f(A)-(A+2) j(A)=5.                            (5)
```

This is impossible. A four-parameter common fibre is already impossible
from the nonzero cubic in (4).

The derivative gcd is exactly `gcd(U',V')=t`. The only critical parameter
is therefore zero. Its fibre has `gcd(U,V)=t^2`, so its image has no other
parameter. In local target coordinates `x=U`, `y=V-U+U^2`,

```text
x=t^2+t^4,             y=t^5+3t^6+t^8.                  (6)
```

This is the ordinary analytic `(2,5)` branch. One can set the first
coordinate equal to a squared local parameter and remove the even part
of the second; its first odd exponent is five. Equivalently, the
Weierstrass quadratic completes to the `A4` normal form. The coordinate
change in (6) is used only as a local chart for the resolution. We do not
change the global projective embedding when computing its degree below.

Let `T(s,t)=U'(s)V'(t)-V'(s)U'(t)`. A direct exact Groebner calculation,
independent of the pair-sum quotient, gives

```text
(N,M,T)=(s+t,t^4) in Q[s,t].                             (7)
```

Thus every off-diagonal collision has two distinct tangent lines. By
(5) and the critical-point analysis, they give exactly five distinct
ordinary node images, with no hidden cusp collision or additional
singularity. Away from critical or multiple normalization fibres, an
immersed single branch is smooth. The root `p=1` in (2) gives the named
control: `t^2-t+1=0` maps to the node `(-1,1)`.

## 3. The projective curve and its line at infinity

The homogeneous normalization map is

```text
[S:T] |-> [S^2 T^4+S^4 T^2 : T^6+S T^5+S^4 T^2 : S^6].  (8)
```

It has no common zero: for `S!=0` the last coordinate is nonzero, and
for `S=0` the middle coordinate is nonzero. Birationality from §2
therefore makes the image a projective sextic `Cbar`. There is a unique
point at infinity, `[0:1:0]`, with a single normalization preimage.
Put `w=S/T` there, and use the local target coordinates
`X=U/V`, `Z=1/V`. The line at infinity is `Z=0`, and

```text
X=(w^2+w^4)/(1+w+w^4),       Z=w^6/(1+w+w^4),
X=w^2+O(w^3),                Z=w^6+O(w^7),
Z-X^3=2w^7+O(w^8).                                    (9)
```

Thus the infinity branch has analytic type `(2,7)`. The line at infinity
has contact order six, so it must be tracked during the resolution;
discarding it would compute the wrong complement.

We now resolve the two cusps while leaving all five nodes unchanged.
The following rows are the successive local coordinates of the strict
curve. Each arrow is one blowup at the displayed origin. The orders and
leading coefficients are checked as rational functions in the original
parameter, rather than inferred from a numerical Puiseux fit.

At the finite cusp, start with `(x,y)` from (6):

| Chart | Coordinates | Parameter orders | Boundary through the point |
|---|---|---|---|
| 0 | `(x,y)` | `(2,5)` | none |
| 1 | `(x,y/x)` | `(2,3)` | new exceptional: first coordinate zero |
| 2 | `(x,y/x^2)` | `(2,1)` | new exceptional: first coordinate zero |
| 3 | `(x^3/y,y/x^2)` | `(1,1)` | two coordinate axes |
| 4 | `(x^5/y^2-1,y/x^2)` | `(1,1)` | only the new exceptional: second coordinate zero |

In chart 3 the tangent ratio of the first coordinate to the second is
one. The last blowup separates the strict branch from both old exceptional
curves; it meets only the new exceptional, transversely, since its second
coordinate is `t+O(t^2)`. The multiplicities at the four blowup centres
are `(2,2,1,1)`, reducing self-intersection by ten.

At infinity, put `rho=Z/X^3-1=2w+O(w^2)`:

| Chart | Coordinates | Parameter orders | Boundary through the point |
|---|---|---|---|
| 0 | `(X,Z)` | `(2,6)` | line `Z=0` |
| 1 | `(X,Z/X)` | `(2,4)` | new exceptional and strict line are the two axes |
| 2 | `(X,Z/X^2)` | `(2,2)` | new exceptional and strict line are the two axes |
| 3 | `(X,rho)` | `(2,1)` | only the new exceptional: first coordinate zero |
| 4 | `(X/rho,rho)` | `(1,1)` | two exceptional coordinate axes |
| 5 | `(X/rho^2-1/4,rho)` | `(1,1)` | only the new exceptional: second coordinate zero |

In chart 2 the tangent is the diagonal. After the third blowup the
strict line at infinity is `rho=-1`, disjoint from the branch's point
`rho=0`. It remains disjoint through the last two blowups. In chart 4
the tangent ratio is `1/4`; the final blowup separates both old exceptional
curves. The second coordinate is `2w+O(w^2)`, so the final intersection
with the new exceptional is transverse. The multiplicities at the five
centres are `(2,2,2,1,1)`, reducing self-intersection by fourteen.

The ordinary analytic cusp identification is not doing hidden resolution
work here: the displayed charts directly retain the original line and
all exceptional divisors through the relevant point.

## 4. Apply the global complement theorem

Let `X` be the smooth projective surface after these nine blowups. Let
`D` be the strict transform of `Cbar`, and let `E` be the **reduced union
of every exceptional curve and the strict transform of the line at
infinity**. Then:

- `D` is irreducible and nodal with exactly five ordinary nodes.
- `D` meets `E` only at the two final exceptional curves, once each,
  transversely and away from the crossings of components of `E`.
- Its self-intersection is
  `D^2=6^2-(4+4+1+1)-(4+4+4+1+1)=12>10=2r(D)`.
- Removing `E` removes both the fibre over the finite cusp and the full
  divisor over the line at infinity. The blowdown gives isomorphisms
  `X\E = A2\{(0,0)}` and `X\(D union E)=A2\C`.

The target of Nori's kernel map is simply connected: complex affine
two-space minus a point retracts onto the real sphere
`S^3`. Thus the kernel in Proposition 3.27 is the entire group
`pi_1(A2\C)`, and that group is abelian.

An abelian permutation group acting transitively on a finite set cannot
contain a nonidentity permutation fixing a label. Indeed, if `g` fixes
one label, commutation shows it fixes every translate of that label,
hence every label. Therefore the three-cycle-plus-three-fixed meridian
required by the odd-cusp passport is incompatible with a connected
degree-six cover of this complement. This answers the specific global
representation question: **no such transitive representation exists**.

There is also an independent general way to finish the Keller implication.
For a hypothetical polynomial Keller map with this whole support, the
proper-locus cover is connected, and a generic smooth-support meridian
fixes at least one actual affine sheet. The latter follows from the
positive generic actual count proved in the
[infinity supplier](planar_jc48_sep06_infinity.md). By abelianness and
transitivity its action is the identity. Meridians normally generate
the complement group because the ambient `A2` is simply connected:
fill a loop by a disk and perturb the disk to meet the smooth part of
the curve transversely, avoiding its finite singular set. Removing
small disks around those intersections expresses the loop as a product
of conjugate meridians. Thus the entire monodromy is trivial, forcing
degree one, a contradiction to nonautomorphy. This finish does not need
the odd-cusp degree spectrum, although that spectrum was the source of
the target question.

## 5. Retained sidecars, controls, and limits

For future braid work the source independently computes the full implicit
sextic `F(u,v)` and retains its vertical discriminant:

```text
Disc_v F=-16 u^5 (u+1)^2 (4u+1)^2
               (4u^4+40u^3+120u^2+125u+25)^2.           (10)
```

The special first coordinates are the cusp value `u=0`, the simultaneous
two smooth vertical tangencies at `u=-1/4`, and the five node values.
No braid word is inferred just from (10). The successful connection is
instead:

| Coordinate | Content |
|---|---|
| Source | the actual projective normalization, including its line at infinity |
| Target | Nori's nodal strict-transform criterion |
| Map | the nine explicitly charted blowups, with all exceptional curves retained in `E` |
| Preserved predicate | the affine complement is unchanged after removing `D union E` |
| Lost information | individual local braid words and the precise nonabelian presentation before the criterion applies |
| Required sidecar | `X\E` and its fundamental group, transverse `D/E` intersections, and the strict numerical inequality |
| Cheapest hostile | replacing five nodes by six would give equality `12=2r`, outside the cited theorem |

The separate local five-cycle passport remains valid abstractly; this
curve-specific global obstruction does not retract it. Resolving the
curve while forgetting `E` or the original line would change the space,
so a projective-complement calculation alone would not suffice. We prove
no statement here about every intrinsic `(4,6)` curve, another cusp
parameter family, or the full planar Jacobian conjecture.

The source
[planar_jc48_sep06_global_curve.py](../../04-computation/planar_jc48_sep06_global_curve.py)
has a single explicitly declared curve as its universe. It checks the
ordinary pair quotient and independent elimination, no triple fibres,
critical and tangent controls, projective coordinates, every blowup map
and leading jet, the two multiplicity lists, and the final strict
intersection inequality. The genus identity `2+3+5=10` is an independent
inventory check, not a substitute for singularity classification.
There is no numerical root search, braid sampling, or finite census of
Keller maps. All checks remain active under optimized Python.

```sh
python3 04-computation/planar_jc48_sep06_global_curve.py
python3 -O 04-computation/planar_jc48_sep06_global_curve.py
```

Normal and optimized output are byte-identical, with **56 always-active
gates**. Source and output are frozen.

```text
source SHA256 6a5ef57de9491090c013ced57355d888e5fb7959cd45c8f78fb7af7fabf0478d
output SHA256 c7016a7d0d45b28f777da16c784373712507128f5e7a1894361a262c0fe768cc
```

The saved output is
[planar_jc48_sep06_global_curve.out](planar_jc48_sep06_global_curve.out).
The topological theorem, the two complement identifications, and the actual
Keller consequence require the analytic proof above; finite symbolic
checks alone do not establish them.

The complete independent
[analytic and source audit](planar_jc48_sep06_global_curve_audit.md)
passes the actual primary Nori statement, the pair/triple/critical
inventory, both resolution chains including the line at infinity,
transversality, the two complement identifications, and the exact
monodromy consequence. Independent normal/optimized/frozen replay
matches all 56 gates and both raw hashes. No correction was requested.
Root separately checked the infinity chain and the self-intersection
calculation before the final audit.

```text
independent audit SHA256 510ddb6a3f78b097e39164affe44a1da646fe0180a01939ea808ba7785dddd44
```
