---
id: THM-4421
title: "Exceptional-quartic integral residue gap is nongeometric"
status: >
  PROVED ELEMENTARY RELATIVE TO THM-4411 + VERIFIED-EXACT. In one fixed
  integral exceptional-quartic compiler model, a parity aggregate gives a
  sharp 2Z-minus-6Z collision-period gap and the tangent-motion cokernel has
  3-torsion. A visible mixed direction, actual transgression divisibility,
  and a lawful rational target rescaling show that neither the gap nor the
  torsion is a characteristic-zero geometric obstruction. JC(2) remains open.
source: root + jc_transversality_transfer / Planar Jacobian wildcard session, 2026-09-05
depends_on:
  - THM-4411-first-order-collision-transgression-seminormal-tradeoff
related:
  - THM-4381-exceptional-quartic-seminormalization-and-conductor-fibre-classification
  - THM-4404-exceptional-quartic-descended-two-form-seminormal-cokernel
  - THM-4412-exceptional-quartic-seminormal-suspension-compensator-firewall
script: 04-computation/exceptional_quartic_integral_residue_gap_nongeometric_thm4421.py
output: 05-knowledge/results/exceptional_quartic_integral_residue_gap_nongeometric_thm4421.out
script_sha256: 934a1b8f916bbe7e966354e297a5e083d3a2ada6aa9d27efc75d734982d6409e
output_sha256: be73681ce81b7c794c82fbc9d3a59cd10cbc13287ef02da812842df82a7fd85e
hash_basis: raw LF bytes
audit: >
  PASS. A standard-library exact verifier computes all determinantal divisors,
  bounded lattice memberships, polynomial evaluation identities, collision
  motions, residue controls, and the rescaled Smith forms. All 137,041
  explicit gates remain active under optimization; normal and optimized
  outputs byte-match.
---

# THM-4421 -- exceptional-quartic integral residue gap is nongeometric

**PROVED ELEMENTARY RELATIVE TO THM-4411 + VERIFIED-EXACT.** The positive
result is conditional on one chosen integral compiler model; the invariant
result is a no-go for promoting its numerical residue gap to the
characteristic-zero geometry. This proves no chart entry, Keller pair,
`JC(2)`, or `DC(2)` statement. Both conjectures remain **OPEN**.

## 1. Integral tangent lattices

At the retained exceptional-quartic branches `x=(-1,0,1)`, use the THM-4411
tangent rows

```text
t_-=(3,-9),       t_0=(3,4),       t_+=(3,9).             (1)
```

For a common target velocity `(c,e)`, their wedge periods form

```text
       [ 9  3 ]
A  =   [-4  3 ],
       [-9  3 ]                                             (2)
```

whose column cross product is

```text
3 ell,          ell=(5,-18,13).                            (3)
```

The vector `ell` is primitive and `ell A=0`. The determinantal divisors of
`A` are `(1,3)`, hence

```text
coker(A:Z^2 -> Z^3)=Z direct_sum Z/3,
[ker(ell):im(A)]=3.                                      (4)
```

More explicitly,

```text
v in im(A)  iff  ell dot v=0 and v_-=0 (mod 3).           (5)
```

Keeping the three integral branch reparametrizations gives the full `6 x 5`
collision-motion matrix. Its determinantal divisors are

```text
(1,1,1,3,9),                                             (6)
```

so its nonzero Smith factors are `(1,1,1,3,3)` and its cokernel is
`Z direct_sum (Z/3)^2`.

## 2. Exact conditional residue gap

Let

```text
h(x)=a_0+a_1x+...+a_dx^d in Z[x],
O(h)=sum_(j odd) a_j,
E(h)=sum_(j>=2 even) a_j.                                (7)
```

THM-4411 says that the retained triple persists to first order under
`Q -> Q+s h` exactly when

```text
F(h)=5h(-1)-18h(0)+13h(1)=0.                             (8)
```

Since

```text
h(-1)=a_0+E-O,  h(0)=a_0,  h(1)=a_0+E+O,
```

direct substitution gives

```text
F(h)=8O(h)+18E(h)=2(4O(h)+9E(h)).                        (9)
```

Therefore, in this integral model,

```text
3 does not divide O(h)
  implies F(h) in 2Z minus 6Z
  implies |F(h)|>=2.                                    (10)
```

The gap is sharp at

```text
h=x^2-2x,       (O,E)=(-2,1),
(h(-1),h(0),h(1))=(3,0,-1),       F(h)=2.               (11)
```

The complete collision boundary in the two parity aggregates is

```text
F(h)=0  iff  O=9k and E=-4k for some k in Z.             (12)
```

Thus individual nonconstant monomial rays have periods `8` or `18`, but
mixed rays can cancel. This is an exact filtration fact, not a Newton-polygon
or coefficient-positivity theorem.

## 3. Visible collision-preserving hostile

The primitive mixed cancellation in `(12)` is already geometrically visible:

```text
h=9x-4x^2,
(h(-1),h(0),h(1))=(-13,0,5),
F(h)=0.                                                  (13)
```

It is not a constant and not evaluation-invisible on the retained fibre. The
literal common-motion solution is

```text
(c,e)=(-12,-16),
(v_-,v_0,v_+)=(14/3,-4,-2/3).                            (14)
```

Hence the unrestricted characteristic-zero compiler already contains a
two-ray signed cancellation preserving the collision. There is no analogue
of the LRC hypothesis that forbids the divisible residue class.

## 4. Actual transgressions kill the first torsion

The THM-4411 two-form transgression is

```text
W_h=12(h(-1),h(0),h(1)).                                 (15)
```

If `F(h)=0`, then `ell dot W_h=0`, while the first coordinate of `W_h` is
divisible by three. Criterion `(5)` therefore puts every actual integral
zero-period transgression in `im(A)`. The apparent `Z/3` wedge obstruction is
never seen by a polynomial transgression. One integral common target velocity
is

```text
c=2(h(-1)-h(1))/3,       e=2(h(-1)+h(1)).               (16)
```

The first expression is integral because `(8)` implies
`h(-1)=h(1) (mod 3)`.

On the full collision-motion lattice, write `b=a_0`. Along `(12)`, values are

```text
(h(-1),h(0),h(1))=(b-13k,b,b+5k),                       (17)
```

and the remaining integral-motion residue is

```text
rho(h)=b-k (mod 3).                                      (18)
```

It is not an obstruction on the allowed additive module. The two lawful
collision-preserving directions

```text
h_1=1:              rho=1,
h_2=9x-4x^2:        rho=-1                               (19)
```

cancel, and their sum has the integral solution

```text
h=1+9x-4x^2,
(c,e,v_-,v_0,v_+)=(-12,-12,4,-4,0).                     (20)
```

## 5. The Smith residue is not geometric

Modulo three, the endpoint tangent rows in `(1)` vanish and the middle row
has rank one:

```text
t_-=(0,0),       t_0=(0,1),       t_+=(0,0)  (mod 3).    (21)
```

Thus the only torsion prime is also a bad tangent-reduction prime. More
decisively, the lawful characteristic-zero target-coordinate change

```text
C_tilde=C/3                                               (22)
```

replaces the tangents by `(1,-9),(1,4),(1,9)`. The wedge Smith factors become
`(1,1)` and the full collision-motion factors become `(1,1,1,1,1)`: all
torsion disappears. The free relation `(3)` and the collision predicate `(8)`
remain, while the numerical size of the period and its integral residue do
not.

Consequently the `2Z minus 6Z` gap and the Smith `3`-torsion are artifacts of
the chosen integral target lattice. Only vanishing versus nonvanishing of the
free period is geometric over the characteristic-zero category used by
`JC(2)`. This rules out this specific transfer mechanism, not every possible
arithmetic or moving-graph obstruction.

## 6. Reproduction and next test

Run

```powershell
python -B 04-computation/exceptional_quartic_integral_residue_gap_nongeometric_thm4421.py
python -B -O 04-computation/exceptional_quartic_integral_residue_gap_nongeometric_thm4421.py
$env:PYTHONHASHSEED='2718'
python -B -O 04-computation/exceptional_quartic_integral_residue_gap_nongeometric_thm4421.py
```

The exact verifier checks all minors for both Smith ledgers, 117,649 bounded
wedge-lattice memberships, 16,807 polynomial formula cases, 211 actual
zero-period wedge vectors, and 425 polynomial-kernel motion solutions.

The honest next test is to compute the actual admissible moving-graph/source-
normal tangent module and ask whether its initial module contains the visible
kernel direction `(O,E)=(9,-4)` or the integral compensator `(20)`. Excluding
them is useful only if forced by the chart equations and stable under lawful
target-coordinate changes. `JC(2)` remains **OPEN**.
