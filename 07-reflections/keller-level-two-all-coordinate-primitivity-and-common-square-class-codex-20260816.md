# Level-two coordinate primitivity and the common grade-nine sign class

**Status: VERIFIED-EXACT PROOF CANDIDATE; AWAITING INDEPENDENT AUDIT.**
This note concerns only the fixed sporadic Keller map and its second iterate.
It does not identify discriminant multiplicities, extend to higher iterates or
the outside family, or change the status of `JC(2)`.

## 1. Inheritance pass and the open coordinate debt

The closest proved mechanism is
[THM-2582](../01-canon/theorems/THM-2582-odd-block-discriminant-tower-and-composite-jelonek-square-class.md):
the generic degree-nine **x-coordinate** eliminant of `F o F` has

```text
[Disc_x(F o F)] = [H] in Q(a,b,c)^*/Q(a,b,c)^{*2}.       (1)
```

The canonical hostile is that two primitive elements of different finite
algebras may have unrelated discriminants, while a nonprimitive coordinate
may have a lower-degree minimal polynomial and an eliminant with repeated or
extraneous structure.  The corrected near miss was therefore to infer the
other two classes merely from the atom-level equality of the three cubic
classes.  Composition does not preserve that inference automatically.

The least-used sidecar is the basis discriminant of the *common degree-nine
algebra*.  If `x`, `y`, and `z` are each primitive elements of that algebra,
then their power bases are three bases of one vector space and their
discriminants differ by squares of basis-change determinants.  This turns the
missing coordinate calculation into a primitivity gate.

## 2. The exact specialization

Set the final target to `(a,b,c)=(1,1,1)`.  The outer inverse cubic is

```text
25 X^3 + X - 2,                                           (2)
```

so form

```text
K_1 = Q[X]/(25X^3+X-2).                                  (3)
```

Apply the exact rational inverse formulas of THM-2582 to the outer root `X`.
They produce the middle point `q=(A,B,C)` in `K_1^3`.  Its inverse cubic

```text
xi^3 + p(q) xi + r(q)                                    (4)
```

defines a rank-nine algebra

```text
E_0 = K_1[xi]/(xi^3+p(q)xi+r(q)).                         (5)
```

The script reconstructs the remaining source coordinates `eta,zeta` by the
same `N/D` inverse section used in the norm proof.  Substitution in all three
displayed coordinates of the sporadic map gives exactly `(A,B,C)` inside
`E_0`; this is a type check that the three candidate elements lie on the same
nine-sheet fibre, not three unrelated resultants.

The norm control is independently visible inside the outer algebra:

```text
Norm_(K_1/Q)(L(q))
 = 190265288239/320
 = H(1,1,1)/(64 L(1,1,1)).                               (6)
```

## 3. Three nonzero power-basis determinants

For `v` equal to `xi`, `eta`, or `zeta`, write the columns of

```text
P_v = (1,v,v^2,...,v^8)                                  (7)
```

in the fixed tower basis

```text
1,X,X^2, xi,xi X,xi X^2, xi^2,xi^2 X,xi^2 X^2.           (8)
```

Exact rational determinants give

```text
det(P_x) != 0,
det(P_y) != 0,
det(P_z) != 0.                                            (9)
```

The three characteristic polynomials are squarefree of degree nine.  Each
has rational factor-degree pattern `(3,6)` at this target.  This splitting is
lawful: the specialization is a nine-dimensional etale algebra rather than
a field, and a primitive element of an etale algebra may have a reducible
squarefree characteristic polynomial.  What matters for the generic claim is
the nonvanishing of (9).

The pinned coefficient hashes are

```text
x  af1933db4a9e9720c71ef56628ad9eb646662ba020929a726a77387fb8e04929
y  393eb575a2b74eb0662c74c6e8c390de656f7a3b84d7b0930fc3f4088e2a688a
z  6ac03f6767ef062f7b0ea76dd3f158fb4c225dd29608c74dfbf142bd908646c0.
```

As a hostile, the intermediate outer root `A=X`, viewed in `E_0`, has a
singular rank-nine power-basis matrix: it generates only the degree-three
outer algebra.  Thus the test does not certify every visible coordinate by
dimension bookkeeping alone.

## 4. Why one good fibre proves generic primitivity

Let

```text
K = Q(a,b,c),
E = K(x,y,z),                                             (10)
```

where `E/K` is the generic degree-nine function-field extension of `F o F`.
In a fixed rational tower basis of `E`, each determinant `det(P_v)` is a
rational function of `(a,b,c)`.  All denominators used above are nonzero at
the chosen target.  Since its specialized value is nonzero, the generic
rational function is not identically zero.  Hence

```text
K(v)=E for v=x,y,z.                                       (11)
```

In particular all three generic minimal polynomials have degree nine.  Any
degree-nine coordinate eliminant obtained with a different nonzero leading
coefficient differs from the monic minimal polynomial by a scalar.  For
degree nine,

```text
Disc(cP)=c^(2*9-2) Disc(P)=c^16 Disc(P),                  (12)
```

so this normalization cannot alter its square class.

## 5. The common square class

For two bases `B` and `B M` of the same finite separable algebra,

```text
Disc(B M)=det(M)^2 Disc(B).                               (13)
```

Taking the three power bases in (7) gives

```text
[Disc_x]=[Disc_y]=[Disc_z].                               (14)
```

The exact specialization also checks the stronger numerical identities

```text
Disc_y/Disc_x = (det(P_y)/det(P_x))^2,
Disc_z/Disc_x = (det(P_z)/det(P_x))^2,                    (15)
```

with ratio-ledger hash

```text
4ef93fd8a0384bf6b9bb44340194a7d8a159a8da4d13268173663f800b94c001.
```

Combining (14) with THM-2582's proved x-coordinate class (1) yields the proof
candidate

```text
[Disc_x(F o F)]
 = [Disc_y(F o F)]
 = [Disc_z(F o F)]
 = [H].                                                   (16)
```

Thus the atom-level diagonal `[-L],[-L],[-L]` survives composition as a
diagonal, but its common entry mutates to the newest grade-two image prime
`[H]`.  This is stronger than three accidental discriminant calculations:
the mechanism is one common extension with three primitive views.

## 6. Information retained and destroyed

| field | exact content |
|---|---|
| source | the generic nine-sheet algebra of the fixed `F o F` |
| target | one quadratic sign class `[H]` with three coordinate views |
| map | primitive-element power basis followed by discriminant modulo squares |
| preserved | permutation sign quotient and odd divisor class |
| destroyed | labelled sheets, full monodromy, discriminant multiplicities, boundary effectivity |
| needed sidecar | basis-change determinant for coordinate transport; valuation data for effectivity |
| hostile | the degree-three intermediate root has singular rank-nine power basis |

No tournament is intrinsic among the three coordinates: they agree in one
binary quotient.  The tetrahedral/K4 analogy becomes relevant only after a
fourth labelled view or three independent binary characters are constructed;
equation (16) supplies neither.

## 7. Reproduction and boundary

```text
python3 04-computation/keller_level_two_all_coordinate_primitivity_probe_20260816.py
python3 -O 04-computation/keller_level_two_all_coordinate_primitivity_probe_20260816.py
```

Both routes emit identical output.  Promotion should wait for an independent
reconstruction of the tower or a direct resultant/finite-field primitivity
audit.  Even after promotion, (16) is a fixed-map, level-two, square-class
theorem only.  It does not determine exact discriminant exponents, prove the
same diagonal at levels three and four, classify all Keller counterexamples,
or settle the surviving two-dimensional Jacobian conjecture.
