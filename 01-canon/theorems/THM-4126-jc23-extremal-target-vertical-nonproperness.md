---
id: THM-4126
title: "JC(2) extremal target vertical nonproperness classification"
status: >
  PROVED RELATIVE TO THM-4103/4120/4122 + CITED NGUYEN PURITY on the smooth
  nonresonant theta-only maximum-weight-eight reduced 2:3 seam. Full
  source-boundary exhaustion makes the generic affine-fibre map finite of
  degree 21; clearing
  denominators spreads finiteness away from finitely many pencil values.
  Every nonproperness component is therefore vertical and THM-4122 forces it
  to be one of the two nodal I1 affine target fibres. Its intrinsic pole pair
  is (2,3), so rho=1. The global mapping degree is 21 and the polynomial
  degrees are (2G,3G) with G>=14, hence at least (28,42); the higher-DE-weight
  branch has G>=15 and degrees at least (30,45). THM-4130 later excludes this
  smooth seam, and THM-4134/4138 empty the Delta_V wall. JC(2), the other two
  collision walls, other depth cells, and maximum residual weight at least
  nine remain OPEN.
source: planar-jacobian-squeeze / 2026-08-25
audit: >
  PASS. Two independent agents checked the full-boundary restriction,
  function-field degree, monic-equation spreading lemma, vertical-component
  classification, smooth-elliptic normalization hostile, nodal-fibre pole
  pair, and degree-floor arithmetic. The proof uses THM-4103/4120's complete
  boundary equality, not generic finiteness alone.
depends_on:
  - THM-4103-jc23-theta-boundary-ramification-and-degree-response
  - THM-4120-jc23-extremal-target-degree-twenty-one-response
  - THM-4122-planar-keller-asymptotic-width-and-resonant-shear-contraction
related:
  - THM-3996-etale-node-address-balance-cycle-and-nonproperness-dichotomy
  - THM-4124-planar-keller-integral-degree-ratio-all-vertex-shear
  - THM-4130-theta-only-extremal-seam-critical-monodromy-obstruction
  - THM-4134-delta-v-collision-wall-strict-transform-and-high-branch-exclusion
  - THM-4138-delta-v-horizontal-carrier-monodromy-exclusion
external:
  - "Nguyen Van Chau, Non-proper value set and the Jacobian condition, Ann. Polon. Math. 84 (2004), Corollary 2, arXiv:math/0305088."
---

# THM-4126 -- the extremal seam has only vertical nodal nonproperness

**PROVED RELATIVE TO THM-4103/4120/4122 + CITED NGUYEN PURITY on one reduced
seam; JC(2) OPEN.**
Work over `C` and retain every hypothesis and notation of THM-4120. In
particular, `(A,C):A2_(x,t)->A2_(U,V)` is a hypothetical nonautomorphic Keller
pair in the smooth nonresonant theta-only maximum-weight-eight survivor, and

```text
E(U,V)=V^2-U^3+(3a^2/4)U+a^3/4,             a!=0.        (1)
```

Let `S_F` be its target nonproperness set. Define the two affine nodal fibres

```text
N_0={E(U,V)=0},                 N_1={E(U,V)=a^3/2}.       (2)
```

Then

```text
S_F=N_0, N_1, or N_0 union N_1.                           (3)
```

Every irreducible component has normalization `A1` and intrinsic coordinate
pole pair `(2,3)`. In THM-4122's notation,

```text
rho_C=1,                       (d,e)=(2,3).               (4)
```

Moreover the generic mapping degree of the global polynomial map is

```text
[C(x,t):C(A,C)]=21.                                      (5)
```

Writing its polynomial degrees as

```text
(deg A,deg C)=(2G,3G),                                   (6)
```

one has

```text
G>=14,                       (deg A,deg C)>=(28,42).      (7)
```

In THM-4120's higher-maximum-`DE`-weight branch this sharpens to

```text
G>=15,                       (deg A,deg C)>=(30,45).      (8)
```

The inequalities between ordered pairs in `(7),(8)` are coordinatewise.

## 1. Full boundary exhaustion makes the generic affine map finite

Put `R=C[q]`, `K=C(q)`, and introduce the target and source pencil rings

```text
Bcal=R[U,V]/(E(U,V)-q),
Acal=R[x,t]/(E(A(x,t),C(x,t))-q).                         (9)
```

The polynomial pair induces `Bcal->Acal`. In fact

```text
Bcal isomorphic to C[U,V],              Acal isomorphic to C[x,t], (10)
```

with `R` acting through `q=E(U,V)` and `q=E(A,C)`, respectively. Their
`K`-localizations are therefore exactly the two affine generic fibres.

Write `X_Kbar=Spec(Acal tensor_R Kbar)` and base-change the target fibre
similarly. THM-4103 supplies their smooth projective completions and the
geometric extension

```text
bar(phi)_Kbar:C_q -> E_q.                                (11)
```

THM-4120 proves both

```text
deg bar(phi)_Kbar=21,
bar(phi)_Kbar^*(O)=P_AB+2(P_BC,+ +P_BC,-)+3P_CD+7P_DE
                     +3(P_EF,+ +P_EF,-).                (12)
```

The seven displayed punctures are the complete source boundary from
THM-4103. Since `A,C` are polynomial, no affine source point maps to the
projective target origin `O`. Thus `(12)` has the support equality

```text
|bar(phi)_Kbar^(-1)(O)|=C_q-X_Kbar.                      (13)
```

Removing `O` and its full inverse image from the finite projective curve map
`(11)` gives a finite morphism

```text
F_Kbar:X_Kbar -> (E_q-{O})_Kbar                          (14)
```

of degree `21`. Finiteness is fpqc-local on the target, so faithful-flat
descent along `Kbar/K` makes `Acal tensor_R K` finite over
`Bcal tensor_R K`; degree is unchanged by this base change. The induced
extension of function fields is exactly `C(x,t)/C(A,C)`: adjoining `q`
changes neither field because `q=E(A,C)` on the source and `q=E(U,V)` on the
target. This proves `(5)`.

Full boundary exhaustion is load-bearing. Merely knowing that `(11)` is a
finite projective curve map would not exclude an omitted source puncture with
finite target image, and such a puncture could support horizontal
nonproperness.

## 2. Finiteness spreads away from finitely many pencil values

The `Bcal tensor_R K`-algebra `Acal tensor_R K` is finite by `(14)` and
descent. The algebra `Acal` is generated over `Bcal` by the images of `x,t`.
Hence `x,t`
satisfy monic equations with coefficients in

```text
Bcal tensor_R K=(R-{0})^(-1)Bcal.                         (15)
```

Choose one nonzero `s(q) in R` clearing every coefficient denominator. The
same equations are monic over `Bcal_s`, so both generators are integral over
`Bcal_s`. A finite-type algebra generated by finitely many integral elements
is finite; therefore

```text
Acal_s is finite over Bcal_s.                             (16)
```

Finite morphisms are proper. Equation `(16)` implies

```text
S_F subset {s(E(U,V))=0}.                                (17)
```

Thus no irreducible component of `S_F` dominates the pencil line: every one
is contained in a fibre `E(U,V)=q_0`.

## 3. Only the two nodal fibres have affine-line normalization

After a generic source-linear change making both components monic in one
variable, Nguyen's Corollary 2 gives the load-bearing purity statement: a
nonempty planar Keller nonproperness set is a closed pure one-dimensional
curve. Let `Z` be one of its irreducible components. By `(17)`, `Z` lies in a
target fibre. Both have
dimension one, so `Z` is an irreducible component of that fibre.

THM-4120 computes

```text
Delta(q)=-432q(q-a^3/2),                 c4=36a^2!=0.     (18)
```

Every other finite fibre is a smooth irreducible affine elliptic curve. Its
smooth projective completion has genus one, so its affine normalization is
not `A1`, contradicting THM-4122. At the two zeros in `(18)`, nonzero `c4`
gives irreducible nodal `I1` fibres. Hence every `Z` is `N_0` or `N_1`.

The pair is assumed nonautomorphic, so `S_F` is nonempty: a proper Keller map
would be a finite unramified cover of simply connected `C^2`, hence an
automorphism. Purity excludes additional isolated nonproper values, so `(3)`
follows.

## 4. The nodal normalization fixes rho and the degree scale

The projective normalization of either `I1` fibre is `P1`. Its target origin
`O` is smooth and has one preimage; removing it gives `A1`. The Weierstrass
factorizations and explicit normalizations are

```text
N_0: V^2=(U+a/2)^2(U-a),
     U=w^2+a,     V=w(w^2+3a/2);
N_1: V^2=(U-a/2)^2(U+a),
     U=w^2-a,     V=w(w^2-3a/2).                         (19)
```

Thus the coordinates have pole orders

```text
(poleord_O U,poleord_O V)=(2,3).                          (20)
```

THM-4122 identifies the same pair as `(rho_C d,rho_C e)`. Equation `(20)`,
coprimality of `(d,e)`, and positivity force `(4)`. Thus `(6)` holds for an
integer `G`. THM-4120 gives `deg A>=28`, so `2G>=28` and `(7)` follows. In
its higher-weight branch the stronger `deg A>=30` gives `(8)`.

Notice the four different carriers:

```text
21       = global/function-field and generic-fibre mapping degree;
(14,21)  = pole orders at the single DE source puncture;
(42,63)  = total pullback pole-divisor degrees on the source curve;
(2,3)    = poles on a normalized target nonproperness component. (21)
```

Identifying `rho_C` with `7` or `21` is a type error.

## 5. Boundary

This theorem is an implication under hypothetical entry into THM-4120's
complete smooth seam. It does not treat the three collision walls, maximum
residual weight at least nine, other pole-depth cells, or arbitrary planar
Keller pairs. A triangular target move can change the displayed pencil and
embedded curves; THM-4124 does not by itself put this local gauge in a
target-orbit-minimal form. The theorem classifies the only possible
nonproperness set inside the seam. THM-4130 later proves that neither nodal
option is realizable by a nonautomorphic Keller pair in this seam. THM-4134
and THM-4138 later empty the separate `Delta_V=0` collision wall; the other
two collision walls remain open.
