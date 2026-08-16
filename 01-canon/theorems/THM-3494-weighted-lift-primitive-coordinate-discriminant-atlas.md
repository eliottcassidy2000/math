---
id: THM-3494
title: "Weighted-lift primitive-coordinate discriminant atlas"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.
  For every THM-3438 weighted atom of degree n>=3, the invariant inverse
  discriminant is a constant unit times one reduced irreducible branch
  polynomial, and the actual source coordinates x and y are primitive
  elements whose eliminant discriminants differ from the full inverse
  discriminant only by nonzero squares.  Degree 3 also has an
  exact z-coordinate check; Sympy verifies degrees 3,4,5 and an independent
  FLINT route verifies degrees 3,4, including z at degree 3 and the flat-view,
  unit, and XOR hostiles.  This is an atlas for one explicit family, not a
  classification of Keller maps or a composition law.
source: codex2 derivation, 2026-08-16
audit: >
  codex independent all-degree proof audit plus a disjoint python-flint
  resultant/discriminant implementation; repaired the provisional claim that
  the branch polynomial alone always represents the full field square class,
  repaired the claim that every unsquared volume ratio has trivial square
  class, and corrected the primary output hash from CRLF replay bytes to
  LF-normalized stored bytes.
depends_on:
  - THM-3438-weighted-lift-keller-degree-spectrum
related:
  - THM-1300-jacobian-counterexample-dixmier-A3-explicit
  - THM-2473-sporadic-keller-branch-tower-depressed-trisection-anatomy
  - THM-2546-integral-coordinate-dichotomy-and-parity-lens-scope
  - THM-2582-odd-block-discriminant-tower-and-composite-jelonek-square-class
  - THM-3487-two-twenty-four-state-fibonacci-bundles-cycle-type-obstruction
  - THM-3495-level-three-sporadic-keller-norm-divisor-and-three-component-nonproperness
  - THM-3517-weighted-odd-family-m3-three-coordinate-quintic-and-sign-blind-jelonek-component
script: 04-computation/jc_weighted_primitive_coordinate_discriminant_atlas_thm3494.py
output: 05-knowledge/results/jc_weighted_primitive_coordinate_discriminant_atlas_thm3494.out
independent_script: 04-computation/jc_weighted_primitive_coordinate_discriminant_atlas_thm3494_independent_audit.py
independent_output: 05-knowledge/results/jc_weighted_primitive_coordinate_discriminant_atlas_thm3494_independent_audit.out
script_sha256: 5a84d7ae2a9c89a1a8c70eb39612f7eda41fb0f0d787366dda767a5390ae3e8c
output_sha256: 6402ddf6c7993fcd6dce0390d5191b0da3edb149259632e6d03de8a925e9fb1f
independent_script_sha256: c5326c6820d5cc955653cc3c3a9e50aab909f455313fdf025e1c8470cf168e13
independent_output_sha256: 6dda1bb7f64b8a7bec6bdca5b591e56c1a505f9e64600a2ed3222d63c356ed45
semantic_sha256: 5d327ad97829f806a16fbe9116329c32c3941fe1def2de6234155a9f8e31e075
hash_basis: LF-normalized bytes
---

# THM-3494 -- weighted-lift primitive-coordinate discriminant atlas

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

This theorem explains which part of the fixed sporadic map's three
coordinate-cubic discriminants is special to that map and which part belongs
to every separable primitive-element atlas.  It also gives the promised first
three rows of the THM-3438 all-degree family.

**Closure sidecar (THM-3517).**  The third-coordinate all-degree gate is now
proved for THM-3448's explicit cyclic weighted subfamily: an exact remainder
coefficient shows `z notin K` in every grade, and `S_n` maximality makes it
primitive.  This does not prove `z` primitive for the distinct canonical seed
used below or for every weighted seed.

## 1. Statement

Let `k` be a characteristic-zero field.  Use THM-3438's canonical seed with
`d=n-1>=2`, `b=c=1`:

```text
p_d(w)=2w-3w^2+w(1-w)(w^(d-2)-6/[d(d+1)]),
R_d(w)=integral_0^w p_d(s) ds.                         (1)
```

On the target chart `C!=0`, put

```text
P=BC,                 Q=AC^2,
K=k(P,Q,C),           T_n(w)=R_d(w)-Pw+Q.             (2)
```

Then `T_n` is irreducible of degree `n`, and THM-3438 proves that its
geometric monodromy is `S_n`.  Define

```text
D_n(P,Q)=Disc_w(T_n).                                  (3)
```

The claims are:

1. Up to a nonzero scalar in `k`, `D_n` is one reduced irreducible
   polynomial `B_n(P,Q)`.  Thus the invariant inverse equation has one
   reduced branch carrier in every grade.
2. In `K(w)`, the actual source coordinates reconstruct as

   ```text
   gamma=P-p_d(w),       x=C/gamma,
   y=(w-gamma)/C=(w-P+p_d(w))/C.                       (4)
   ```

   Both `x` and `y` are primitive: `K(x)=K(y)=K(w)`.
3. Let `M_x,M_y in K[U]` be their monic minimal polynomials.  There are
   nonzero `I_x,I_y in K` such that

   ```text
   Disc(M_x)=I_x^2 Disc(T_n^mon),
   Disc(M_y)=I_y^2 Disc(T_n^mon),                       (5)
   ```

   where `T_n^mon` is the monic normalization of `T_n`.  Equivalently, the
   raw coordinate resultants

   ```text
   E_x(X)=Res_w(T_n,X(P-p_d(w))-C),
   E_y(Y)=Res_w(T_n,CY-w+P-p_d(w))                      (6)
   ```

   have degree `n` and satisfy

   ```text
   Disc_X(E_x)/D_n in K^{*2},
   Disc_Y(E_y)/D_n in K^{*2}.                           (7)
   ```

Consequently every primitive coordinate view has the same full discriminant
square class `[D_n]=[u_n B_n]`.  At the divisor-parity level the unique odd
carrier is `B_n`; the constant class `[u_n]` must be retained in
`K^*/K^{*2}` unless it is known to be a square.  Each view's additional
squared factor is a coordinate-projection index, not another sheet or another
branch carrier.

The statement is over the generic chart and the explicit THM-3438 family.  It
does not assert that every source coordinate of every Keller map is primitive,
nor that `B_n(BC,AC^2)` is already a globally saturated Jelonek equation on
the chart `C=0`.

## 2. The branch carrier is irreducible and reduced

Differentiate `(2)`:

```text
partial_w T_n=p_d(w)-P=-gamma.                         (8)
```

The branch incidence is therefore

```text
T_n=partial_w T_n=0
  iff (P,Q)=(p_d(w),q_d(w)),
q_d(w)=w p_d(w)-R_d(w).                                (9)
```

It is the image of an irreducible affine line, so its closure is an
irreducible plane curve.  The derivative `p_d'` is nonzero in characteristic
zero.  On its dense nonvanishing locus,

```text
dQ/dP=q_d'(w)/p_d'(w)=w.                               (10)
```

The tangent slope is a rational function on the image curve and recovers the
parameter.  Hence `(9)` induces equality of function fields, is generically
one-to-one, and supplies the unique irreducible support of the resultant
`Res_w(T_n,partial_w T_n)`.

At a generic point of this curve, generic one-to-one-ness leaves exactly one
critical root, and `T_n` has exactly one double root because

```text
partial_w^2 T_n=p_d'(w)!=0.                            (11)
```

The independent coefficient `Q` occurs with derivative one, so it gives a
transverse smoothing of that double root.  The polynomial discriminant
therefore vanishes to order one.  Since the leading coefficient of `T_n` is a
nonzero constant, there is no coefficient-at-infinity factor.  This proves

```text
D_n=u_n B_n,       u_n in k*,       B_n irreducible.   (12)
```

This is the discriminant version of THM-3438's generic-transposition proof of
`S_n` monodromy.

## 3. Why the two source coordinates are primitive

First put `K_0=k(P,Q)`.  THM-3438 gives an `S_n` normal closure over `K_0`.
The extension `K_0(C)/K_0` is purely transcendental and therefore regular;
it is linearly disjoint from every finite algebraic extension of `K_0`.
Consequently irreducibility, the normal closure, and its Galois group survive
the base change to `K=K_0(C)`.  If `N/K` is that base-changed normal closure,
then

```text
Gal(N/K)=S_n,        Gal(N/K(w))=S_(n-1).              (13)
```

The point stabilizer `S_(n-1)` is maximal in `S_n`, so `K(w)/K` has no proper
intermediate fields.

Neither coordinate in `(4)` belongs to `K`.  Indeed, if `x` belonged to `K`,
then `gamma=C/x` and hence `p_d(w)=P-gamma` would belong to `K`.  That would
give a nonzero equation for `w` of degree `d=n-1`, contradicting the degree-`n`
minimal polynomial `T_n`.  Similarly, if `y` belonged to `K`, then the
nonconstant polynomial

```text
w+p_d(w)-P-Cy                                          (14)
```

would annihilate `w` and have degree `d<n`, another contradiction.  Maximality
in `(13)` now forces

```text
K(x)=K(w)=K(y).                                        (15)
```

This step is exactly what fails for an arbitrary coordinate of an imprimitive
composition tower.

## 4. The square-index identity

For a primitive element `xi` of a separable degree-`n` extension, compare the
two power bases

```text
(1,w,...,w^(n-1)),       (1,xi,...,xi^(n-1)).          (16)
```

If the change-of-basis determinant is `I_xi`, the trace Gram matrices are
congruent, so

```text
Disc(xi)=I_xi^2 Disc(w).                               (17)
```

Applying `(17)` to `x` and `y` proves `(5)`.  More explicitly, away from the
generic discriminant the conjugate values `C/(P-p_d(w_i))` and
`(w_i-P+p_d(w_i))/C` are distinct because `x` and `y` are primitive; their
denominators are nonzero because `P-p_d(w_i)=-T_n'(w_i)`.  Thus each raw
resultant in `(6)` is a nonzero scalar in `K*` times the corresponding monic
minimal polynomial, rather than a power of one.  Scaling a degree-`n`
polynomial multiplies its discriminant by the even power `2n-2`.  The same
observation removes the constant leading coefficient of `T_n`.  This proves
`(7)`.

The primitivity gate is essential.  The exact companion replaces the
coordinate by the flat base-field view `xi=P`; its resultant is

```text
(X-P)^n,                                                (18)
```

whose discriminant is zero.  A symmetric pairwise comparison without the
primitive-element sidecar can therefore lose the entire extension.

## 5. Exact degrees 3, 4, and 5

The primary resultant companion works over `Q[P,Q,C]`, reconstructs every
displayed polynomial, factors over `Q`, and obtains:

| `n` | `D_n=u_n B_n` | `B_n` terms / `(deg_P,deg_Q)` | `I_x` core | `I_y` core |
|---:|---|---|---|---|
| 3 | `u_3=-1` | 5 / `(3,2)` | 3 terms, degree 1 | 1 term, degree 1 |
| 4 | `u_4=-1/16` | 9 / `(4,3)` | 10 terms, degree 3 | 8 terms, degree 3 |
| 5 | `u_5=1/4000000` | 14 / `(5,4)` | 28 terms, degree 6 | 26 terms, degree 6 |

Here an index core removes the forced factor
`C^(n(n-1)/2)` from the square root of `(7)`.  Every listed `B_n` is one
exponent-one factor, and each displayed index core is coprime to `B_n`.
The constants are not cosmetic in the field square class: over `Q`, the first
two rows have `[D_3]=[-B_3]` and `[D_4]=[-B_4]`, while the displayed positive
scalar in the fifth-degree row is a rational square.  All three rows have the
stated odd divisor carrier `B_n`.

For `n=3`, the formulas are small enough to print:

```text
p_2(w)=-3w^2+2w,
T_3(w)=-w^3+w^2-Pw+Q,

B_3=4P^3-P^2-18PQ+27Q^2+4Q,                           (19)

Disc(E_x)/Disc(T_3)
  =C^6(9P-27Q-2)^2,
Disc(E_y)/Disc(T_3)
  =(27C^3Q)^2.                                         (20)
```

The actual `z` reconstruction is also checked in this row.  Its cubic
eliminant has 44 terms and

```text
Disc(E_z)/Disc(T_3)=I_z^2,                             (21)
```

with a 15-term, total-degree-seven index core coprime to `B_3`.  Thus the
canonical weighted cubic has three coordinate cubics with one square class,
by a calculation independent of the fixed sporadic map.

The coefficient-ledger hashes are:

```text
B_3  37562bfd2b6d06cedcfabcb4b4904db42f503f410ffa86a6534c6b24e6d7697e
B_4  62c7064c3132534306fa01216b99211a1a1bd5651226758277fa12edab971978
B_5  54742f630852a01c2b4d693e48702a83ee22f37901cefd61a63f4637f3753551

I_x: ce9f8dad...d49c454, e4cf75a5...cbc6b36, ff542f9b...e3ba674
I_y: 1c9df71a...4785acb, a0ce0c6d...88f19c1, b0ca5895...0d5c1
I_z(n=3): b4542923...be1a00
```

The full exact ledgers, eliminant hashes, signs, scalar denominators, gcd
checks, and ordinary/optimized controls are in the stored output.  The
independent audit uses FLINT's `fmpq_mpoly` representation, imports neither
Sympy nor the primary script, and reproduces the complete `n=3,4` branch,
`x`, and `y` rows, the 44-term `n=3` `z` eliminant, all displayed coefficient
hashes in those rows, branch/index coprimality, and the flat hostile.

## 6. The fixed three cubics and the infinite family

For the fixed sporadic cubic map, THM-2546 proves

```text
Disc_x=-4 S_x^2 L,
Disc_r=-4 S_r^2 L,
Disc_z=-4 S_z^2 L.                                    (22)
```

The invariant content is `[Disc]=[-L]`; `4` and the three `S_i^2` are square
index/normalization data.  Sections 3--4 explain why three primitive views of
one cubic extension must differ by squares.  What remains special in `(22)` is
the particular carrier `L`, its boundary effectivity, and the fact that all
three named source coordinates were proved primitive—not the existence of a
common square class itself.

THM-3438 then supplies the exact all-grade continuation: one explicit
`S_n`-monodromy atom and the primitive `x/y` atlas for every `n>=3`.  This does
**not** classify all maps of degree `n`.  In product grades, a weighted `S_n`
atom coexists with decomposable maps whose imprimitive block systems introduce
norm and cross-block discriminant factors.  THM-2582's fixed-map degree-nine
calculation is the canonical hostile: composition changes the odd carrier.

THM-3495's proved fixed-map level-three theorem finds a new prime image in a
norm numerator and records the full square class `[-2J]`, independently
reinforcing the constant-unit guardrail above.  This is evidence for a
divisor-orbit tower under
composition, not evidence that the weighted `B_n` are powers or iterates of
the cubic `L`.  The two constructions remain separately typed.

## 7. Four views, six edges, and the tournament/XOR boundary

Choose any primitive views `xi_1,...,xi_m` of the same separable extension and
let `v_i` be the determinant taking one reference power basis to the `xi_i`
power basis.  Label the oriented complete-graph edge by

```text
g_ij=v_j/v_i.                                          (23)
```

Then

```text
g_ji=g_ij^(-1),       g_ij g_jk=g_ik,
Disc(xi_j)=g_ij^2 Disc(xi_i).                          (24)
```

Thus the edge labelling is the exact multiplicative coboundary of the vertex
volumes `v_i`.  For four views there are exactly six edges of `K_4`, matching
the proposed size-four/six-edge carrier.  But a tournament orientation is only
a gauge choice: reversing an edge inverts its label.

There are three distinct quotients, and the distinction is load-bearing:

```text
unsquared square class: [g_ij]=[v_j]-[v_i] in K^*/K^{*2},
discriminant ratio:     Disc(xi_j)/Disc(xi_i)=g_ij^2,
chosen XOR character:   chi([g_ij])=chi([v_j])-chi([v_i]) in F_2.  (25)
```

The first is an exact square-class coboundary but its individual edges need
not vanish.  The second has trivial square class on every edge.  The third is
an exact XOR coboundary only after a homomorphism
`chi:K^*/K^{*2} -> F_2` is chosen.  The bare indicator "square versus nonsquare"
is not a homomorphism when the square-class group has rank greater than one.
The minimal rational hostile uses vertex volumes `(1,2,6)`: the edge ratios
`(2,3,6)` are all nonsquares, so the bits `(1,1,1)` violate
`b_12 XOR b_23=b_13`, although `g_12 g_23=g_13`.  In contrast the
2-adic-valuation character gives `(1,0,1)`, and the discriminant ratios
`(4,9,36)` are edgewise square-trivial.

The preserved targets are the common full discriminant square class and its
odd branch-divisor carrier; quotienting by constant units preserves only the
second.  The destroyed data are the individual projection-collision divisors
and boundary valuations.  The needed sidecar is the unsquared
basis-volume/index cochain.  This particular
cochain on the complete-graph 1-skeleton has zero `H^1` class because it is a
coboundary; its edge labels, including their square classes, may still be
nontrivial.  This does not say that every graph cochain is exact, nor does it
transport a physical LRC current.

## 8. Exact verification and open boundary

Run

```bash
python3 -B 04-computation/jc_weighted_primitive_coordinate_discriminant_atlas_thm3494.py
python3 -B -O 04-computation/jc_weighted_primitive_coordinate_discriminant_atlas_thm3494.py
python3 -B 04-computation/jc_weighted_primitive_coordinate_discriminant_atlas_thm3494_independent_audit.py
python3 -B -O 04-computation/jc_weighted_primitive_coordinate_discriminant_atlas_thm3494_independent_audit.py
```

The companion checks, without `assert`:

1. the canonical seed identities and `T_n'=p_d-P`;
2. the resultant/discriminant sign and constant leading coefficient;
3. irreducibility and exponent one of `B_3,B_4,B_5` over `Q`;
4. exact `x/y` resultants, degree `n`, nonzero discriminants, and rational
   perfect-square ratios;
5. coprimality of every low-degree branch/index pair;
6. the full `z`-coordinate cubic at `n=3`; and
7. the flat-coordinate hostile `(18)`.

Ordinary and optimized transcripts agree byte-for-byte after the platform
newline convention is fixed.  The all-degree proof is Sections 2--4; the
finite computations are independent checks, not extrapolations from the low
rows.  The audit also verified that pure-transcendental base change preserves
the `S_n` closure, that point-stabilizer maximality is used only after proving
`x,y notin K`, that raw-resultant scaling contributes only even powers, and
that the constant unit `u_n` cannot be discarded from a field square class.
The originally pinned primary output hash was the Windows-CRLF replay hash;
the frontmatter now records the actual LF-normalized stored-output hash.

Still open or outside scope:

- global saturation across `C=0` and the full Jelonek divisor for every seed;
- third-coordinate primitivity for this canonical seed or every weighted
  seed (THM-3517 closes it only for the cyclic THM-3448 subfamily);
- composition-norm recurrences for the weighted atoms;
- classification up to tame conjugacy, stable equivalence, or composition;
- any implication for `JC(2)`, `DC(1)`, `DC(2)`, LRC(14), or a physical
  tournament current.
