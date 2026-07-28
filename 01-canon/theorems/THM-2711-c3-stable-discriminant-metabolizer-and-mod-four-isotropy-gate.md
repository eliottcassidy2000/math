---
id: THM-2711
title: "C3-stable discriminant metabolizer and mod-four isotropy gate"
status: >
  PROOF-COMPLETE CANDIDATE + VERIFIED-EXACT; AWAITING INDEPENDENT HOSTILE
  AUDIT.  A full-rank integral lattice L inside a unimodular lattice P
  determines the stable metabolizer H=P/L in the discriminant group
  A_L=L*/L, and every stable metabolizer arises from a stable unimodular
  overlattice.  For Gram matrix M, A_L[2] is ker(M mod 2), but its
  discriminant pairing is the genuinely mod-four form
  lambda(x,y)=x^T M y/2 mod 2.  Consequently a quartic V4 character plane
  on a full-rank independent rational-surface boundary must be a totally
  isotropic standard plane contained in a C3-stable metabolizer.  D4
  triality has the required standard plane but its pairing is symplectic,
  so it cannot occur in such a completion.  A doubled pair of -2 triples
  gives a sharp stable-metabolizer control and the unimodular overlattice
  I_(1,6).  Geometric realization, reflection completion, non-full-rank
  boundaries, and general A4/S4/JC2/DC2 conclusions remain open.
source: a4-resolvent-next-gate-scout-2026-07-28
depends_on:
  - THM-2655-quartic-keller-resolvent-v4-quasietale-torsor-and-kummer-class-group-gate
related:
  - THM-2695-secondary-kummer-bockstein-picard-divisibility-spectrum-and-prime-alignment-boundary
  - THM-2700-danielewski-s3-resolvent-standard-plane-exclusion
  - THM-2703-c3-boundary-tree-arm-determinant-standard-plane-gate
  - THM-2708-c3-hermitian-gain-holonomy-discriminant-gate
  - THM-2672-slope-seven-carry-nerve-exact-eleven-simplex-and-root-zero-cap
script: 04-computation/c3_discriminant_metabolizer_mod4_thm2711.py
output: 05-knowledge/results/c3_discriminant_metabolizer_mod4_thm2711.out
script_sha256: 74781f52aaaffdcea2b10afcf18738ed2627518d0caa8230200bd6189cdeb96b
output_sha256: dcac0db2a240601ee2ce32e12ecf573a2220e8b31eae26e39eee3eb2e304832d
hash_basis: LF-normalized bytes
---

# THM-2711 -- C3-stable discriminant metabolizer and mod-four isotropy gate

**PROOF-COMPLETE CANDIDATE + VERIFIED-EXACT; AWAITING INDEPENDENT HOSTILE
AUDIT.**

THM-2708 detects the quartic standard plane in the two-torsion of a
boundary discriminant group.  That is necessary, but it is not yet the
condition imposed by a smooth rational completion.  The ambient Picard
lattice is unimodular.  Its finite-index position inside the boundary dual
lattice selects a self-orthogonal subgroup of the discriminant group, and
the quartic plane has to lie inside that subgroup.

This exposes a coordinate discarded by reduction modulo two:

```text
M mod 2 identifies the order-two classes and their C3 representation;
M mod 4 identifies which of those classes can coexist in an integral
unimodular overlattice.                                             (1)
```

The distinction is exact.  In the smallest hostile, `D4` triality carries
the right `C3` plane, but the plane is symplectic rather than isotropic.  In
a positive control, two copies of a three-coordinate `-2` packet pair
diagonally and do extend to an odd unimodular lattice.

## 1. Stable metabolizers are stable unimodular overlattices

Let `L` be a nonsingular integral lattice, let

```text
L*={x in L tensor Q : (x,L) is contained in Z},
A_L=L*/L,                                                   (2)
```

and equip `A_L` with its discriminant bilinear pairing

```text
b_L(x+L,y+L)=(x,y) mod Z.                                  (3)
```

Suppose a finite group `Gamma` acts on `L` by isometries.

### Metabolizer--overlattice theorem

There is a bijection

```text
{Gamma-stable integral unimodular overlattices P,
       L subset P subset L*}

       <---->

{Gamma-stable subgroups H subset A_L with H=H^perp},       (4)
```

given by

```text
P |-> H=P/L,                                               (5)
H |-> P_H={x in L*:x+L belongs to H}.                       (6)
```

In particular, every such `H` satisfies

```text
|H|^2=|A_L|=abs(det L).                                    (7)
```

### Proof

Let `P` be an integral unimodular overlattice.  Integrality makes `P/L`
isotropic in `A_L`.  More precisely,

```text
(P/L)^perp=P*/L.                                           (8)
```

Since `P` is unimodular, `P*=P`, and hence

```text
H=P/L=H^perp.                                              (9)
```

Conversely, let `H=H^perp`.  The preimage `P_H` in `(6)` is integral
because `H` is isotropic.  Applying `(8)` to this preimage gives

```text
P_H*/L=H^perp=H=P_H/L,                                    (10)
```

so `P_H*=P_H`.  Stability is preserved in both directions because every
map in `(2)`--`(6)` is `Gamma`-equivariant.  Nondegeneracy of `(3)` then
gives `(7)`.

This is stronger than asking for a particular subgroup of `A_L`: a
metabolizer remembers the entire finite-index unimodular completion.

## 2. The order-two pairing is a mod-four invariant

Choose a basis of `L` and let `M` be its symmetric integral Gram matrix.
There is a canonical isomorphism

```text
Phi:ker(M mod 2) -> A_L[2],
Phi(x)=[x/2].                                              (11)
```

Here `x` on the right denotes any integral zero-one lift.  Indeed,
`x/2` belongs to `L*` exactly when `Mx` is even.  Changing the lift by
`2u` changes `x/2` by an element of `L`; injectivity and surjectivity follow
immediately.

For `x,y in ker(M mod 2)`, formula `(3)` becomes

```text
b_L(Phi(x),Phi(y))=x^T M y/4 mod Z.                        (12)
```

The numerator in `(12)` is even.  Therefore define

```text
lambda_M(x,y)=x^T M y/2 mod 2.                             (13)
```

It is a well-defined symmetric bilinear form on `ker(M mod 2)`, invariant
under every integral symmetry of `M`, and

```text
b_L(Phi(x),Phi(y))=lambda_M(x,y)/2 mod Z.                  (14)
```

Consequently an order-two subgroup `W_0 subset A_L[2]` is totally isotropic
if and only if

```text
lambda_M|_(Phi^-1 W_0)=0.                                  (15)
```

The kernel and its `C3` module structure depend only on `M mod 2`, while
`lambda_M` depends on `M mod 4`.  This is the first coefficient layer at
which the unimodular-completion gate can be seen.

## 3. The strengthened quartic boundary gate

Let `U` be a smooth complex surface with a `C3` action and a smooth
projective rational equivariant completion

```text
U subset Xbar,                 D=Xbar\U.                   (16)
```

Assume the irreducible components of `D` have independent classes and span
a finite-index `C3`-stable lattice

```text
L subset P=Pic(Xbar).                                      (17)
```

The rational-surface Picard lattice `P` is torsion-free and unimodular.
Divisor localization gives

```text
Pic(U)=P/L=:H,                                             (18)

O(U)^*/C^*=ker(Z^{components of D} -> Pic(Xbar))=0.        (19)
```

By Section 1, `H` is a `C3`-stable metabolizer in `A_L`.

Now suppose `U` is the regular Galois-resolvent surface of a quartic
`A4` or `S4` Keller candidate in the scope of THM-2655.  Restrict the
quotient action to its `C3`.  The canonical `V4` character module is the
irreducible standard plane

```text
W=Hom(V4,C2).                                              (20)
```

Equation `(19)` kills the unit-squareclass alternative in THM-2655, so
that theorem forces a `C3`-equivariant injection

```text
W -> Pic(U)[2]=H[2].                                       (21)
```

Since `H` is isotropic, the image in `(21)` is a totally isotropic standard
plane for `lambda_M`.  Thus the exact necessary conditions are

```text
A_L has a C3-stable metabolizer H,
H[2] contains a standard plane W,
lambda_M|_W=0.                                             (22)
```

This sharpens the THM-2708 determinant test:

```text
det(B)=0                         detects some W in A_L[2];
lambda_M|_W=0                    detects local isotropy;
W subset H=H^perp                detects completion compatibility. (23)
```

The second and third lines are genuine additional gates.  Existence of an
isotropic standard plane need not imply that it extends to a stable
metabolizer.  Conversely, a stable metabolizer might exist but contain only
trivial `C3` constituents.  Hence `(22)`, rather than any one line of `(23)`,
is the exact lattice condition forced by this completion.

Two useful contrapositives are immediate:

```text
every standard plane in ker(M mod 2) is nonisotropic
  ==> no quartic carrier in this completion;

A_L has no C3-stable metabolizer containing a standard plane
  ==> no such equivariant rational completion of the quartic carrier. (24)
```

## 4. `D4` triality passes THM-2708 but fails the new gate

Take the negative `D4` matrix, with `C3` cycling the last three coordinates,

```text
M_D4=
[-2  1  1  1]
[ 1 -2  0  0]
[ 1  0 -2  0]
[ 1  0  0 -2].                                            (25)
```

It has

```text
abs(det M_D4)=4,
Smith(M_D4)=(1,1,2,2),
A_D4[2]=W.                                                 (26)
```

In the basis given by two adjacent even leaf pairs, `(13)` is

```text
lambda_D4|_W=[0 1]
                  [1 0].                                  (27)
```

Thus the unique standard plane is symplectic, not totally isotropic.  Its
three order-two lines are permuted transitively by triality.  A metabolizer
of the order-four discriminant group would have order two, so no such line
is `C3`-stable.  Therefore

```text
D4 passes det(B)=0 but has no C3-stable metabolizer.        (28)
```

This repairs the interpretation of the sharp `D4` controls in THM-2700,
THM-2703, and THM-2708.  They prove that `D4` realizes the abstract standard
module.  They do not make its root lattice the full-rank boundary lattice of
an equivariant rational completion carrying the canonical quartic plane.

## 5. Mod two is sharply insufficient

On one free `C3` orbit compare

```text
M_symp=-2 I_3,

M_iso=
[-2  2  2]
[ 2 -2  2]
[ 2  2 -2].                                               (29)
```

Both matrices reduce to zero modulo two.  Hence they have the same kernel,
the same permutation module, and the same unique standard plane

```text
W={x in F2^3:sum x_i=0}.                                  (30)
```

Nevertheless

```text
lambda_symp|_W=[0 1],       abs(det M_symp)=8,
                     [1 0]

lambda_iso|_W=0,             abs(det M_iso)=32.            (31)
```

Thus no invariant computed from `M mod 2`, including the Hermitian kernel
of THM-2708, can decide isotropy.  The second matrix is only a local-isotropy
control: its determinant is not a square, so its discriminant group has no
metabolizer at all.  This also shows why the local second line of `(23)` is
not sufficient for completion.

## 6. A sharp stable-metabolizer positive control

Let

```text
L=<h> direct_sum <a_0,a_1,a_2> direct_sum <b_0,b_1,b_2>,
Gram(L)=diag(1,-2,-2,-2,-2,-2,-2),                        (32)
```

with `C3` fixing `h` and cyclically permuting both triples.  Its
discriminant group is `(Z/2)^6`.  Put

```text
H={<(u,u)>:u in F2^3} subset A_L[2].                      (33)
```

Every two elements of `H` have zero discriminant pairing: their two triple
contributions occur twice.  Moreover

```text
|H|=8=sqrt(64)=sqrt(abs(det L)),
H=H^perp,
H is C3-stable.                                           (34)
```

The diagonal copy of the sum-zero plane in `F2^3` is a standard `W` inside
`H`.  Adjoining

```text
(a_i+b_i)/2,                i=0,1,2,                      (35)
```

gives an index-eight overlattice.  In the basis

```text
h, (a_i+b_i)/2, (a_i-b_i)/2,                              (36)
```

its Gram matrix is

```text
diag(1,-1,-1,-1,-1,-1,-1).                               (37)
```

So the overlattice is integral, odd, unimodular, `C3`-stable, and contains
the required plane in its metabolizer.  This proves the gate is algebraically
attainable in the correct rational-surface signature.  It does not realize
the seven displayed vectors by irreducible boundary curves, build a quartic
cover, or satisfy a Keller equation.

## 7. Relation to the secondary Kummer Bockstein

Both this theorem and THM-2695 require information beyond a mod-two carrier,
but their obstruction maps are different:

```text
THM-2695:  can a mu_2 torsor lift through mu_4?
           measured by Pic(U)[2]/2 Pic(U)[4];

THM-2711:  can boundary classes sit inside a unimodular Picard lattice?
           measured by the discriminant pairing M mod 4.              (38)
```

Neither invariant determines the other without an additional geometric
identification.  Likewise the switching/gain language of THM-2708 and the
carry nerve of THM-2672 share a quotient-with-sidecar pattern, but no
physical LRC homology class is produced here.  In particular, the coarse
`boundary Delta^12` of THM-2672 remains virtual after carry is forgotten;
it is not a discriminant metabolizer.

## 8. Exact companion

Run

```text
python 04-computation/c3_discriminant_metabolizer_mod4_thm2711.py
python -O 04-computation/c3_discriminant_metabolizer_mod4_thm2711.py
```

and compare both transcripts with

```text
05-knowledge/results/c3_discriminant_metabolizer_mod4_thm2711.out.
```

The script performs only exact arithmetic and uses explicit failure raises.
It checks:

- all `4^7=16,384` symmetric `C3`-invariant matrices on two free triples
  whose two diagonal circulant blocks have parameters `(a_i,b_i)` and whose
  cross block has parameters `(c_0,c_1,c_2)` in `{0,1,2,3}`;
- all `15,058` nonsingular matrices in that box, finding `5,048` with at
  least one standard plane but only `2,560` with at least one isotropic
  standard plane;
- `7,992` standard planes in total, of which only `3,840` are isotropic;
- invariance of `(13)` under `C3` on every enumerated plane;
- the mod-two-identical pair `(29)`--`(31)`;
- the determinant, Smith form, symplectic plane, and lack of stable
  metabolizer in `D4`;
- the stable self-orthogonal subgroup `(33)`, its contained standard plane,
  and the exact index-eight unimodular change of basis `(35)`--`(37)`.

Normal and optimized transcripts byte-match the stored output.

## 9. Connection contract and next decisive test

```text
source object:
  a full-rank C3-stable boundary lattice L in a rational surface;

target object:
  the discriminant module A_L together with its mod-four pairing and a
  stable metabolizer H;

map:
  identify A_L[2] with ker(M mod 2), compute lambda_M from M mod 4, and
  identify H with Pic(Xbar)/L;

preserved predicate:
  whether the quartic V4 character plane can lie in the Picard quotient of
  a unimodular completion;

destroyed information:
  reflection, effective boundary realization, non-full-rank unit relations,
  the graph-quartic map, and every physical LRC endpoint/current;

needed sidecars:
  an equivariant geometric completion theorem and the C2 reflection action;
  a saturation/unit analysis when boundary classes are dependent;

cheapest decisive next test:
  compute the full integral boundary matrix modulo four for each surviving
  cyclic-cubic or S3 resolvent model, enumerate its C3-stable metabolizers,
  and test whether any contains an isotropic standard plane with the required
  reflection action.                                                (39)
```

The key lesson is therefore not merely that two and three coexist.  The
`C3` action finds the `F2^2` standard packet, while the next two-adic digit
decides whether that packet can inhabit a unimodular surface boundary.

QED (candidate; independent hostile audit pending).
