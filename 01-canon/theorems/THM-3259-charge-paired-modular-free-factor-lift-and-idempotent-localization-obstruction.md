---
id: THM-3259
title: "Charge-paired modular free-factor lift and idempotent-localization obstruction"
status: "PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED"
source: root/modular-charge-pair/2026-08-03
depends_on:
  - THM-3250-charged-heisenberg-blowup-address-intertwiner-and-pointed-multiplicity-gate
  - THM-3243-contact-deformation-blowup-equivariance-and-full-orbit-resolution
  - THM-3245-pointed-divisor-median-cubes-saturation-band-no-go-and-z219-supplier-support
related:
  - THM-2595-modular-v4-affine-lift-dichotomy-and-six-vertex-tournament-no-go
  - THM-2596-modular-free-factor-farey-gram-owner-cocycle
  - THM-2598-quartic-v4-resolvent-torsor-and-universal-cusp-boundary
  - THM-2820-boolean-idempotent-rigidity-and-norm-top-cotangent-jet-no-go
  - THM-2839-prime-power-unit-mass-full-spectrum-and-q11-response-provenance
  - THM-2852-prime-power-orbit-spectrum-harvest-and-cayley-tournament-nonsingularity
  - THM-3253-positive-owner-mass-newton-cyclicity-and-maximal-common-heisenberg-module
script: 04-computation/heisenberg_charge_paired_modular_free_factor_thm3259.py
output: 05-knowledge/results/heisenberg_charge_paired_modular_free_factor_thm3259.out
script_sha256: 98b522c9661b4e4b07f09845e24bde1dac659ec4ed068961b1f7db1c50752458
output_sha256: bf8d9826268a765f77a5f730067668bc3acd06b6621e738bd422408fce16e96c
semantic_sha256: 0017e80fd7a3793e2712d7359ce9cba371a9993584db1708db5b9a93af0963e2
hash_basis: LF-normalized bytes
audit: >
  An independent hostile audit rederived the GL2 similitude action, generated
  group and scalar kernel, inverse-determinant charge transport, whole-block
  minimality and 26-dimensional boundary, flag formulas, THM-3250 transported
  action, and route-specific idempotent obstruction.  It also checked the
  repaired THM-3253 scope and found no mathematical defect.  Normal and
  optimized exact transcripts agree, with no assertions, floats, randomness,
  or optimization-sensitive checks.
---

# THM-3259 -- charge-paired modular free factors need a torsor, not a restoration band

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

The binary/ternary modular analogy has an exact realization on the
nonvertical Heisenberg flag torsor, but not on THM-3245's divisor restoration
band.  The order-two generator necessarily reverses nonzero central charge;
therefore the smallest whole-charge-block carrier supporting both free factors
is the pair of conjugate charge blocks.  Transport through THM-3250 gives the same
normalizing action on the paired exact-address carrier, but it remains a
cyclotomic linear construction rather than a positive physical clutch.

## 1. Idempotent restoration cannot become a nontrivial free factor

Let `J` be any restoration operator in THM-3245's internal saturation band.
It satisfies

```text
J^2=J.                                                     (1)
```

Every algebra homomorphism, invariant subquotient, or conjugation by an
invertible operator preserves `(1)`.  If the image of `J` is itself
invertible, then multiplying `J^2=J` by `J^(-1)` gives

```text
J=1.                                                       (2)
```

Thus localizing a restoration band at its idempotents collapses the inverted
restorations to the identity.  It cannot turn them into a nontrivial
involution, an order-three map, or a `C2*C3` action.  Fourier conjugation does
not change this conclusion.

THM-2839/2852's group-algebra units do not contradict `(2)`: there a Boolean
mask is reinterpreted under group convolution.  That changes the
multiplication and the source object; it is not a localization of the
pointwise/join restoration operator.  Hence the partial-cube and divisor
shadows remain exact state spaces, but idempotent localization or Fourier
conjugation cannot by itself make them modular free-factor carriers.

## 2. The Heisenberg similitude action and a uniform modular lift

Let `p` be a prime with `p=1 mod 4`, and choose `i in F_p` with `i^2=-1`.
Use symmetric coordinates on the order-`p^3` Heisenberg group:

```text
(v,z)(v',z')=(v+v',z+z'+omega(v,v')/2),
omega((x,y),(x',y'))=xy'-yx'.                             (3)
```

Every `M in GL_2(F_p)` acts by the automorphism

```text
alpha_M(v,z)=(Mv,det(M)z),                                (4)
```

because `omega(Mv,Mv')=det(M)omega(v,v')`.  Take

```text
S0=[0 -1],             S=i S0,             R=[ 0  1].     (5)
   [1  0]                                      [-1 -1]
```

Direct calculation gives

```text
S^2=1,       R^3=1,       det(S)=-1,       det(R)=1,
SR=i[1 1],  SR!=RS,       ord(SR)=4p.                      (6)
     [0 1]
```

Therefore `s->S,r->R` defines an honest finite quotient of

```text
C2*C3 = PSL_2(Z).                                         (7)
```

The projective classes of `S0,R` are the standard modular generators and
generate `PSL_2(F_p)`.  If `U=[[1,1],[0,1]]`, then

```text
(SR)^p=(iU)^p=iI.
```

It follows that the generated group is exactly

```text
Gamma_p=SL_2(F_p)<iI>,       |Gamma_p|=2p(p^2-1),
ker(Gamma_p -> PSL_2(F_p))=<iI>=C4.                        (8)
```

Indeed `(SR)^p` supplies `iI`, hence also `S0=i^(-1)S`; conversely the
displayed generators lie in the central product on the right, whose two
factors meet in `{I,-I}`.

At `p=13`, choose `i=5`.  Then

```text
S=[0 8],       R=[ 0  1],       |Gamma_13|=4,368,
  [5 0]         [12 12]
ker(Gamma_13 -> PSL_2(F_13))={1,5,8,12}I,                 (8a)
```

and `SR=5[[1,1],[0,1]]` has order `52`.  The distinguished kernel `(8)`
is cyclic; it is not the quartic matching kernel `V4` in `S4 -> S3`.

## 3. Why charge pairing is forced

Let `V_kappa` be a central-character block, so `(0,z)` acts as
`zeta^(kappa z)`.  On either permutation module write

```text
U_M [h]=[alpha_M(h)].
```

Then `rho(h)U_M=U_M rho(alpha_M^(-1)(h))`, and equation `(4)` sends charge by

```text
kappa -> det(M)^(-1) kappa.                               (9)
```

Thus `R` fixes `kappa`, while `S` swaps `kappa` and `-kappa`.  At `p=13`
the six nonzero charge orbits are

```text
{1,12},{2,11},{3,10},{4,9},{5,8},{6,7}.                  (10)
```

Consequently a single nonzero charged block cannot carry both modular free
factors.  The smallest canonical carrier assembled from whole nonzero
central-isotypic permutation blocks is

```text
V_kappa direct_sum V_-kappa.                              (11)
```

The paired central-isotypic carrier has dimension `338` on THM-3250's
permutation modules.  Abstractly, the two Stone--von Neumann representation
types have total dimension `26`; selecting one multiplicity copy in each
charge requires another frame and is not a canonical invariant subcarrier of
either permutation module.  Thus `338` is not a lower bound for arbitrary
stable subrepresentations.

This doubling is necessary, not an artifact of `(5)`.  For every odd `p`, if
a center-preserving `A in SL_2(F_p)` is an honest involution, then its minimal
polynomial divides `(X-1)(X+1)`.  It is diagonalizable, and determinant one
forces `A=I` or `A=-I`; both are projectively trivial.  The standard
projective modular involution instead squares to `-I`, hence to parity on
the Heisenberg plane.  Multiplying it by `lambda` with `lambda^2=-1` makes
the lift square to `I`, but changes its determinant to `-1` and therefore
forces the charge swap `(9)`.

The congruence assumption is sharp for this lift.  If `p=3 mod 4`, no scalar
multiple of the standard `SL_2` representative both has the same projective
class and squares to `I` over `F_p`; one must extend scalars or replace the
`PSL_2` involution by a `PGL_2` reflection.  The `p=13` construction works
over the base field precisely because `-1` is a square.

For a real packet the two coefficient matrices satisfy
`A_-kappa=conjugate(A_kappa)`.  Hence THM-3250's cyclicity determinant is
nonzero on one half exactly when it is nonzero on the other: charge pairing
doubles the carrier, not the independent determinant debt.

## 4. Exact action on the nonvertical flag torsor

Return to THM-3250's coordinates

```text
(x,y,c)(x',y',c')=(x+x',y+y',c+c'-yx').                  (12)
```

The change to `(3)` is `z=c+xy/2`, so `(4)` becomes

```text
alpha_M(x,y,c)=(X,Y,det(M)(c+xy/2)-XY/2).                 (13)
```

For `(5)` this reads

```text
alpha_S(x,y,c)=(8y,5x,-c-xy),
alpha_R(x,y,c)=(y,-x-y,c+xy+y^2/2).                       (14)
```

Identify THM-3243's regular nonvertical flag carrier by

```text
j(r,w,u)=(r,-u,w) in H_13.                                (15)
```

Transporting `(14)` through `j` gives the explicit permutations

```text
S_R(r,w,u)=(5u,ru-w,8r),
R_R(r,w,u)=(-u,w-ru+u^2/2,r-u).                           (16)
```

They satisfy

```text
S_R^2=1,       R_R^3=1,
S_R(h.x)=alpha_S(h).S_R(x),
R_R(h.x)=alpha_R(h).R_R(x).                               (17)
```

Thus the same `2/3` free factors genuinely act on one finite torsor.  This
is a based semilinear normalizing automorphism of the regular `H_13`-set,
with the group action twisted as in `(17)`.  It is not thereby an
ambient contact deformation, a Singer symmetry, or a physical LRC endpoint
map.

## 5. Transport to the paired exact-address carrier

For a fixed nonzero pair `{kappa,-kappa}`, let

```text
T_pair=T_kappa direct_sum T_-kappa                        (18)
```

be the unitary-normalized THM-3250 intertwiner from the two exact-address
blocks to the two nonvertical flag blocks.  If `U_S,U_R` are the permutation
operators from `(16)`, define

```text
S_G=T_pair^(-1) U_S T_pair,
R_G=T_pair^(-1) U_R T_pair.                               (19)
```

Then

```text
S_G^2=1,                 R_G^3=1,                         (20)
```

and these operators normalize the Heisenberg action through
`alpha_S,alpha_R`.  This is the exact charge-paired address-side modular
lift promised by `(7)`.

The construction is nevertheless cyclotomic, signed, and dependent on the
THM-3250 multiplicity frame.  Promoted THM-3253 proves that the positive
source-side owner-mass packet is cyclic in every Singer gauge and spans the
maximal common `2,041`-dimensional Heisenberg module.  What remains is not a
source cyclicity debt: no target/common physical cyclic point, owner-to-plane
relocation, positive intertwiner, terminal phase, or physical current is
supplied here.

## 6. Relation to the binary/ternary and quartic pictures

THM-2596 already proves that the binary Farey grammar consists of parabolic
words in the modular torsion generators, while the ternary Pythagorean tree
is a reduction cross-section rather than a `C3` orbit.  Equations `(5)--(20)`
give a different, literal co-occurrence object: the paired charged
Heisenberg torsor.  They do not retroactively identify branch counts with
group orders.

Likewise, THM-2595/2598's quartic packet has normal translation kernel
`V4=F_2^2` and quotient `S3`; equality of quartic and cubic-resolvent
discriminants preserves branch-divisor data but not the missing `V4` origin.
The scalar kernel `(8)` is instead `C4`, and the projective quotient is
`PSL_2(F_13)`, not `S3`.  Therefore this theorem supplies no graph-quartic,
resolvent-cubic, degree-four Keller, or Feuerbach transfer.

## 7. Evidence and scope

The exact companion verifies `(3)--(17)`, all `57,122` symplectic-pair
identities, all `2,197` flag formulas and orders, the `4,368`-element matrix
closure, its `1,092` projective classes and scalar kernel, all `2,184`
elements of `SL_2(F_13)`, the six charge pairs, and the idempotent
localization control.  Independent generic-group controls at
`p=5,13,17,29` reproduce `|Gamma_p|=2p(p^2-1)`, the scalar `C4`, and
`ord(SR)=4p`.  Normal and optimized transcripts must byte-match.

This theorem is exact finite representation and torsor anatomy.  It proves
neither a positive Boolean clutch nor an LRC row exclusion, ledger decrement,
projected cap, final rung, or LRC(14).
