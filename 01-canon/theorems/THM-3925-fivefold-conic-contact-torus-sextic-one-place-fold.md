---
id: THM-3925
title: "Fivefold conic-contact torus sextics: one place forces the monogenic coordinate fold"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED. In the complete
  marked torus-sextic slice whose smooth conic and cubic meet in the divisor
  5P+Q and whose infinity line is tangent at P, unique projective infinity
  support leaves a three-parameter normal form. On its irreducible sextic
  seam the normalization denominator is exactly the pullback of the
  coefficient map's Jacobian. A nonconstant Jacobian has two zeros and gives
  exactly three normalization places at the unique infinity point; the
  resulting three-puncture affine branch is excluded by the deleted-divisor
  Jelonek gate. The only one-place seam has constant Jacobian and an explicit
  polynomial inverse; its cubic is the universal monogenic A2 cover, and
  deleting its different creates a forbidden nonconstant unit. This closes
  the marked fivefold-contact family, not all torus sextics or JC(2).
source: root / post-THM-3879 conic-contact and normalization-pole lane, 2026-08-23
audit: >
  TWO INDEPENDENT HOSTILE AUDITS PASS (jc_degree6_one_place and
  incoming_thm3925_audit/root, 2026-08-23). The
  conic-ideal quotient proves completeness of the linear ambiguity. The
  unique-infinity iff, cusp-normalization fibre product, resultant/reducible
  equivalence, absence of collapsed extra components, basepoint-free
  projective normalization, birational inverse, exact three-versus-one
  infinity-address count, polynomial fold, and monogenic different-unit
  obstruction were all independently rederived. The ordinary cusp supplies
  the reducedness sidecar in the resultant argument. The second audit found
  and repaired one routing gap: the `ell!=0` addresses are projective
  punctures meeting, not violating, THM-3920's cubic cap, so THM-3841's
  deleted-divisor/Jelonek lemma supplies the no-plane conclusion. See
  MISTAKE-467. LF-normalized normal and optimized streams match the frozen LF
  output in all 36 active gates; raw hashes and documentation checks pass.
depends_on:
  - THM-3801-cubic-etale-normalization-nonmonogenic-and-companion-sheet-gate
  - THM-3841-deleted-ramification-three-puncture-jelonek-nonentry
related:
  - THM-3844-two-cusp-polynomial-branch-quadratic-resolvent-design-gate
  - THM-3879-rational-torus-sextic-c3-packet-one-place-tradeoff
  - THM-3882-rational-dual-one-place-wronskian-projection-criterion
  - THM-3920-affine-plane-boundary-unibranch-depressed-cubic-chart-obstruction
  - THM-3928-split-affine-conic-one-place-fold-degree-barrier
  - THM-3932-infinity-component-linear-conic-torus-sextic-fold-classification
script: 04-computation/jc2_fivefold_conic_contact_torus_sextic_thm3925.py
output: 05-knowledge/results/jc2_fivefold_conic_contact_torus_sextic_thm3925.out
script_sha256: 422fd040f49fd3002c714bca4c29f24898c42ed661f243a5c19b3e769af399fb
output_sha256: 26e96492098674290a36ea956fab4b82b035d0b883aaf1b3bdbb706c0ed330bc
semantic_sha256: ba2e5cf7ec97eef34af09b2977a58f0b7fd1c51d309505100599c4f0a3365263
hash_basis: raw LF bytes
---

# THM-3925 -- Jacobian zeros are normalization poles

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.** Work over an
algebraically closed field `k` of characteristic zero. This theorem studies
one complete marked slice of torus sextics. It does not classify every torus
sextic.

Put

```text
Q2=YZ-X^2,
Q30=Z^2(Z-X),
Q3=c Q30+Q2(ell X+mY+nZ),                    c!=0,          (1)
Delta=4Q2^3-27Q3^2.                                           (2)
```

The conic `Q2=0` is smooth and has parametrization

```text
[S:T] |-> [ST:T^2:S^2].                                     (3)
```

Along it,

```text
Q3=c S^5(S-T).                                               (4)
```

Thus the conic and cubic have intersection divisor `5P+Q`, where
`P=[0:1:0]` and `Q` is the point `S=T`. Conversely, after fixing `(3)-(4)`,
every cubic in this marked contact slice has the form `(1)`: two cubics with
the same restriction to the conic differ by `Q2` times a linear form.

Assume that `Delta` has degree six and is irreducible. Then its projective
curve has only the infinity point `P` if and only if

```text
m=0,                         4+27 ell^2 !=0.                (5)
```

On this seam its normalization is rational. If `ell!=0`, exactly three
normalization places lie above `P`, and THM-3841 excludes the resulting
three-puncture branch from a dominant polynomial-plane atlas. If `ell=0`,
exactly one place lies above `P`; this case is a polynomial-coordinate
pullback of the universal depressed-cubic cusp, and the associated cubic
completion is globally monogenic. Its maximal etale open admits no dominant
polynomial-plane atlas by THM-3801.

The load-bearing bridge is the exact identity

```text
normalization denominator = pulled-back Jacobian.           (6)
```

Thus the extra infinity addresses are not accidental poles of a chosen
formula. They are the zeros of the coefficient map's Jacobian.

## 1. The unique-support normal form

On the line `Z=0`, direct restriction gives

```text
Delta(X,Y,0)
 =-X^4(4X^2+27(ell X+mY)^2).                               (7)
```

If `m!=0`, the quadratic factor has a zero away from `X=0`, because its
value at `[0:1]` is nonzero. Hence there is another projective infinity
point. If `m=0`, equation `(7)` becomes

```text
-(4+27ell^2)X^6.                                           (8)
```

The coefficient in `(8)` must be nonzero for a genuine degree-six curve;
if it vanishes, the homogeneous sextic contains the infinity line. This
proves `(5)`.

In the affine chart `Z=1`, write

```text
p=y-x^2,
q=c(1-x)+p(ell x+n).                                       (9)
```

Then `(2)` is exactly

```text
Delta=4p^3-27q^2.                                         (10)
```

The coefficient map `Phi=(p,q)` has Jacobian

```text
Jac(Phi)=c-ell p.                                          (11)
```

At the unique finite intersection `p=q=0`, this Jacobian is `c`, so the
coefficient coordinates are regular parameters and `(10)` has an ordinary
`A2` cusp there.

## 2. Rational normalization and the reducible seam

Normalize the cusp in coefficient space by

```text
p=3t^2,                         q=2t^3.                    (12)
```

Solving `(9)` for `x` gives

```text
D(t)=c-3ell t^2,
N(t)=c+3n t^2-2t^3,
x=N/D,                         y=x^2+3t^2.                 (13)
```

Conversely, away from `p=0`,

```text
t=3q/(2p).                                                  (14)
```

Equations `(12)-(14)` therefore identify a dense open of the branch with

```text
D(t)x-N(t)=0.                                               (15)
```

The only possible common factor in `(15)` is detected by

```text
Res_t(D,N)=c^2(4c-27ell(ell+n)^2).                          (16)
```

If `(16)` is nonzero, `Dx-N` is primitive and linear in `x`, hence
irreducible; `(14)` proves birationality. If `(16)` vanishes, a common root
`t_0` is necessarily nonzero. Then the whole parabola

```text
y=x^2+3t_0^2                                                (17)
```

has the constant coefficient values `(p,q)=(3t_0^2,2t_0^3)` and lies in
`Delta=0`; the sextic is reducible. Thus `(16)!=0` is exactly the
irreducible seam in this family.

## 3. The Jacobian-place identity

Homogenize `(13)` with `t=T/S`:

```text
D_h=cS^2-3ell T^2,
N_h=cS^3+3nT^2S-2T^3,

[X_h:Y_h:Z_h]
 =[N_h S D_h : N_h^2+3T^2D_h^2 : S^2D_h^2].              (18)
```

All three coordinates have degree six and satisfy the homogeneous equation
`Delta=0`. Under `(5)` and `(16)!=0`, they have no common zero. Equation
`(18)` is therefore the projective normalization map.

Its infinity divisor is transparent:

```text
Z_h=S^2D_h^2.                                               (19)
```

At `S=0`,

```text
[X_h:Y_h:Z_h]=[0:(4+27ell^2)T^6:0]=P.                     (20)
```

If `ell!=0`, `D_h` has two distinct roots. At either root, `(16)!=0` makes
`N_h` nonzero, so `(18)` again equals `[0:1:0]=P`. These are the only zeros
in `(19)`. Hence the unique projective infinity point has exactly

```text
three normalization addresses: S=0 and the two roots of D_h. (21)
```

By `(11)-(13)`,

```text
Jac(Phi) on the normalization =c-3ell t^2=D(t).            (22)
```

This proves `(6)`: the two additional addresses are precisely the two
Jacobian-zero poles needed to solve the coefficient pullback.

For `ell!=0`, the affine branch normalization is therefore `P1` minus the
three points in `(21)`.  This does **not** violate THM-3920's cubic address
cap: all three lie over projective infinity, and three meets the degree-three
bound.  Instead apply the reusable part of THM-3841.  In any normal finite
completion whose etale open admitted a dominant polynomial-plane atlas, the
deleted ramification divisor would force this branch to be a component of
the composite's Jelonek set.  Such a component is polynomially uniruled, but
`P1` minus three points admits no dominant polynomial parametrization.  Thus
the entire `ell!=0` seam is excluded from planar Keller geometry.

If `ell=0`, then `D_h=cS^2`, and `(19)` vanishes only at `S=0`; equation
`(20)` remains nonzero. The branch has affine normalization `A1` and one
projective normalization place.

## 4. One place forces the universal monogenic fold

On the one-place seam `ell=0`, the coefficient map is

```text
(x,y) |-> (p,q)=(y-x^2,c(1-x)+n(y-x^2)),                  (23)
Jac(p,q)=c.
```

It has the explicit polynomial inverse

```text
x=1+(np-q)/c,                         y=p+x^2.              (24)
```

Thus no appeal to the Jacobian conjecture is hidden here: `(23)` is visibly
a polynomial automorphism.

The natural cubic completion is the universal Cardano cover

```text
O=k[p,q,u]/(u^3-pu-q).                                     (25)
```

Eliminating `q` gives

```text
O isomorphic to k[p,u],                                    (26)
```

so the completion surface itself is `A2` and is globally monogenic over
`k[p,q]`. Its discriminant and different factor as

```text
disc(u^3-pu-q)=4p^3-27q^2,

(4p^3-27q^2)|_(q=u^3-pu)
 =(3u^2-p)^2(4p-3u^2).                                    (27)
```

The maximal etale open deletes `3u^2-p=0`, so `3u^2-p` becomes a
nonconstant unit. A dominant morphism from `A2` to that open would pull it
back to a scalar and therefore could not be dominant. Equivalently,
THM-3801 excludes every globally monogenic cubic completion of a
constant-unit plane etale open.

This is a sharp three-way trade:

```text
nonconstant coefficient Jacobian -> three projective punctures
                                 -> deleted-divisor Jelonek obstruction;
one place -> polynomial automorphism fold
          -> monogenic different-unit obstruction.                           (28)
```

## 5. Scope and reproduction

The theorem closes the complete marked family `(1)`, including every linear
change of the cubic which is invisible on the fixed conic. It does not
classify other conic contact partitions, non-torus discriminants,
nonmonogenic cubic orders with the same abstract branch, or arbitrary
rational torus sextics. Subsequent THM-3928 closes every affine singular
conic in the full classical sextic grammar: all three distinct-line fold rows
fail, and a double line factors the discriminant. THM-3932 classifies the
infinity-component conic: its fold-three family exists, while the audited
explicit member's natural cubic is monogenic. For arbitrary `q`, chosen
high-fold or double-line components remain open. In particular this theorem
does not settle `JC(2)`.

Run

```bash
python3 04-computation/jc2_fivefold_conic_contact_torus_sextic_thm3925.py
python3 -O 04-computation/jc2_fivefold_conic_contact_torus_sextic_thm3925.py
```

After platform newlines are normalized to LF, both streams must byte-match
the frozen LF file
`05-knowledge/results/jc2_fivefold_conic_contact_torus_sextic_thm3925.out`.
The companion uses exact symbolic arithmetic and no inactive `assert` gates.
