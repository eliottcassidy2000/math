---
id: THM-3250
title: "Charged Heisenberg blowup-address intertwiner and pointed multiplicity gate"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.
  For every nonzero central character of H_p, an explicit two-stage Fourier
  transform identifies THM-3240's p-copy exact-address carrier with
  THM-3243's regular nonvertical blowup-flag orbit.  The equally sized
  permutation carriers are not isomorphic before localization: their entire
  mismatch is the central-neutral character spectrum.  In each charged
  block all intertwiners are a frame on a p-dimensional multiplicity space,
  and
  a pointed current supplies that frame exactly when its p-by-p coefficient
  determinant is nonzero.  The bridge is unitary after normalization but is
  neither positive nor Singer-equivariant and supplies no physical clutch.
source: root/multiscale-newton-flag/low-child-flag-extension/2026-08-03
audit: >
  The assertion-independent exact companion pins THM-3240 and THM-3243;
  exhaustively checks the charged intertwining formula at p=3,5,7,13
  (4,455,516 p=13 state/character tests); verifies the neutral character
  multiplicities, Weyl orthogonality and full matrix-algebra dimensions;
  checks the normalized energy factors and exact cyclic/noncyclic determinant
  controls; and records the p=13 dimensions 169+2028=2197.  Normal and
  optimized runs byte-match the stored transcript and LF-normalized hashes
  below.  An independent hostile audit rederived both actions and their group
  law, every Fourier phase and sign, the neutral character multiplicities,
  unitary normalization, Weyl span, orbit-rank formula, pointed uniqueness,
  and all scope boundaries, and found no defect.  A post-promotion independent
  audit also verified the sharp full-module rank defect and the nonnegative
  neutral-mass obstruction in Section 4.
depends_on:
  - THM-3240-exact-address-heisenberg-clutch-on-carrier-imbalance
  - THM-3243-contact-deformation-blowup-equivariance-and-full-orbit-resolution
related:
  - THM-3234-singer-owner-compactification-and-pointed-heisenberg-carrier-gate
  - THM-3236-contact-spectrum-primitive-element-and-root-reconstruction-gate
script: 04-computation/heisenberg_charged_blowup_address_intertwiner_thm3250.py
output: 05-knowledge/results/heisenberg_charged_blowup_address_intertwiner_thm3250.out
script_sha256: 52df74d7468c78bd791eb02f5bb8f9555edab8286e5f7cb14024b5b74a61008b
output_sha256: 889b1a0cf29f5ecf883fadff41cb8c6443430998362cbffba7ecf53e9faebee6
hash_basis: LF-normalized bytes
---

# THM-3250 -- charged Heisenberg blowup-address intertwiner and pointed multiplicity gate

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-3240's exact-address carrier and THM-3243's nonvertical exceptional
flags both have `p^3` points, but they are not isomorphic `H_p`-sets.  The
first is `p` affine orbits with stabilizer `p`; the second is one regular
orbit.  This theorem identifies the precise sense in which the two carriers
nevertheless agree: every nonzero central-Fourier sector is explicitly and
unitarily isomorphic.  All of the representation mismatch lives at central
charge zero.

The surviving ambiguity is a multiplicity frame.  One pointed current kills
that ambiguity exactly when a finite coefficient matrix is invertible.  This
is a linear charged-sector bridge, not a positive physical identification.

## 1. The two carriers

Use the THM-3240 convention

```text
H_p={(x,y,c):x,y,c in F_p},
(x,y,c)(x',y',c')=(x+x',y+y',c+c'-yx').              (1)
```

The `p`-copy exact-address carrier is

```text
G={(s,t,delta):s,t,delta in F_p},                    (2)
```

with action

```text
(x,y,c).(s,t,delta)
 =(s+x,t,delta+c-ys).                                (3)
```

For `p=13`, this is exactly THM-3240's `G_delta`; `t` is its invariant
spectator coordinate.

The nonvertical part of THM-3243's full-orbit blowup boundary is

```text
R={(r,w,u):r,w,u in F_p},                            (4)
```

where `u` is the affine slope of the direction `span(1,u)`.  Differentiating
the affine action gives

```text
(x,y,c).(r,w,u)
 =(r+x,w+c-yr,u-y).                                  (5)
```

This is a regular `H_p`-set.  In contrast, `(2)` is the disjoint union of
the `p` fixed-`t` affine planes.

## 2. Central Fourier bases

Let `K=Q(zeta)` with `zeta` a primitive `p`th root of unity.  For
`kappa in F_p`, let the `kappa`-charged block mean the eigenspace on which
the central element `(0,0,c)` acts by `zeta^(kappa c)`.

In `K[G]` define

```text
E_(s,t)^kappa
 =sum_(delta in F_p) zeta^(-kappa delta)[s,t,delta].  (6)
```

For `kappa!=0` these `p^2` vectors form a basis of `K[G]_kappa`, and `(3)`
gives

```text
(x,y,c)E_(s,t)^kappa
 =zeta^[kappa(c-ys)] E_(s+x,t)^kappa.                (7)
```

In `K[R]`, Fourier-transform both the central coordinate `w` and the slope
coordinate `u`:

```text
F_(r,v)^kappa
 =sum_(w,u in F_p) zeta^(-kappa w-vu)[r,w,u].        (8)
```

Equation `(5)` gives

```text
(x,y,c)F_(r,v)^kappa
 =zeta^[kappa c-kappa yr-vy] F_(r+x,v)^kappa.        (9)
```

## 3. The explicit charged intertwiner

Fix `kappa!=0`.  Define

```text
T_kappa(E_(s,t)^kappa)
 =F_(s-kappa^(-1)t,t)^kappa.                         (10)
```

Substituting `r=s-kappa^(-1)t` and `v=t` into `(9)` gives

```text
kappa c-kappa y(s-kappa^(-1)t)-ty
 =kappa(c-ys).                                       (11)
```

The shifted first coordinate also satisfies

```text
(s+x)-kappa^(-1)t=(s-kappa^(-1)t)+x.                 (12)
```

Equations `(7)--(12)` prove

```text
T_kappa:K[G]_kappa -> K[R]_kappa                     (13)
```

is an `H_p`-equivariant isomorphism.

After embedding `K` in `C`, put

```text
Ehat=p^(-1/2)E,             Fhat=p^(-1)F.             (13a)
```

These are orthonormal bases.  The rescaled complex intertwiner

```text
That_kappa=p^(-1/2)T_kappa                              (13b)
```

sends `Ehat_(s,t)` to `Fhat_(s-kappa^(-1)t,t)`, so it is a basis
permutation and hence unitary.  Thus the algebraic map `(10)` admits a
canonical unitary normalization preserving the complete charged `l^2`
energy exactly; the unscaled `K`-linear map has squared norm factor `p`.

## 4. The neutral obstruction is complete

Write

```text
chi_(a,b)(x,y,c)=zeta^(ax+by)                         (14)
```

for the `p^2` one-dimensional characters of `H_p`.

At `kappa=0`, summing over `delta` in `(6)` kills the `y,c` action and leaves
only translation of `s`.  Hence

```text
K[G]_0
 = direct_sum_(a in F_p) chi_(a,0)^[direct_sum p].   (15)
```

On `R`, use `(8)` with `kappa=0` and then Fourier-transform `r`.  The pair
of dual variables runs over all `(a,b)`, each once:

```text
K[R]_0
 = direct_sum_(a,b in F_p) chi_(a,b).                (16)
```

Thus `(15)` and `(16)` are not isomorphic for `p>1`.  Conversely, Section 3
identifies every nonzero central block.  The entire linear representation
mismatch between the two equal-cardinality carriers is therefore the
central-neutral spectrum.

This also proves that no extension of the direct sum of the `T_kappa`,
`kappa!=0`, can be an `H_p`-isomorphism on the full permutation modules.

More sharply, every full `H_p`-linear map

```text
L:K[G] -> K[R]                                           (16a)
```

has

```text
rank(L)<= (p-1)p^2+p = p^3-p^2+p.                       (16b)
```

Indeed the charged blocks contribute at most `(p-1)p^2`.  In the neutral
block, each source character `chi_(a,0)` has multiplicity `p` but the target
contains it once, so the neutral rank is at most `p`; the target characters
`chi_(a,b)` with `b!=0` receive nothing.  Taking the isomorphisms `(13)` on
all charged blocks and a rank-one multiplicity map for every `chi_(a,0)`
attains `(16b)`.  Thus every full intertwiner has kernel and cokernel of
dimension at least

```text
p^2-p.                                                   (16c)
```

For `p=13`, the sharp maximum rank is `2041` and the exact defect on both
sides is `156`.

## 5. Multiplicity space and the exact cyclic-current gate

For fixed `kappa!=0`, let `pi_kappa` be the `p`-dimensional representation
on the `s` coordinate in `(7)`.  Translation in `s` and multiplication by
the `p` distinct phases `zeta^(-kappa ys)` generate the full matrix algebra:

```text
image(K[H_p] -> End(pi_kappa))=Mat_p(K).              (17)
```

Consequently

```text
K[G]_kappa = pi_kappa tensor K^p,
K[R]_kappa = pi_kappa tensor K^p,                     (18)
```

and every `H_p`-linear map between them is `T_kappa` followed by an arbitrary
linear map on the multiplicity factor.  In particular the isomorphisms form
a `GL_p(K)` torsor.  Formula `(10)` is one frame selected by the labelled
spectator and slope coordinates; it is not an unlabelled canonical frame.

Now let a charged current on `G` be

```text
v_A=sum_(s,t) A_(s,t) E_(s,t)^kappa,
A=(A_(s,t)) in Mat_p(K).                              (19)
```

Because of `(17)`, the linear span of its `H_p` orbit is

```text
span_K(H_p.v_A)=pi_kappa tensor rowspan(A),
dim span_K(H_p.v_A)=p rank(A).                        (20)
```

Hence the following are equivalent:

```text
v_A is cyclic in K[G]_kappa;
span_K(H_p.v_A)=K[G]_kappa;
rank(A)=p;
det(A)!=0.                                            (21)
```

Use the aligned target basis `T_kappa(E_(s,t)^kappa)` and let `v_B` have
coefficient matrix `B`.  With the convention that a multiplicity map acts
by right multiplication, a pointed intertwiner carrying `v_A` to `v_B`
must solve

```text
A U=B.                                                (22)
```

If `det(A)!=0`, `(22)` has the unique solution `U=A^(-1)B`.  This pointed
map is an isomorphism exactly when `det(B)!=0`.  Thus two cyclic pointed
currents canonically remove the multiplicity ambiguity; one noncyclic
current cannot.

The gate is sharp.  `A=I_p` is cyclic.  The all-ones matrix has every
coefficient nonzero and sees every spectator label, but has rank one, so
its orbit spans only one copy of `pi_kappa`.  Nonzero energy, full support,
and nonzero central charge do not imply cyclicity.

The determinant in `(21)` is the multiplicity analogue of THM-3236's
Krylov determinant: both decide when one observer reconstructs every copy
hidden behind a finite quotient.  Here, however, it supplies only a linear
representation frame.

## 6. Exact p=13 dimensions

At `p=13`, each central block has dimension `13^2=169`.  Therefore

```text
neutral dimension =169,
charged dimension =12*169=2028,
total dimension   =169+2028=2197.                    (23)
```

For each of the twelve nonzero central characters, `(10)` identifies the
thirteen THM-3240 affine copies with the thirteen multiplicity copies inside
the regular nonvertical flag representation.  The required pointed-current
test is one `13` by `13` determinant.

The full THM-3243 exceptional boundary has one additional vertical affine
orbit.  Its `kappa`-charged block is one extra copy of `pi_kappa`, so the
full boundary has multiplicity `14`, not `13`, at each nonzero charge.
The theorem deliberately targets only the invariant nonvertical `H_13`
stratum.

## 7. Boundaries: positivity, Singer mixing, and physical typing

The transform `(10)` is a cyclotomic linear combination of point masses.
Although its normalized form is unitary, it is not a positive kernel, a map
of sets, a Boolean packet map, or a Markov intertwiner.  It cannot by itself
transport a positive physical current or a common-support chart.  The
positive-cone problem survives intact.

In fact every nonzero nonnegative point function has positive total mass,
so its central-neutral projection is nonzero.  The charged bridge therefore
acts on a centred/signed component of any such packet and cannot carry the
whole nonnegative packet.  This is independent of the determinant gate in
Section 5: even a cyclic charged part does not remove the neutral mass.

The nonvertical flag stratum and the central polarization are `H_p`-invariant
but not Singer-invariant.  THM-3243's Singer cycle moves all `p+1`
directions transitively, hence mixes the vertical and nonvertical strata;
THM-3234 also shows that it does not normalize the chosen Heisenberg group.
Therefore `(10)` is not Singer- or `AGL_2`-equivariant and does not identify
an owner cycle.

Finally, the source coordinates in THM-3240 already require a labelled
target axis and Bezout section, while the target coordinates in THM-3243
require a contact basis, origin, and affine slope chart.  The explicit map
is canonical only relative to all those choices.  No canonical endpoint
current is proved to have nonzero determinant `(21)`, and no target-side
physical pointed current is constructed.  This theorem proves no endpoint
survival, no all-`91`-unit preservation, no scalar-row exclusion, and no
`LRC(14)` decrement.

## 8. Exact companion

Run from the repository root:

```bash
python3 04-computation/heisenberg_charged_blowup_address_intertwiner_thm3250.py
python3 -O 04-computation/heisenberg_charged_blowup_address_intertwiner_thm3250.py
```

The companion uses exact modular, cyclotomic-exponent, and rational linear
algebra only.  It checks the two actions independently before testing the
intertwiner, records the complete neutral multiplicity defect, and verifies
the Weyl/cyclicity controls without floating point or optimization-sensitive
assertions.  Ordinary and optimized runs must byte-match the pinned
transcript.

QED.
