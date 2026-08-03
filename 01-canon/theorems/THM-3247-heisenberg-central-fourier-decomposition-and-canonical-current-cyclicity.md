---
id: THM-3247
title: "Heisenberg central-Fourier decomposition and canonical current cyclicity"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  Over any
  characteristic prime to p splitting field, the standard p^2-point affine
  H_p permutation module is the multiplicity-free sum of p center-blind
  characters and the p-1 p-dimensional charged Schrodinger modules.  For
  every nonzero q in F_13^2, all thirteen scalar-character projections of
  the canonical THM-2790 endpoint current survive in both certified finite
  fields.  Together with THM-2790's twelve charged blocks, the current is a
  cyclic vector of orbit-span dimension 169.  In fact every one of its 169
  two-dimensional Fourier modes is nonzero, so the regular translation
  subgroup alone has orbit rank 169.  This closes an algebraic rank debt on
  the coefficient carrier, not physical high-sheet descent, a same-ancestry
  intertwiner, an owner map, or LRC(14).
source: root/creative-reframes/2026-08-03
audit: >
  The assertion-independent exact companion pins THM-2779, THM-2790,
  THM-3240, and the inherited THM-2790 script/output; exhaustively checks the
  decomposition and orbit-rank controls over exact splitting fields for
  p=3,5,13; rebuilds both certified endpoint-factor fields; replays all
  26,208 charged coordinates; and checks the new 2,184 scalar modes, all
  28,392 two-dimensional Fourier modes, all 4,200 irreducible-block
  projections, and all 168 full-rank H and translation directions per field.
  Normal and optimized runs byte-match the stored transcript.  An
  independent hostile audit rederived the group/action signs, basis versus
  pullback conventions, charged irreducibility and multiplicity-free
  cyclicity criterion, determinant-coordinate chart and section-change
  phases; replayed both exact runs and all pinned hashes; confirmed all
  scalar/charged/block censuses and characteristic-zero inference; and
  verified that the exact-address target semantics and physical intertwiner
  remain outside the theorem.  A second independent all-direction scout
  recomputed all 56,784 two-dimensional modes across the two fields, found
  the per-field rank histogram {169:168}, and supplied the fully mode-dependent
  direction-scaling hostile recorded in the synthesis.
depends_on:
  - THM-2779-bockstein-symplectic-decoder-frame-torsor-and-heisenberg-root-degree-gate
  - THM-2790-universal-depth-two-central-response-and-carry-wall-spectrum
related:
  - THM-2625-canonical-endpoint-current-full-transvection-sector-survival
  - THM-2782-semantic-arm-right-wing-local-unit-and-endpoint-deck-boundary
  - THM-2806-literal-fixed-sheet-central-allocation-scalar-law-and-endpoint-translation-no-go
  - THM-3240-exact-address-heisenberg-clutch-on-carrier-imbalance
  - THM-3243-contact-deformation-blowup-equivariance-and-full-orbit-resolution
  - THM-3250-charged-heisenberg-blowup-address-intertwiner-and-pointed-multiplicity-gate
script: 04-computation/lrc14_heisenberg_central_fourier_current_cyclicity_thm3247.py
output: 05-knowledge/results/lrc14_heisenberg_central_fourier_current_cyclicity_thm3247.out
script_sha256: b33f88d792135c44b0c8c7ddb7982e40f11af905f7856a71ceb12f2f298303c3
output_sha256: fb5f05ad46ee71d56c097feeda18eceb8cdc8848a9876a8395e553a8f0ad4fd6
hash_basis: LF-normalized bytes
---

# THM-3247 -- Heisenberg central-Fourier decomposition and canonical current cyclicity

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-2779 supplies the normalized Heisenberg endpoint carrier, and THM-2790
proves that every nontrivial central Fourier coordinate of the canonical
endpoint current survives.  The apparently missing thirteen dimensions are
not more charged modes.  They form the center-invariant regular quotient and
split into thirteen scalar characters.  An exact transverse transform shows
that every one of those characters survives too.

Thus the canonical current has a nonzero projection to every irreducible
summand of the `169`-state permutation carrier.  Its Heisenberg orbit spans
the full carrier module.  This is a representation-rank theorem on the
coefficient endpoint plane; it does not identify that plane with a physical
same-ancestry current.

## 1. Standard affine carrier

Let `p` be prime and let `K` be a field whose characteristic does not divide
`p` and which contains a primitive `p`th root `zeta`.  Use

```text
H_p=F_p^3,
(x,y,c)(x',y',c')=(x+x',y+y',c+c'-y x').             (1)
```

On `X=F_p^2`, with coordinates `(v,w)`, let

```text
(x,y,c):(v,w) |-> (v+x,w+c-yv).                      (2)
```

This is the THM-2779 convention and the fibre action of THM-3240 after its
label and section choices.  Let `M=K[X]` have basis `e_(v,w)`, permuted by
`(2)`.  For `lambda in F_p`, put

```text
E_v^lambda=sum_(w in F_p) zeta^(-lambda w)e_(v,w),
V_lambda=span_K{E_v^lambda:v in F_p}.                 (3)
```

Reindexing `w'=w+c-yv` gives the exact action law

```text
(x,y,c)E_v^lambda
  =zeta^[lambda(c-yv)] E_(v+x)^lambda.                (4)
```

Consequently each `V_lambda` is invariant and

```text
M=direct_sum_(lambda in F_p) V_lambda.                (5)
```

This is central Fourier decomposition: the center `(0,0,c)` acts on
`V_lambda` by `zeta^(lambda c)`.

The convention here is **basis permutation**, `h.e_R=e_(hR)`.  Acting on
coefficient functions by pullback is the contragredient convention and
replaces `(lambda,a)` by `(-lambda,-a)`; every support, rank and cyclicity
statement below is unchanged.

## 2. Multiplicity-free irreducibles

When `lambda=0`, the `y`-line and center act trivially and the `x`-line acts
regularly.  Define

```text
U_a=sum_(v in F_p) zeta^(-av)E_v^0.                   (6)
```

Then

```text
(x,y,c)U_a=zeta^(ax)U_a,                              (7)
```

so

```text
V_0=direct_sum_(a in F_p) K U_a.                      (8)
```

For `lambda!=0`, the `y`-line has the `p` distinct weights

```text
y |-> zeta^(-lambda yv),              v in F_p.       (9)
```

Because `char(K)` does not divide `p` and `zeta` lies in `K`, the exact
`y`-line Fourier projectors onto these weight spaces exist.  Any invariant
subspace is therefore the span of some of the one-dimensional weight spaces.
The `x`-translation in `(4)` is transitive on those spaces, so a nonzero
invariant subspace contains all of them.  Hence `V_lambda` is irreducible of
dimension `p`.  Distinct nonzero `lambda` have distinct central characters,
and the characters `(7)` are distinct.  Therefore

```text
M = (direct_sum_(a in F_p) K U_a)
    direct_sum (direct_sum_(lambda in F_p^*) V_lambda),

dimensions: p*1 +(p-1)*p = p^2,                       (10)
```

is multiplicity-free.

By semisimplicity, a vector `m in M` is cyclic exactly when its projection to
every summand in `(10)` is nonzero.  The central Fourier projectors isolate
the `V_lambda`; inside `V_0`, the `x`-line Fourier projectors isolate the
thirteen `U_a`.  Every nonzero vector in an irreducible summand generates that
summand.  Thus its orbit span has dimension

```text
#{active scalar characters}
 +p*#{active nonzero central characters}.             (11)
```

Three exact controls make the boundary sharp:

```text
constant vector:             orbit rank 1,
one charged Fourier fibre:   orbit rank p,
one point delta:             orbit rank p^2.           (12)
```

In particular every one-dimensional representation of `H_p` kills the
center, because the center is the commutator subgroup.  A response state
with a nontrivial central character has dimension at least `p`: translating
one nonzero `y`-weight vector by the `x`-line produces `p` distinct weights.

## 3. The thirteen missing projections

Now set `p=13`.  For every nonzero direction

```text
q in F_13^2,                                            (13)
```

THM-2790 defines the unnormalized canonical endpoint current

```text
J_q(R)=P_(R+q)Q_R.                                     (14)
```

Choose any `t` with `det(q,t)=1` and write uniquely

```text
R=wq+vt.                                               (15)
```

The value `v=det(q,R)` is canonical.  The choice of `t` only chooses origins
on the central `q`-cycles.  Put

```text
g_q(v)=sum_(w in F_13) J_q(wq+vt),
Ghat_q(a)=sum_(v in F_13) g_q(v)zeta^(-av).            (16)
```

Up to the harmless Fourier sign convention, `Ghat_q(a)` is exactly the
projection of `J_q` to the scalar line `K U_a`.  The companion reconstructs
the two certified THM-2625 finite-field specializations and checks

```text
Ghat_q(a)!=0
for all q!=0 and all a in F_13.                        (17)
```

In each field this is

```text
168*13=2,184/2,184 nonzero scalar-character modes.     (18)
```

THM-2790 already proves, more strongly, that every coordinate of every
charged projection is nonzero:

```text
168 directions *13 determinant fibres *12 characters
 =26,208/26,208 nonzero charged coordinates.           (19)
```

Therefore every direction has

```text
13 scalar blocks +12 charged blocks =25/25
```

active irreducible summands.  Across all directions the block census is

```text
168*25=4,200/4,200.                                    (20)
```

The same computation gives a stronger abelian statement.  Put

```text
Jtilde_q(a,b)
 =sum_(v,w in F_13) J_q(wq+vt) zeta^(-av-bw).         (20a)
```

For every `q!=0` and every `(a,b) in F_13^2`, both certified fields satisfy

```text
Jtilde_q(a,b)!=0,
168*169=28,392/28,392 nonzero modes per field.         (20b)
```

The subgroup

```text
T={(x,0,c):x,c in F_13}=F_13^2                       (20c)
```

acts regularly by translations on `(v,w)`.  Its orbit-span rank is exactly
the number of nonzero characters in `(20a)`.  Hence every `J_q` is already a
cyclic vector for `T`, with rank `169`; no shear generator is needed for the
rank conclusion.

The two finite-field specializations are certified cyclotomic ring maps.
Nonzero images in either one already prove characteristic-zero
nonvanishing; agreement of both is the maintained hostile control.  Combining
`(11)`, `(17)`, and `(19)` gives the same conclusion through the irreducible
`H_13` ledger:

```text
dim_K span{h.J_q:h in H_13}=13+12*13=169               (21)
```

for every `q!=0`.  Thus every canonical `J_q` is a cyclic vector for its
normalized affine Heisenberg permutation module.

Changing `t` in `(15)` changes charged phases but not their nonvanishing,
and leaves the determinant coordinate and scalar transform `(16)` unchanged.
Hence cyclicity is independent of the transverse section.

## 4. What rank debt is now closed

The carrier response-state ledger is exact:

```text
center-blind scalar line:                   1,
one charged central-character channel:    13,
all charged channels:                    156,
full affine carrier:                     169.           (22)
```

THM-2790 had already filled the `156`-dimensional charged sector.  Equations
`(17)--(21)` prove that its complement is not a hidden degeneracy: the
canonical current also fills all thirteen center-blind characters.  No
additional coefficient support or representation rank is needed on this
endpoint carrier.  Equation `(20b)` strengthens this: even the regular
translation subgroup sees every character.

This also turns Heisenberg orbit evaluation into a finite exact response
compiler.  After one central/transverse Fourier transform, each scalar block
evolves by one scalar character and each charged block by the explicit
`13`-state law `(4)`.  The minimal state size is one for a center-blind lane,
thirteen for one charged lane, `156` for all charged lanes, and `169` for a
generic full-carrier current.  This is an exact closed action formula, not a
claim that arbitrary input sequences become constant-width recurrences.

## 5. Hostiles and the remaining physical bridge

The following implications remain false.

1. Separate endpoint-factor edge gates do not imply current cyclicity.
   THM-2790's exact one-cycle hostile passes the separate central-edge gates
   while its product has central-translation orbit rank one and no charged
   central modes.  Independently, the constant-vector control `(12)` has full
   `H_p`-orbit rank one.
2. A central response does not imply physical quotient descent.  THM-2782's
   semantic segment retains a load-bearing high address digit and fails that
   descent.
3. A fixed-sheet allocation does not automatically carry the affine action.
   THM-2806 gives the raw-flat/allocation obstruction.
4. An exact-address Heisenberg carrier is not automatically a physical
   current carrier.  THM-3240 requires a labelled target axis and Bezout
   section and explicitly leaves the current intertwiner open.

THM-3241 supplies a contact-deformation copy of the same abstract `169`-set,
but its derivative/basis choices and Singer multiplier do not preserve the
old endpoint-target semantics; moreover Singer conjugation enlarges `H_13`
to the affine group.  That bridge belongs to the synthesis ledger, not to the
dependency graph of this rank theorem.

The precise open connection contract is therefore:

```text
source:
  canonical coefficient current J_q on the endpoint plane;
target:
  a physical same-ancestry or exact-address endpoint-current state;
known map:
  abstract coordinate identification (v,w) and H_13 action;
preserved:
  central character, orbit rank, determinant fibre;
destroyed/not yet typed:
  ancestry, high-sheet address, exact-address Bezout/physical-section
  compatibility, physical sign;
needed sidecar:
  a lawful physical current plus an H_13-equivariant descent/intertwiner;
cheapest decisive test:
  one exact full-orbit covariance audit with ancestry and section changes
  retained.                                                   (23)
```

Until `(23)` is supplied, `(21)` proves no Boolean row exclusion, owner or
root map, physical high-sheet descent, or instance of LRC(14).

QED.
