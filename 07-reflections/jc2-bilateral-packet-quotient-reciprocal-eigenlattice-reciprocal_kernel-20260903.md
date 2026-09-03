---
title: "JC2 bilateral packet quotient and reciprocal eigenlattice"
status: >
  PROVED elementary relative to THM-4368, THM-4369, and THM-4375;
  VERIFIED-EXACT primary companion, INDEPENDENTLY VERIFIED-EXACT recurrence
  referee, and independent scratch-audit signal.
  JC(2), DC(2), bracket entry, and seam control remain OPEN.
source: reciprocal_kernel / JC2 continuation session, 2026-09-03
artifacts:
  - 01-canon/theorems/THM-4378-bilateral-packet-quotient-reciprocal-eigenlattice.md
  - 04-computation/jc2_bilateral_packet_quotient_reciprocal_eigenlattice_thm4378.py
  - 05-knowledge/results/jc2_bilateral_packet_quotient_reciprocal_eigenlattice_thm4378.out
  - 04-computation/jc2_bilateral_packet_quotient_reciprocal_eigenlattice_independent_referee_thm4378.py
  - 05-knowledge/results/jc2_bilateral_packet_quotient_reciprocal_eigenlattice_independent_referee_thm4378.out
script_sha256: b93dec50de6e75b571c3f1b1cdc1649472504c19b62dfb0c09baf88c7ddaa51b
output_sha256: c23b9afc7682476e1ce0f1e8d9fb1e5e9939395821811871df971386a8d58b6a
independent_referee_script_sha256: 70801c2c943c48258dc6b57781edc4f33248c36b09000662eb3387d93346c31e
independent_referee_output_sha256: 60127767ac254b6069d45e80acdb0466c6927dfd9074a88619f6dc5aca0ac73f
hash_basis: raw LF bytes
related:
  - THM-4368-diagonal-boundary-valuation-triangular-address-and-simplex-stream-rank
  - THM-4369-source-packet-pascal-circuit-kernel-and-boundary-basis
  - THM-4375-bilateral-source-cone-reciprocal-orbits-and-fibre-imbalance
---

# The reciprocal band collapses to a width-one spine

## 1. Inheritance and board

The nearest mechanisms were THM-4368's two-boundary trace, THM-4369's local
Pascal kernel, and THM-4375's bilateral path-power band. The first question
was whether coordinate swap acts on the packet quotient. It does not: the
trace normalization carries a checkerboard sign.

The live board was

```text
bilateral band | Pascal circuit | central spine | checkerboard action
eigenlattice | Smith saturation | norm/difference | source fibre.       (1)
```

## 2. The useful compression

After writing `u=rho+i`, `v=rho+j`, every bilateral trace has the common
factor

```text
q^(rho-1),                  q=z(1-z),                    (2)
```

and residual polynomial `z^i(1-z)^j` of degree at most the excess prefix
depth `D`. Once `ell>=3`, the band contains the alternating central spine

```text
1,z,q,zq,q^2,zq^2,... .                                  (3)
```

Its leading coefficients are units degree by degree, so `(3)` is already a
unimodular basis of `Z[z]_(<=D)`. The width of the path-power band vanishes
from the quotient: loops and nearest-neighbour arcs suffice integrally.

This is positive structure, but also a route warning. Any additive argument
which only reads the complete diagonal trace cannot gain strength by adding
longer reciprocal arcs. It must constrain the spine or restore an erased
coordinate such as source-fibre identity, bracket data, or a seam label.

## 3. The repair and the two-primary seam

For the THM-4369 circuit

```text
C_(u,v)=e_(u,v)-e_(u,v+1)+e_(u+1,v),                    (4)
```

plain coordinate swap obeys

```text
R_0 C_(u,v)=2e_(v,u)-C_(v,u),                            (5)
```

so it leaves the relation kernel. The lawful action is

```text
iota(e_(u,v))=(-1)^(u-v)e_(v,u),                         (6)
```

which intertwines with `f(z)->f(1-z)`. An independent scratch audit recorded
in broadcast commit `c99e9177b5` found the same repair and independently ruled
out hidden four-torsion.

Every polynomial has a unique form `A(q)+zB(q)`. Its invariant part is
integrally `A(q)` and its anti-invariant part is `(2z-1)A(q)`. Each sector is
saturated by itself, but their sum has one Smith factor `2` for every paired
block `(q^k,zq^k)`:

```text
Q/(Q^+ direct_sum Q^-)=(Z/2)^(ceil(D/2)).                (7)
```

The separate norm obstruction is smaller: it is one `Z/2` only when the top
degree is even. Difference maps onto the whole anti-invariant lattice.

THM-4377 supplies the first exact source-sidecar coupling. For a reciprocal
orbit, the two fibre-weighted packet classes have the same nonzero gluing
vector modulo two, so their weighted sum is nonzero exactly when the two
source-fibre sizes have opposite parity. In the postclock range this is the
parity of the growing boundary-jet rank; it vanishes permanently when the
full rank `2d` is reached. This detects the sidecar inside the packet quotient,
but no bracket or seam is yet known to detect that gluing class.

## 4. Smallest hostile and boundary

At `(ell,R)=(3,2)`,

```text
C=e_(1,1)-e_(1,2)+e_(2,1)
```

has zero trace, while its unsigned reflection has constant trace `2`.
Meanwhile `Q=Z{1,z}` and the eigenbasis `Z{1,2z-1}` has index two. This is
simultaneously the first descent failure and the first eigen-gluing class.

At `ell=2`, the band is only the diagonal ray, so `J=id`. There is norm
two-torsion but no eigen-gluing defect. Keeping this endpoint separate avoids
mistaking a missing anti-sector for a split one.

## 5. Next sharp move

The quotient classification shifts the next target away from wider packet
geometry. A genuinely stronger JC2 consumer must couple `(3)` to one of:

- THM-4375's imbalanced source-monomial fibres;
- a bracket coefficient not factoring through the complete trace;
- a boundary jet or seam attachment which detects orientation before the
  two-primary gluing is saturated.

No such coupling is proved here, and no `JC(2)` or `DC(2)` conclusion follows.
