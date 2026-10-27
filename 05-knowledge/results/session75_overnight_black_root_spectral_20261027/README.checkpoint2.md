# Session75 checkpoint 2 — core targets, source invoices and native pairing debt

**Audited Phase1 results; overnight work continues.** This additive checkpoint
preserves every [checkpoint1](README.md) native byte and its spectral audit.
Source pin: `96e35642110104cc547516a63f37eebcf930f707`.

## Black-core geodesics: retain the target and named starts

**PROVED, independently audited:** for regular tournaments on `2m+1>=7`
vertices, black native distance to a singleton plus one strong core is exactly
`m`. Every balanced even-order tournament at order `2m>=4` has exact distance
`m−1` (structural proof for orders at least8, complete finite proofs at4 and6). The construction retains actual Hamilton words, failed
anti-automorphism equations and two distinct named starts.

A different all-colour theorem gives exact distance `lambda(T)` to the entire
nonstrong locus, where `lambda` is the minimum directed-cut size. This is not
an exact black core-attachment formula. An order-ten witness has cut distance
`1` but core distance `2`; repeated starts can leave source–R5–sink instead
of one strong core. The pure-blue small atoms remain excluded.

Read [proof and literal routes](session75-black-core/BLACK_CORE_GEODESICS.md)
and [independent adversarial audit](session75-black-audit/ADVERSARIAL_REPORT.md).
The regular statement is reconstructed inherited work, not a novelty claim.

## Native ROOT: finite norm is not a well-defined operator

**PROVED, independently audited:** for every cutoff `N>=2`, the prescribed
source-weighted quotient has a positive-mass source mapping outside its mass
support; native pullback and forcing difference therefore do not descend to
that quotient. A finite seminorm alone does not repair the domain.

First-hit transport of the prescribed charge has an explicit initial-energy
invoice. Its factor satisfies `B_s^(-2)<=(1+3d)^(2/9)` on a first-hit word;
nonrepeating exterior prefixes and visit-multiplicity records provide the
sharper scoped controls. This is not a uniform source-summed bound.

**Convention firewall:** zero forcing in `e'=epsilon*b*e+F` differs from zero
forcing in `e'=b*e+f`. Separate source-labelled boundary transports also differ
from one common ROOT boundary value. The C3 map preserves labels/pair ownership
but multiplies both energy and debt by three; its index-three lattice and
marked observer do not pay the source invoice.

Read [native ROOT report](session75-root-transport/REPORT.md) and
[independent audit and qualifications](session75-root-audit/REPORT.md), plus
[final reviewed-byte acceptance](session75-root-c2-acceptance/ACCEPTANCE.md).

## Dyadic debt: the averaged ledger is not a realized pairing

**FINITE-EXACT, independently reconstructed:** the literal `N=512` first-arrival
family has 768 sources, source energy `1244`, terminal energy `1313`, whole-level
debt `79`, positive parent-bin interactions `293`, and positive labelled-pair
interactions `433`. This authenticates that exact instance, not the provenance
of an unlocated older display.

**PROVED:** three-pairing mean identities do not guarantee zero debt in an
actual pairing. Native three-source witnesses retain an empty bin and disprove
that shortcut. A charge-dependent minimum-matching construction does give an
alternative compatible balanced tree with debt at most `(1−1/n)P`; its bound
uses the target terminal energy `P` and does not control prescribed low-bit
chronology. Generic fixed-tree debt/source-norm ratios are unbounded; those
generic adversaries are not asserted to be full native source families.

The revision2 receipt fixes a source-record serialization overwrite and checks
the final marked-hostile records `799/1065` exactly. Original and corrected
manifest lineage is preserved. Read the [dyadic report](session75-dyadic/README.md)
and [independent audit](session75-dyadic-audit/AUDIT_REPORT.md).

## What remains open

No all-source ROOT arrival, Collatz convergence, universal native debt bound,
all-order black connectivity or LRC(14) closure is claimed. Bounds for an
alternative pairing tree are not silently transferred to the fixed native tree.
The next overnight rounds target these precise unpaid conditions, not merely
larger finite scans.

## Reproduction

See [verification scope](VERIFICATION.md) and the explicit
[curation scope](CURATION.checkpoint2.md). The newest frozen manifest retains
checkpoint1 and all approved current normal/optimized jobs:

```sh
python3 -I -B verify_packet.py --replay
python3 -I -B -O verify_packet.py --replay
python3 -I -B verify_provenance.py
python3 -I -B test_packet_corruption.py .
python3 -I -B test_provenance_corruption.py .
```
