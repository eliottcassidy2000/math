# The missing phase is real, and the full endpoint spectrum is still a cut

**Status: research reflection for proved THM-3518.**  The theorem and exact
companion are the proof sources.  This file records what the second audit
adds to THM-3515 and why it does not move the ancestry/current boundary.

## Outcome

The all-role candidate is mathematically sound.  Its stored pointwise census
has one transcript-only defect: the `product_zeros_by_chart` tuple contains
71 written zeros although the code and reconstructed bank contain all 72.
The candidate files remain untouched, and the independent output records the
complete tuple.  During the audit, THM-3515 arrived from a concurrent
independent route and promoted its endpoint weighted-tree closure.
The second route did not duplicate that theorem mechanically.  It isolated
the normalization and phase identities on which the bucket formula depends,
then checked the cohomological type on every chart instance.

The decisive equality is

```text
sum_tau g_C(a+tau)g_D(b+tau)zeta^(-k tau)
 =zeta^(ka) sum_s g_C(s)g_D(s+b-a)zeta^(-ks).        (1)
```

The left-sheet phase in `(1)` has positive sign.  A direct tau-slice
contraction, formed before any guard kernel, agrees entrywise with the
canonical-kernel contraction on all five classes and all 117 types.  The
complete reconstructed tensor regenerates the candidate semantic digest,
not only its four nonvanishing censuses.

## What the second audit adds

THM-3515 already proves that all five class rows occupy all 52 active
addresses, have rank five, and make every bridge/tree/product factor nonzero
in all four representations.  THM-3518 adds four independent certificates.

1. **Phase certificate.**  Equation `(1)` is checked at all
   `39^2*2=3042` atom-pair/frequency instances.  Reversing the target,
   tau, or drift-DFT sign changes an explicit coefficient.
2. **Normalization certificate.**  The inverse is genuinely three
   dimensional.  Replacing `13^-3` by `13^-2` changes the first active
   coefficient.
3. **Rank-minor certificate.**  The first five canonical address columns
   already have nonzero determinant; every drift and Fourier rank-four
   statement has a displayed minor certificate.
4. **Cut-cycle certificate.**  Six explicit fundamental cycles pair to zero
   with every one of the 56,592 tested chart cochains, while independent
   sixteen-tree, block-Laplacian, and full-graph determinants agree.

These are proof-quality robustness upgrades.  They do not add a sixth role
class, a new address, or a physical relation between the two endpoint legs.

## The tree/cycle paradox is the mechanism

The private-support graph has eight vertices, thirteen edges, incidence rank
seven, and cycle rank six.  Its two tetrahedra are joined to a leaf by one
forced bridge.  At every address or transformed mode the edge response is

```text
w_(uv)=p_u-p_v.                                      (2)
```

Thus `w` is a cut, or equivalently a graph coboundary.  Every cycle integral
telescopes to zero.  The forced bridge is in no cycle at all.

At the same time every spanning tree must contain that bridge and one tree
from each tetrahedron, so

```text
det L_red=(H-q5) Tau_left(w) Tau_right(w) != 0.       (3)
```

Equations `(2)` and `(3)` are compatible because the matrix-tree polynomial
is not a function on absolute `H^1`.  It sees the chosen cut representative.
This is the concrete U_full realization of the mechanism abstracted in
THM-3482 and corrected in MISTAKE-409.

## Concept board after the audit

| object | exact invariant | operation | lost coordinate | cheapest next test |
|---|---|---|---|---|
| five-class endpoint tensor | support `5*52`, rank five | guard restoration and inverse characters | common base and root | recover the frozen bridge after one lawful joint restriction |
| guard-sheet torsor | `a -> a+tau`, phase `(1)` | translation/Fourier transform | physical root identity | construct one common equivariant gauge inside a Boolean stalk |
| two-K4 tree graph | 72 charts, all factors nonzero | matrix-tree determinant | chronology and ancestry | evaluate only after bridge recovery |
| graph cut space | incidence rank seven | vertex-potential gradient | absolute flux | attach a non-gradient holonomy sidecar lawfully |
| six-dimensional cycle space | all pairings zero | cycle quotient | endpoint magnitude | test a supplied ancestry/phase cochain, not the current endpoint cut |

The board separates two questions that had shared the word “spectral.”
Weighted-tree spectral closure is finished in this endpoint bank.  Absolute
cycle flux has not begun there.

## Connection contract for the next transplant

```text
source:
  the 26 active left atoms and 26 active right atoms before marginalization;

target:
  one THM-2471 Boolean ancestry fibre with a common base and root gauge;

map needed:
  a character-independent relation R subset Omega_L x Omega_R;

preserved:
  both endpoint factors, guard covariance, chamber/drift labels, role class;

destroyed by the current Cartesian product:
  common base, physical root, owner branch, word/source/arrival sheet,
  horizon and chronology;

sidecar:
  one common torsor gauge plus those ancestry fields before Fourier summation;

first decisive test:
  the restricted inverse transform must reproduce the frozen H-q5 bridge.
```

The two tetrahedral factors should remain downstream of this test.  Their
nonvanishing is now so overdetermined that testing them on an unlawful pairing
would add confidence to the wrong object.

## Tournament boundary

The four chamber pairs form `F_2^2`, and each private block is a `K4`, but no
intrinsic pairwise comparison is present.  Owner-order edge orientation is a
sign gauge for the incidence calculation, not a tournament observable.  The
faithful carriers are the address tensor, its Walsh/Fourier transforms, the
cut space, and the cycle space.

## Honest frontier

Direct progress is a second independent, phase-sensitive proof certificate
for the complete endpoint spectrum.  The main LRC obstruction is unchanged:
there is still no lawful character-independent common-stalk relation before
endpoint marginalization.  Consequently there is no physical current,
nonzero absolute `H^1`, LRC `7 x 13` bispectrum, scalar-row exclusion, or
LRC(14) conclusion.
