# All-modulus layer transport, divisor ancestry, and the harmonic atlas

**Research reflection / provenance, not a truth source.**  Exact claims are
routed to repaired audit-pending THM-3472, MISTAKE-407, and proved
dependencies THM-3405, THM-3415, and THM-3453.  LRC(14) remains open.

## The audit changed the theorem, not merely its wording

The first THM-3472 candidate said that an odd fixed-zero primitive cover
became a primitive half cover after doubling its owners.  That was false in
the repository's augmented sense:

```text
Q=15, S=(1,2,3,4,5,7): gcd(15,S)=1,
2S=(2,4,6,8,10,14):   gcd(30,2S)=2.
```

The failed coordinate was `2Q`, not `Q`.  Yet the rank proof asks only for an
ordinary literal half cover.  Once that type is kept honest, the even hostile
at `Q=8` becomes informative rather than terminal: `ell -> 2ell+1` is not a
permutation, but it still maps every target sheet into the odd source coset.
The self-opposite fixed owner is empty on that coset and can be deleted.

The repaired transport is therefore

```text
fixed-zero cover on Q
  -- sample on image(ell -> 2ell+1); double owners --> half cover on Q
  -- divisor dilation --------------------------------> half cover on q.
```

For odd `Q`, the first arrow is a conjugacy.  For even `Q`, it is a one-way
cover morphism.  This distinction is exactly why the stronger all-modulus
rank equality survives while augmented primitivity does not.

## The divisor interface is an equation, not an ancestry slogan

With `U=dV`, `g=gcd(q,d)`, `q=gQ`, and `d=gd_0`, THM-3405's two centres give

```text
u(c_epsilon+ell/q)
 =epsilon v/(2Q)+d_0 v ell/Q.
```

This formula names every retained coordinate:

```text
object       = labelled mask packet on the active quotient Q,
operation    = quotient followed by divisor dilation,
grade        = cover cardinality,
unit data    = d_0 in (Z/QZ)^x and the Boolean layer epsilon,
forgotten    = primitive ancestry after taking the minimum,
not kept     = augmented prime-breaker gcd under owner doubling.
```

Divisor dilation still composes multiplicatively, so the family is a
degree-graded monoid rather than a list.  But the corrected objects are
ordinary covers at the transport stage; attaching the word “primitive” would
change the category and invalidate the map.

## A complete all-natural-number harmonic atlas

THM-3453 supplies fifteen literal atoms.  The repaired equality imports all
of them into the global zero-cochain rank.  Their word has period

```text
2^3*3^2*5^2*7*11*13*17*19*23*29*37
 =14,362,718,970,600.
```

Its exact natural and harmonic coefficients are

```text
rank 4:   2/9,
rank 5:   4/45,
rank 6:   25/207,
rank 7:   165741596/1554406815,
rank >7:  717341696/1554406815.
```

This is the precise sense in which these subsets of the natural numbers are
subsets of the harmonic series: periodicity supplies logarithmic
coefficients.  An arbitrary subset has no density or harmonic coefficient
for free.  The old odd word remains a conditioned subatlas rather than a
discarded calculation.

## Boolean derivatives meet factorial face descent

Incoming THM-3466 gives the factorial face identity

```text
B_I(h)=L_n(product_(i in I)(1-partial_i)h).
```

The private-sheet atlas being developed in THM-3473 uses the incidence
derivative

```text
delta_i=I_i product_(j!=i)(1-I_j).
```

These are not the same operator, but they share a useful grammar:

- a Boolean coordinate deletion is expanded by commuting local operators;
- the surviving term lives on a boundary/private stratum;
- the scalar quotient forgets which coordinate supplied that stratum;
- a current theorem needs the labelled boundary response, not only its sum.

This suggests a disciplined D5-map route.  The source should be a labelled
Boolean deletion complex, the target a boundary-current cochain, and the map
must commute with the coboundary.  Static cover incidence alone cannot create
the Keller flux in THM-3466, but its private atoms may provide the missing
basis on which such a map is specified.

## Why the generalized tournament is only a quotient

Pairwise coactivity is symmetric, so the intrinsic directed encoding has
both arcs or no arcs.  It is a tournament with missing and both-way edges in
the user's generalized sense, not an oriented tournament.  The quotient
retains whether two owners ever meet.  It loses sheet multiplicities,
private singleton supports, divisor ancestry, and the even-map kernel.

The `Q=8` correction is a canonical warning: a nonbijective sheet map can
still preserve the target cover predicate.  Pairwise graph isomorphism would
have been unnecessarily strong; image coverage plus an empty-kernel sidecar
is the right invariant.

## Updated concept board

1. **All-modulus re-audit.**  Independently derive the active-gcd formula,
   self-opposite deletion, all fifteen atom strata, and the huge-period CRT
   counts before promotion.
2. **Irredundant versus minimum.**  Complete THM-3473's exact eleven-support
   Boolean atlas.  Rank-four members can still have an irredundant eight-owner
   presentation, proving that grade forgets deletion ancestry.
3. **Current-aware transport.**  Attach mode endpoints and boundary labels to
   the layer map.  Test whether Boolean private atoms lift to nonzero current
   coordinates or die under the same common-centre collapse.
4. **Factorial-face/D5 bridge.**  Seek a cochain map from owner-deletion faces
   to THM-3466's Keller boundary-current block.  The cheapest hostile is a
   cover with positive private atoms but zero transported current.
5. **Spectral closure.**  Rank equality remains upstream of the LRC
   `7 tensor 13` bispectrum.  A bridge must preserve the contraction, not
   merely cover grade or pairwise coactivity.

The central research move was to weaken the map while strengthening the
theorem: from primitive conjugacy on odd moduli to ordinary cover transport
on every modulus.  The lost coordinate became explicit, and the global
conclusion expanded.
