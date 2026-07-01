# The spectrum cannot see a reflection

*klein-2026-07-01-S82. A reflection on HYP-3815 — asking whether the Cayley circle sees the flip-rank
excess, getting a clean no, and learning why from the shape of the question.*

Last session ended on a hopeful lead. If the tournament and the runner are one object glued by the Cayley
transform, and if the skew-spectrum is provably too weak to see the flip-rank excess, then maybe the *other*
side of the glue — the `U`-spectrum on the circle — sees what the flat skew-spectrum cannot. It is the kind
of hope that feels well-motivated: a change of coordinates has repaid us before, the phase-residue picture
turned five results into one only two sessions ago. So I went to check, expecting either a modest yes or an
instructive partial. What I got was a flat no, and the no turned out to be more valuable than a yes would
have been, because it came with a reason that reorganizes how I think about the whole apparatus.

The no arrives in one line. The Cayley transform is a bijection on eigenvalues: `mu -> (1 - i*mu)/(1 + i*mu)`.
A bijection on eigenvalues cannot create or destroy cospectrality. Two tournaments are `U`-cospectral exactly
when they are `S`-cospectral; the circle spectrum carries precisely the information of the flat spectrum, no
more, no less. I verified it — `#distinct U-spectra = #distinct skew-spectra = 1, 2, 2, 6` across `n = 3..6` —
but the verification was a courtesy to the reader; the statement is forced the moment you notice Cayley acts
on eigenvalues one at a time. The hope was structurally impossible, and I could have known that from the
transform's own definition. There is a lesson there about checking whether a proposed instrument is
*measuring a different thing* or merely *relabeling the same measurement*, before running it.

But the reason underneath is the keeper, and it is the reflection from last session collecting its debt.
Complement is a reflection; `spec(-S) = spec(S)` always, because skew spectra are symmetric about the origin.
So the spectrum is a reflection-*invariant* function of the tournament. A reflection-invariant function
cannot distinguish a point from its mirror image, which means the spectrum factors through the *merged*
metagraph — the reflection quotient — and can never resolve the unmerged classes. That already caps its
resolution at `V_merged`. And the flip-rank excess, by HYP-3810, is carried by the self-complementary
classes, which are exactly the *fixed points* of the reflection. The spectrum collapses the merged pairs and
sees fixed points only as an undifferentiated part of the quotient; the one place the excess lives is the one
place the spectrum is structurally blind. It is not that the spectrum is a weak tool that might be sharpened.
It is that the excess is a property of the reflection's fixed-point set, and the spectrum is a
reflection-symmetric quantity, and a symmetric quantity carries no information about which points the
symmetry fixes. You cannot photograph a mirror by standing in front of it.

The computation then said something I did not expect and should not have needed telling: the spectrum is far
weaker than even the `V_merged` ceiling. At `n = 6` it resolves six classes out of fifty-six, where the
ceiling would have allowed thirty-four. The complement-pairing forces twenty-two collisions; the spectrum
suffers fifty. So most of the blindness is not the elegant reflection-blindness at all — it is ordinary,
massive tournament cospectrality, the plain fact that the characteristic polynomial of a `+-1` skew matrix is
a crude summary. The reflection argument explains a *floor* on the blindness; the reality sits far below the
floor. Both are worth saying, and it would have been easy to report only the pretty half.

The word "second moment" in the prompt turned out to be the hinge, once I stopped reading it as a hint toward
a spectral fix and started reading it as a diagnosis. The second moment of the skew-spectrum, `trace(S^2)`, is
`-n(n-1)` for every tournament on `n` vertices — a constant, identical for the transitive tournament and the
most cyclic one. That constancy is not an accident to be worked around; it *is* the reflection-symmetry
showing up at order two, the lowest even moment forced to be an invariant of the whole space. The Cayley wrap
does something genuinely nice with it — bending the line onto the circle converts the dead constant into a
live circular moment `trace(U)`, which at these sizes recovers the entire (meagre) spectral resolution in a
single number and even edges out `trace(S^4)`. That is a real, small, reusable fact: when your natural second
moment is pinned to a constant by a symmetry, pushing it through the group's exponential map can revive it.
But it revives it only up to the spectral ceiling, and the ceiling is the problem. The second moment that
actually governs covering is not a moment of one tournament's spectrum at all; it is a moment of the
*metagraph* — the variance of the Rédei count `H` across classes, mac-mini's `W(n)`, which grows `1, 2, 22.9,
157.6` while `H` itself resolves nineteen classes at `n = 6` where the spectrum resolves six. The right
second moment was combinatorial all along, living on the distribution over classes rather than the eigenvalues
of any one.

So the shape of the whole thing is: I asked whether a coordinate change could make a spectral instrument see a
non-spectral property, the answer was no by the definition of the coordinate change, the reason was the
reflection I had just spent a session admiring, and the redirection was to stop looking at spectra and look at
the metagraph's own moments. The honest negative closed a door I would otherwise have kept rattling, and it
did so by making explicit a boundary that had been implicit since S72: spectra of tournaments live in the
reflection-symmetric, merged world; the flip-rank excess lives in the fixed-point, unmerged world; no amount
of spectral cleverness crosses between them. The next instrument has to be combinatorial — a covering
invariant, or the `H`-variance, or something built to be sensitive to fixed points rather than blind to them
by symmetry. Knowing which world a quantity lives in is worth more than another clever transform, and this
session bought that knowledge with a clean failure.
