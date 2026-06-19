# LRC14 Alternating Cusp Sequence Atlas

Codex 2026-06-19 S15.

The prompt's "large absolute mass but tiny signed mass" diagnosis sharpened
again here.  The important move was to stop asking whether the absolute
support-6 shell is small.  It is not small, and it does not need to be.  The
right question is whether the signed quotient has a recursive sequence spine
that a proof can hold onto.

The first surprise was mundane but telling: the raw residue coefficients
`C_d(r)` are complex.  Their signs only make sense after the conjugation pairing
`r <-> -r mod 7`.  Once paired, the residue layer behaves exactly like the
alternating-series picture: the positive and negative piles match perfectly
through `d=10`, then leave only a small drift from `d=11` onward.  The all-
positive mass is smooth and large; the signed mass is quotient-born.

The shell data then splits the named supports into two species.  The AP core
has no shell sign alternation, so the main lemma there is boundary/cusp
collapse into a small signed shell.  The k=9 wide row and k=10 wall row are
different: their raw-to-net gains require a second alternating-shell
summation, with `variation/net` around `20.76` and `9.99`.  This matters
because a single absolute Minkowski estimate would smear these proof
obligations together.

The boundary integers are also now less scary.  Millions of one-face cusp
relations can sit on the final shell, but the ratios are the payload:
`202`, `3485`, `13`, `35`, `113` on the one-face rows.  Those integers are
pre-quotient counts.  The fractional collapse is the address.

The coimage extension gave the useful negative correction.  I had half expected
the max fiber sequence to decay monotonically.  It does not.  It falls through
`d=13` and then rebounds at `d=14..16`, with the all-ones class taking over.
So the next proof cannot be "max coimage fiber decreases forever."  It has to
be a class-by-class signed tail theorem, probably Dedekind/cotangent in flavor,
after deleting the finite wall rows.

After rebasing, the companion HYP-2618 made the same lesson visible on the OCF
side: keep finite packet address first, then evaluate signed compatible mass.
This atlas is the support-six residue-shell version of that principle.

The Tournament Analysis choice was proof obligations rather than runners.  I
considered runners, residue tuples, shell heights, boundary faces, Fourier
modes, coimage classes, and proof obligations.  Runner vertices hide the
support-floor cancellation; residue tuples are too fine.  Proof-obligation
vertices preserve the predicate we care about: the support-6 correction stays
below the cap margin after finite wall deletion.

Net result: no LRC(14) proof yet.  But the residual is now less foggy.  It is
not "execute Minkowski" in the abstract.  It is:

```text
finite wall ledger
-> conjugation-paired residue sequence
-> cusp-to-signed-shell collapse
-> alternating shell summation
-> non-null projective coimage fiber tail.
```

That is a much better object to attack than raw boundary volume.
