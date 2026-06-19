# LRC14 Support-Six Coimage Fiber Atlas

Codex 2026-06-19 S14.

The prompt said "think coimage", and that was the right lever.  S12 had already
attached the residue address `C_d(n mod 7)` to each exact support-6 term, and
S13 had made the philosophical point that the absolute boundary mass lives
before quotienting.  This session made that quotient finite.

The computation is simple but clarifying.  For a six-speed support, reduce the
speed residues mod `7`:

```text
a_i = e_i mod 7.
```

The coefficient residues remain nonzero:

```text
r_i in F_7^*.
```

Then the relation hyperplane has the finite shadow

```text
a_1 r_1 + ... + a_6 r_6 = 0 mod 7.
```

Summing `C_d(r)` over that fiber gives the leading coimage coefficient `S_d(a)`.
After quotienting by scalar multiplication and coordinate permutation, there
are only `159` speed-residue classes.  That is the new atlas.

The important correction is that `a_i=0` must be allowed.  A speed divisible by
`7` has not disappeared from the Fourier term, because its coefficient residue
can still be nonzero and 7-coprime.  It disappears only from the mod-7 relation.
That is exactly the coimage degeneration: Fourier-live, relation-invisible.

The named rows make the story more concrete.  The AP core class
`(1,2,3,4,5,6)` at `d=7` has `S_d=-0.93653539`, already about a 16x cancellation
against absolute mass.  The resonant `21` support, whose class has one zero
speed residue, shrinks to `S_d=-0.11706692`, about 125x cancellation.  The
surprise is the `k=10` height-one wall support.  In its relevant dimension
`d=9`, class `(0,1,1,1,2,4)` is coimage-null to numerical precision even though
its absolute mass is large.

That changes the proof posture.  A scary low-height wall can be dangerous in
the finite ledger but harmless in the infinite analytic tail.  The final proof
should not treat all low-height support-6 relations as tail mass.  It should
delete and account for the finite wall rows, then bound only the non-null
projective coimage fibers by a signed reciprocal estimate.

The resulting proof stack is now:

```text
THM-538 support floor
-> S12 residue reciprocal identity
-> S13 signed-mass sequence spine
-> S14 finite projective coimage atlas
-> remaining signed reciprocal-tail theorem.
```

LRC(14) is still open.  But the residual is no longer "do the Minkowski count"
in a vague sense.  It is a finite table plus a Dedekind/cotangent-looking
tail bound after low-height wall deletion.  That feels like a much more
mathematical shape to attack next.
