# Audit: one edit of a nonprimitive mechanical repetition

20 October 2026. **Verdict: valid, including crossing-cut edits and full
prime-power content.** This elementary lemma allows the completed primitive
positive one-swap theorem to drop its coprimality hypothesis.

Let `0<p<q` and `d=gcd(p,q)>1`. Write `p=d p0`, `q=d q0`. The zero-intercept
mechanical word is `c=v^d`, where `v` is the mechanical word of shape
`(p0,q0)`. This follows directly from the floor definition and its period
`q0`. Define

    D0=2^q0-3^p0,
    H=sum_(j=0)^(d-1) 2^(q0(d-1-j))*3^(p0 j).

Then

    D=2^q-3^p=H D0,
    B(c)=H B(v).

The second formula follows by composing the block affine map repeatedly.
The first is the difference-of-powers identity.

The integer `H` is greater than `1`. It is odd, because only its final term
is odd, and it is nonzero modulo `3`, because only its first term survives
modulo `3`. Therefore `gcd(H,6)=1`.

For an edited edge inside a chronological cut, exchanging `01` with `10`
changes the carry by exactly `-2^a3^b`, where `a` is the first edge position
and `b` counts the ones strictly after the pair. The reverse exchange changes
it by `+2^a3^b`. Contributions before and after the pair do not otherwise
change. Consequently

    B(w)=H B(v) +/- 2^a3^b,
    gcd(H,B(w))=1.

Since `H` is coprime to the edited carry, including at every prime-power depth,

    g(w)=gcd(|H D0|,B(w))
        =gcd(|D0|,B(w))
        <=|D0|
        <|D|.

The last inequality is strict because `H>1` and `D0!=0`. The latter follows
from unique factorization of positive powers of `2` and `3`.

For an edit crossing the original cut, first rotate the cut away from the
edited edge. A rotation of `v^d` is still the `d`-fold repetition of a rotated
block, so the same argument applies. It also gives the exact identity at the
original cut, not merely properness at the new one: under one left rotation,
with initial bit `e`, carries satisfy

    2 B(rotated w)=3^e B(w)+e D.

Both `H` and `D0` divide `D` and are coprime to `6`. Thus rotation preserves
the gcd with either factor, and with `D`, transferring every claimed gcd
statement back to the original word.

This proof does not use the sign of `D`. It excludes integer fixed cycles of
either sign for the nonprimitive one-swap family. Together with the audited
primitive positive exclusion, it excludes the positive family without a
coprimality assumption, subject to the same stated length and clock conditions.
It does not assert anything about two edits, arbitrary nonprimitive words, or
general Collatz convergence.

The separate finite checker `verify_nonprimitive_repetition.py` constructs
every actual cyclic unequal-adjacent swap of every noncoprime shape
`0<p<q<=100`. It checks the repetition and carry identities, exact local edit
difference, original-cut and rotated gcd statements, and strict properness.
It is corroboration of the symbolic lemma, not the source of its quantifiers.

The check passed in normal and optimized Python with byte-identical receipts:
1907 nonprimitive shapes, 67,672 actual swaps, 1907 crossing-cut swaps, 49,192
positive-clock occurrences, and 18,480 negative-clock occurrences.

Reproduce from the workspace:

    python3 session13-collatz-resultants/refinement/audit/verify_nonprimitive_repetition.py
    python3 -O session13-collatz-resultants/refinement/audit/verify_nonprimitive_repetition.py
