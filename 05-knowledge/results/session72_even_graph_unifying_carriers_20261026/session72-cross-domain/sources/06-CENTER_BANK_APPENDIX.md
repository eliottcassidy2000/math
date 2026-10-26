# Two exact center-bank limitations

**Status: PROVED elementary lemmas.** Independent appendix requested by the
root audit; it is not part of the producer's original sealed proof packet.
The accompanying finite controls illustrate, but do not establish, the
universal statements. No undecidability theorem or root-27 nonhalting theorem
is claimed.

## Definitions and domains

For a nonempty balanced word `P`, let `L=|P|`, `A=3^L`, `B=2^L`,
`C=C_P`, `D=A-B`, and

    F_P(z)=(Az+C)/B,   G_P(t)=F_P^{-1}(t)=(Bt-C)/A,
    r_P=-C/D.

A balanced word has nonnegative prefix credit increments, increment `+2`
for bit 1 and `-1` for bit 0, and total increment zero. In particular `L>=3`,
`C,D` are positive odd integers, and `r_P` is a negative rational with odd
denominator. All native execution statements below enforce the guard before
the update: at `q=1` and even parity, no macro is executed.

## A. No nontrivial finite bank closes under inverse blocks

**Lemma A.** Let `Bcal` be a nonempty collection of nonempty balanced words.
Let `S` be a finite nonempty subset of the real line. If

    G_P(S) subseteq S  for every P in Bcal,

then `S={r}` is a singleton and `r=r_P` for every `P` in `Bcal`.
Conversely such a common-fixed-center singleton is closed.

**Proof.** Each `G_P` is injective. Hence finite cardinality upgrades
`G_P(S) subseteq S` to `G_P(S)=S`. For any one `P` in the nonempty collection,

    diameter(S)=diameter(G_P(S))=(2/3)^L diameter(S).

The finite diameter is therefore zero, as `0<(2/3)^L<1`. Thus `S={r}`.
Closure now reads `G_P(r)=r`, whose unique solution is `r=-C/(A-B)=r_P`.
The converse follows by substitution. This proof works for finite real sets;
rationality is not needed. Nonemptiness of the map collection is essential
for the singleton conclusion.

**Consequences.** A finite nonempty bank closed under every nonempty balanced
inverse block does not exist: `r_100=-9/19` and
`r_101000=-351/665` differ. For a restricted block family, closure exists
exactly when all words have a common center, equivalently when all are powers
of one primitive word (by the producer's finite-word theorem).

This is a finite *set* obstruction under total inverse-map closure. It does
not exclude an infinite bank, a context-restricted transition system, or a
richer state-dependent rank. It also does not assert that all inverse-block
maps are actually traversed on the marked root-7 trajectory.

## B. The full center bank is exactly the legal-lifetime observable

Write `S_all={r_P : P is a nonempty balanced word}`. For a positive integer
source `z` at initial credit `q=1`, let `h(z)` be the number of legally
executed macros before the first blocked `q=1,c=0` state; set `h(z)=infinity`
if no such state occurs. The blocked symbol is not included in `h(z)`.

**Lemma B.** With the supremum taken in the extended nonnegative integers,

    sup_(r in S_all) v_2(z-r) = h(z).

If `h(z)<infinity`, the supremum is attained. For positive `z`, every
individual depth `v_2(z-r_P)` is finite because `r_P<0`.

### Preliminary cylinder fact on rational centers

Parity and the recurrence `T(t)=(3t+c)/2` extend to rationals of odd
denominator: `c` is the numerator modulo 2 after reduction. Choosing that
parity makes the next denominator odd again. The usual induction proves
that a word of length `n` is the first `n` parities iff

    3^n t + C_word = 0 modulo 2^n

in the odd-denominator rationals, i.e. iff its 2-adic valuation is at least
`n`. The same induction proves that two such sources have the same first
`n` parity bits iff their difference is divisible by `2^n`.

Since `F_P(r_P)=r_P`, and `3^L r_P+C_P=2^L r_P`, the first `L` bits of
`r_P` are exactly `P`. Thus its itinerary is `P P P ...`. At initial
`q=1`, every period is legal and returns credit to 1. These rational
periodic centers never reach a blocked symbol. This is a literal rational
parity/credit statement, not an assertion of a positive integer cycle.

### Lower bound

Let `w` be a legally executed prefix of length `n>=1`, ending with credit
`q`. Then

    P=w 0^(q-1)

is nonempty and balanced: the prefix was legal, so its credits remain at
least 1, and appending `q-1` zeros reduces its terminal credit to 1 without
going below 1. Therefore `r_P` has a periodic itinerary whose first `n`
bits equal those of `z`, and `v_2(z-r_P)>=n`.

The appended zeros are a *formal balanced completion*; they are not claimed
to be the actual next bits of `z`. Only the executed prefix is matched.
Consequently every finite legal prefix supplies that much center precision.
For an infinite legal trajectory, the supremum is infinite.

### Upper bound and attainment for finite lifetime

Suppose the first blocked symbol is at index `h<infinity`. Any center whose
depth were at least `h+1` would share all actual parity bits through that
blocked symbol. Starting from the same credit 1, it would therefore reach
credit 1 with bit 0 at the same index, contradicting the never-blocked
periodic-center property. Thus every depth is at most `h`.

If `h>0`, the whole legally executed prefix `w` ends at credit 1, so it is
itself a nonempty balanced word. Its center matches the first `h` bits;
the upper bound forces `v_2(z-r_w)=h`. Equivalently, its next bit is the
initial 1 of `w`, whereas the actual blocked bit is 0.

If `h=0`, the integer source `z` is even. Every balanced center is odd in
the 2-adic sense, so every depth is exactly 0. For example `P=100` attains
that value. This handles the empty-prefix boundary without introducing an
empty balanced block or an undefined center.

The two bounds prove the identity and attainment assertion. In particular,
a finite lifetime from `q=1` is a multiple of three, because its whole legal
prefix is balanced.

## What the identity does and does not provide

The complete infinite center bank exactly recodes the original legal
lifetime. Showing unbounded precision in that bank at `z=7` is equivalent
to showing all-prefix root-27 survival, not a new sufficient criterion that
can be assumed for free. The identity alone supplies neither an independent
bound on center precision nor a finite termination/nontermination test.
It does **not** prove that the lifetime function is uncomputable, that the
specific root-27 question is undecidable, or that no useful rank exists.

For an initial credit `q>1`, the set of balanced-block centers alone is the
wrong unqualified bank for the same formula. The proof uses initial `q=1`
and a terminal single-tree boundary. For example `(z,q)=(4,2)` executes one
zero macro to `(6,1)` and then blocks, so its lifetime is 1, while every
balanced center is odd and the displayed supremum at `z=4` is 0. The
initial-credit hypothesis is therefore necessary for this formulation.
