# The quadratic character keeps deciding how the involution acts - four sightings

**opus-2026-07-20-S403.** Prompted by an owner probe (Fibonacci/Pisano <-> 1001 <-> Erdos 592
<-> the JC triple collision). Most of that chain does not hold - see HYP-8230 for the triage
- but working it surfaced a motif the repo has now hit **four independent times**, in four
different subfields. Worth naming before a fifth sighting gets mistaken for a discovery.

## The motif

> A **quadratic character** decides whether a natural **involution** acts freely,
> reverses orientation, or splits an object into `1 + 2`.

| # | Setting | The character | What it decides |
|---|---|---|---|
| 1 | Tournaments (repo canon) | Redei sign = discriminant character | orientation of the class |
| 2 | Paley `T_p`, THM-1380 SS3 | `p mod 4` (is `-1` a QR?) | complement is an automorphism vs **anti**-automorphism; in the LRC form, `D_v` is `s`-invariant vs `s`-**free** by the parity of `v` |
| 3 | JC counterexample, THM-1375 (II,III) | fibre discriminant a non-square | Galois-ness of the degree-3 cover; `d=3` is the **unique** degree detectable by a single quadratic character |
| 4 | Fibonacci (new here) | `5 mod p` (is `5` a QR?) | `pi(p) | p-1` when `p = +-1 (mod 5)`; `pi(p) | 2(p+1)` when `p = +-2 (mod 5)` |

Verified for #4 at `p = 2,5,7,11,13`: `pi(11)=10` divides `10` (split); `pi(7)=16` divides
`16`, `pi(13)=28` divides `28`, `pi(2)=3` divides `6` (inert).

The common shape: a rank-2 object (a quadratic field, a `+-1` eigen-split, a two-element
orbit) carries a Z/2, and a Legendre symbol says whether that Z/2 is *internal* (fixed
points, split, orientation-preserving) or *external* (free, inert, orientation-reversing).
THM-1380 is the sharpest instance: there the character's answer was the difference between
Brouwer's hypothesis and Borsuk-Ulam's.

## What this is NOT

**It is not evidence of a unification, and it must not be cited as such.** The four
settings share no object; what they share is that Z/2 has exactly two behaviours and
quadratic reciprocity is the standard instrument for choosing between them. That is close
to a tautology about rank-2 phenomena.

The discipline is the one from THM-1380 SS1: *a structure present for trivial reasons is not
a discriminant.* Worked example from this session - it is tempting to note that `1001 =
7*11*13` splits mod 5 as **one split prime (11) + two inert (7,13)**, rhyming with the JC
triple fibre's **1 sigma-fixed + 2 free**. Reject it: each prime is split-or-inert about
half the time, so a `1+2` pattern among three primes is the **generic** outcome, arising
~3/4 of the time by chance. The rhyme carries no information.

## The one operational use

Treat the motif as a **cheap first question**, not a result: when a new involution appears,
ask which character controls it *before* investing in machinery. In S401-S402 that question
would have saved two sessions - the Borsuk-Ulam route died on `ind = 1`, and the parity
dichotomy that survived is exactly sighting #2. Cost of asking: one Legendre symbol.

**Related:** THM-1380, THM-1385, THM-1375, HYP-8230.
