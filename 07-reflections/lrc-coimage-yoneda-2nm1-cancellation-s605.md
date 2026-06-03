---
source: codex-2026-06-03-S605
status: synthesis / proof-program
tags: [LRC, coimage, Yoneda, resonance, 2n-1, cancellation, HYP-2161, THM-401, THM-406, n14]
---

# Coimage + Yoneda; `2n-1` resonances are the cancellation

The user's sentence has the right compression:

```text
coimage + Yoneda; 2n-1 resonances = the cancellation
```

I read it as a missing bridge between three already-live facts:

1. THM-406: `p_0=sum_j (-1)^j S_j`, so collapse is an all-orders overlap
   cancellation.
2. The resonance-lattice S604 pass: `p_0=sum_{c in L(V)} prod_i kappa(c_i)`,
   so collapse is a Fourier correction by integer relations.
3. THM-401: at the LRC floor, the natural pair-sum/summand-shell modulus is
   `C=2n-1`.

The bridge is: the `C=2n-1` shell ledger is a candidate Yoneda presentation of
the coimage on the floor stratum.

## Coimage tells us what to keep

The coimage is the quotient that forgets everything the observable cannot see.
For LRC, the observable can be presented as the coverage-depth map

```text
N(t)=#{i : ||v_i t|| < delta}.
```

The pushforward of Lebesgue measure by `N` is `{p_k}`.  THM-406 says this is
also the spectral measure of the depth operator, and its factorial moments are
the overlap volumes `S_j`.  The coimage lesson is: do not keep runner identity,
raw phase geometry, or arbitrary pair data unless it factors through the
observable that decides `p_0`.

The resonance-lattice pass gives a second presentation.  The orbit map
`t -> (v_i t)` has annihilator

```text
L(V)={c in Z^m : c.V=0},
```

and Poisson summation writes `p_0` as the free term plus resonance corrections.
This is still the same coimage, now read by characters instead of overlap
moments.

## Yoneda tells us when a probe family is enough

Yoneda is not mystical here.  It is the discipline of asking whether a family of
probes determines the object.  The full functor of probes `Hom(-,Q)` determines
the coimage `Q`; a finite or structured subfamily is useful only if it is
conservative for the predicate we care about.

For LRC, the proposed probe family is:

```text
C=2n-1 summand shells P_a={a,C-a},
unit shell witnesses,
nonunit shell holes,
lift denominators,
endpoint-owner labels,
CRT/determinant states.
```

The conjectural statement is not that these probes recover every phase detail.
They should not.  The statement is that they recover the ground-cell decision:
whether the all-orders cancellation has driven `p_0` to zero while preserving
the floor witness structure.

## Why the cancellation is `2n-1`

The free/independent term in the harmonic formula is positive.  A tight row must
produce a negative resonance correction equal to that free term.  This is
impossible for relation-poor speeds and possible only for relation-rich speeds:
AP rows, additive chains, and composite-shell sporadics.

THM-401 says the first floor-visible resonance quotient is exactly `C=2n-1`.
If a visible unit shell is missed, the clock `a^{-1}/C` gives a witness above
the floor.  If all visible shell probes are silent, the row is already inside a
small arithmetic residual.  At that point the cancellation is no longer a
measure estimate; it is a ledger question over unit shells, nonunit shells, and
lifts.

So "`2n-1` resonances are the cancellation" should be read as:

```text
the all-orders coimage cancellation is first finitely visible at the
THM-401 shell modulus, and the remaining residual is exactly the failure
of those shell probes to separate the coimage.
```

This also explains the sporadics.  Prime `C` leaves mainly unit-shell
transversal flips.  Composite `C` allows nonunit holes: `C=15` at `n=8` and
`C=27` at `n=14` are not nuisances but the ramified part of the same probe
system.

## What this asks for at `n=14`

At `n=14`, `C=27=3^3`.  The target proof program is now clearer:

```text
Show that the C=27 shell/lift/CRT ledger is conservative for p_0 collapse.
```

That has three layers.

1. Unit-shell layer: a missed unit shell gives a `2/27` witness and exits.
2. Nonunit-shell layer: a hole in a gcd-3 shell must lift or be shielded.
3. CRT/determinant layer: if shell probes remain silent, the bounded residue
   automaton should certify emptiness or route to a classified cancellation
   family.

This is the coimage/Yoneda reframing of the n=14 finite ledger.  The shell
ledger is not trying to approximate the overlap sequence.  It is trying to
represent exactly the part of the coimage that decides the ground cell.

## Assumption challenge

The tempting vertex set is runners.  That is probably too concrete.  The better
vertices for this quotient are shell probes, lift denominators, endpoint-owner
obligations, relation-lattice vectors, CRT states, and proof obligations.  The
quotient preserves the `p_0=0` decision and the floor-witness/cap decision; it
destroys phase order and runner labels.

The challenged assumption is that cancellation should be explained by a low
moment, a density estimate, or a runner-level picture.  The alternative is that
`C=2n-1` probes are a conservative Yoneda basis for the floor coimage, and the
cancellation is exactly the resonance ledger seen through that basis.

## Link to HYP-2157

HYP-2157 says Helly rank is "how many overlap orders a proof must keep."  The
remote sibling reflection `coimage-yoneda-2n-minus-1-resonance-s605.md` names
the broad category/number-theory addendum.  This HYP-2161 reflection narrows it
to a theorem-shaped target: when no finite moment order suffices, prove a
conservative probe basis for the coimage.  For LRC at the floor, the candidate
basis is the `2n-1` shell/lift/CRT system.
