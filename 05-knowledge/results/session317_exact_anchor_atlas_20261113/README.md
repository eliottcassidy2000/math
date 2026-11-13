# Session 317: Exact anchors, changing readouts, and the limits of static digit banks

The research question was not whether larger searches keep finding
witnesses. It was which part of a constructive proof is genuinely
missing after the previous sessions' local count switches. The answer
here is more specific: **choose an exact anchor, and control a menu
of legal readouts there, rather than asking one endpoint state to
support every operation.**

## Constructive results

* [The two-apex terminal atlas](TERMINAL_ATLAS.md) reduces **any
  ending after a binary prefix** to four physical path-count quantities.
  Three literal endings give strong tournaments with consecutive odd
  counts $254B-9K(B)+\{91,93,95\}$ at the same prefix $B$, at unbounded
  scale. The coefficient of the nonlinear binary defect is always
  negative and cannot be removed by choosing a longer ending.
* [Exact target inversion](EXACT_INVERSE.md) uses the bit-reversal seam
  to find all prefixes for a fixed terminal in square-root rather than
  linear-in-target search, up to logarithms. In the other direction,
  a prefix near the square root of the target limits the end response
  to a bounded interval. A fixed reservoir count supports only
  polylogarithmically many terminal offsets at fixed response.
* [The one-return-star factorization](ONE_RETURN_STAR.md) separates
  the two path slots and the signed exposed boundary. Cancellation
  gives the strong count family
  $(2^{k+1}+1)(2^k+1)2^m+1$ at fixed intercept 1, escaping a
  one-sided dyadic slope obstruction.
* [Coalescent menus](COALESCENT_MENUS.md) give one exact strong count
  with $r+1$ distinct next $Q(2)$ values for every $r$, by
  preserving $Q(2)/H$ while changing $Q(3)/H$ in complementary
  source factors. These menus are not at arbitrary anchors.

## Design constraints, not just failed examples

* [Rational radix and capacity](RADIX_AND_CAPACITY.md): mixed powers
  of two coprime bases really can represent **every integer** with a
  bounded digit alphabet and overlapping windows. But **every**
  fixed one-sided additive tournament type has sparse absolute
  support, no matter how many endpoint-matched cuts it has. A
  weighted capacity inequality extends to finite changing-type
  banks. This identifies exactly why the old modular selectors
  alone cannot close Chan–Pak.
* [Slot invariants](SLOT_INVARIANTS.md): fixed $Q(2)$ rules out two
  successive $+2$ steps in any grammar, even when order parity
  can change. Multiple cyclic readiness returns force much stronger
  divisibilities. Literal changing-interface arrivals at 1431 show
  how to preserve the required combination instead of each input.
* [Contact capacity](CONTACT_CAPACITY.md) proves marked-stabilizer
  and arithmetic slot-divisor bounds, with sharp strong constructions
  even for rigid blocks. Homogeneous strong multiplication must
  overshoot every fixed factor once the source is large.
  [A second exact two-ended aperture](TWO_END_APERTURE.md) gives
  strong $20\cdot2^L+17,20\cdot2^L+19$ pairs, but changes $Q(2)$.

The older Gold prefix guards, lonely-runner rational-base sheets,
Collatz bit reversal, Pell arithmetic, and the square-sum hinge were
used as ways to choose the retained state. None was treated as an
implication between conjectures. The most concrete next step is an
exact terminal compiler that **changes the reservoir count** while
filling the narrow end-response band, or an interacting construction
that beats the additive capacity bound. The proofs above rule out
trying to do this by more padding of one fixed state.

## Reproduce

Run the companion `replay_*.py` files individually from this folder.
The independent optional `replay_contact_literal.cpp` checks the
19-vertex contact example. The replays test literal tournaments or
integer decoders; the all-order statements use the proofs.

**Chan–Pak remains open.** Consecutive triples and arbitrarily large
same-count menus are exact, but neither fills the gaps between
absolute anchors.
