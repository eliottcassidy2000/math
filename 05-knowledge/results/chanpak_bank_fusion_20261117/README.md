# Fusion of the growing Chan–Pak bank with rooted and sevenfold constructions

The previous variable-height three-contact bank admits all sorts of fine
shifts but does not jointly choose their absolute source. Here two physical
combinations work exactly, and two weighted bounds explain why natural
recursive enlargements still do not prove coverage.

## Positive results

- [Root migration and a neutral pair](migration/PROOF.md): changing the
  distinguished root preserves its unique ending path, so arbitrarily long
  rooted constructions become legitimate inputs to the reflected +2 bank.
  One 5-vertex tournament with two port choices has the same full unmarked
  slot polynomial for both, then produces 15 and 17 under the same kind of
  rooted migration. After optional transitive padding the whole response
  difference is q²[(q+1)^(t+1)−(q−1)^(t+1)]. In the actual reflected bank,
  changing this one rooted port moves both tracks by 294−10u+512w, so every
  prescribed even shift has infinitely many physical witnesses with +2
  intact. The choices still lie on a constrained skew line.
- [Sevenfold fusion](readers/PROOF.md): an equal-order +2 pair cannot have
  equal P₂, by reduction modulo 3. Use the smallest compatible gap, +6.
  A rational-base carry loop pumps a strong, equal-P₂ +6 pair to arbitrary
  order; the prior sevenfold exterior gives the exact four-corner family
  N,N+6,7N,7N+42. There is a certified prime example
  N=2466476862736151977 with its physical 7N. This is not a prime infinitude
  or all-prime theorem.
- [Prime-valued collision and exact subsequence](arithmetic/PROOF.md):
  the earlier joint-placement collision line actually attains the prime
  42284279505304223759 in two continuation-distinct ways. A quadratic
  dyadic subfamily has an exact target inverse, no fixed prime divisor,
  and infinitely many distinct prime factors collectively. It is not
  asserted to take infinitely many prime values.

## Proven limits

- [Recursive rooted enlargement](coupling/PROOF.md): unrestricted depth
  and a rigorously bounded frequency of weak stages preserve the +2, but
  the repair endpoint library is still O(X^(.99) log X). This deliberately
  leaves unrestricted weak stages outside the proof. The
  [weak-tail normal form](migration/WEAK.md) supplies an exact
  divisibility filter and lets a target choose a dyadic tail when
  the preceding root state is known.
- [Positive recursive fusion](closure/PROOF.md): even the original growing
  bank plus a high-count rooted tail, mixed at arbitrary depth by ordinal
  joins, cyclic substitutions and unbounded prescribed transitive clock
  blocks, realizes at most O(X^(.995)) different counts in the stated
  grammar. No such estimate is claimed for arbitrary tagged or signed glues.

These splits suggest a concrete next objective: use the rooted
difference (unseen by all the pure slot measurements) to select the
low error inside a target-selected high count, rather than further
increase the number of positive bank compositions.
**Chan–Pak remains open.**

## Replays

    python3 migration/verify.py
    python3 readers/replay.py
    python3 arithmetic/verify.py
    python3 coupling/replay.py

The arithmetic replay imports the earlier published joint-placement
certificate at its relative repository path. All other replays in
this list use only the standard library and directly encode their
small literal graph checks.
