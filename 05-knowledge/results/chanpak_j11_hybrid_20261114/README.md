# Combining the exact eleven-block orbit with absolute source repairs

**Universal placement, scalar orbit cofiniteness, and Chan–Pak remain open.**
This packet proves that several constructions can repair one another,
and separates those repairs from local congruence saturation.

The eleven are the strict Wave-20 O/Z/J alphabet (not the distinct
13-port extraction or the eleven-offset bank). Its starting run is
capped at four, and every subsequent complete block is J followed by
one of

    empty, O, OO, OZ, OOO, OOZ, OZO, OZZ, OOOZ, OZOZ, OZZZ.

The inherited exhaustive range is 2^26 with 408 gaps, largest 418429.
[Wave-20 source](../chanpak_waves_20261113/README.md).

## Where the strategies actually join

- **[Exact root interface](ROOTED_BRIDGE.md).** A bilateral window
  can feed the eleven readers. Two different reader words have rows
  (9,9,−3) and (9,9,−1), giving a unit odd step on every literal source
  with d=1. Changing both the source and the reader gives the three-odd
  window 3456uv+459u+2430v+{319,321,323} for independent dyadic u,v>=2.
  Three of the four earlier free-height offsets survive a J reader
  exactly; the fourth misses by a computable change of rooted W.
  [Literal replay](rooted_window.py).
- **[Absolute module repair](TWO_VISIT_WINDOWS.md).** Transitive blow-ups
  produce the five consecutive odds 2^k+{1,3,5,7,9} at every k>=6.
  This includes actual eleven-block misses, in particular 262147; a second
  ray 76·2^n+21 includes the miss 38933, and 2^k+21 includes 65557.
  The repeated-module terminal has at most two visits, and the path
  cover count is exact rather than a residue bound.
  [Literal replay and target decoder](replay_power_three.py).
- **[Finite head union](HYBRID_HEAD.md).** Uncapping only the source and
  using at most two of the same readers covers 363 of the 408 gaps;
  the other port-chain alphabet and three literal punctures finish
  the admissible head. This verifies real complementarity, not an
  extension to a larger numerical range. [Replay](head_hybrid.py).
- **[Changed readers at independent heights](CHANGING_READERS.md).**
  A fixed source phase repair can be absorbed into an enlarged reader,
  or transported into the two-port reader modes and through the
  shared J move. The two-port alphabet is marked separately.
  [Replay](replay_readers.py).

## Why the tail is not a congruence problem

- **[Joint local saturation inside the actual eleven](SIMULTANEOUS_LOCAL.md).**
  Seven short seeds followed just by J supply every odd target prime
  to 3 modulo arbitrary simultaneous 2^a3^b5^c7^d. This covers the
  congruence conditions on the critical prime and 7p lanes, regardless
  of precision. The digit decoder finds congruent representatives
  even for the actual orbit miss 67, but the literal height can be
  enormous. [Replay](local_11.py).
- **[Rank separators and ordinal products](RANK_SEPARATORS.md).**
  Four explicit nonlocal rank-one patterns break a candidate count
  into two factors greater than 7, ruling them out on both critical
  lanes. The pure reset on entering d=0 can be removed without
  changing the multiplicative closure. It cannot be thrown away
  merely because it is singular; the actual split and normalized
  source are the reason. [Replay](replay_rank.py).
- **[What the orbit cannot supply](ELEVEN_ORBIT_LIMITS.md).** The
  full-state count through X is O(X^1.85), even after freeing the
  one-apex starting word. This rules out using the orbit as a full
  two-dimensional interface reservoir, not scalar cofiniteness.
  Full-rank products cannot provide an exact positive integer-radix
  self-loop, and a single transitive seam cannot add exactly two to
  an existing large count. [Exact controls](eleven_cone_controls.py).

The next proof target is an **absolute** decreasing selector for primes
and 7p in the reduced orbit, with alternative source or module charts
on its true failures. The local result means any such selector has to
remember height/magnitude or richer arithmetic, not merely more digits
at 2,3,5,7. The sparse absolute windows isolate several real repairs,
but do not cover the unspecified infinite residual.
