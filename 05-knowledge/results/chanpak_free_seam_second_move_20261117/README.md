# Second physical moves for the free +2 Chan–Pak seam

The [original seam](../chanpak_wheel_square_selector_20261117/FREE_SEAM.md)
has an exact +2 reflection on two freely chosen binary contacts, and its
initial absolute placement error was O(N^(2/5)). This packet separates
three issues: whether a second actual move exists, how close the first
placement can get, and whether the remaining error always falls within
the realizable second moves.

## What is now proved

1. **One arc can give an independent second direction.**
   The [one-arc switch](one_arc/PROOF.md) gives an all-height physical
   commuting square with vertical +2 and any prescribed even horizontal
   displacement. For displacement 4 the literal 19-vertex fixture gives
   261163, 261165, 261167 and 261169; each horizontal edge is one apex
   arc reversal. Small displacements require tuned, skew address splits.
   All one-flank, leading-bit- and parity-preserving single-arc repairs
   of size up to E have only O((E+1)sqrt(X)log X) possible endpoints
   below X for E<=sqrt(X). In particular that mechanism alone cannot
   universally correct an error of order N^(1/3).
2. **The baseline can be chosen closer without spending the +2.**
   The [complementary-carry selector](arithmetic/SELECTOR.md) improves
   the uniform exponent to 1/3, and gives 1/4 on explicit infinitely
   many target bands. It uses sign bracketing, not monotonicity.
   The [combined phase selector](phase/PROOF.md) proves that the
   *same* approximate witness can also match any fixed prescribed
   modulus M with error O_M(N^(1/3)); the constants are not uniform in M.
3. **Growing physical banks, and the obstruction to treating one as a
   generic fine dial.** [Paired-digit moves](paired/PROOF.md) shift both
   counts by +512 or +522 at all common prefix heights. The entire
   same-width common-prefix mechanism has sharp minimum nonzero shift
   512, so it cannot simply be iterated as an unrestricted Bezout dial.
   The [middle-collar grammar](collars/PROOF.md) grows a different
   center without losing the signed switch. Its four-state radix-8
   recursion has an exact count of legal words, but a growing spacing
   in the possible translations at fixed height.
4. **Two ways out of that restriction.** The [Farey flank
   construction](extra/PROOF.md) makes the gap a determinant of
   marked endpoint counts; strong interiors of arbitrarily large order
   and reflected guards can preserve it. The [full-profile rotor](bank/PROOF.md)
   instead retags one actual center, giving a three-term progression
   with the *entire* signed flank interface fixed, not just a balanced
   unsigned read. Neither construction yet programs the required
   absolute residual for every target.

The next useful question is not whether any two moves commute: they
do. It is whether a coordinated **multi-arc, or changed-module**
construction has small signed moves at near-balanced address scales,
where the improved baseline selector operates. The one-arc counting
bound and exact common-prefix barrier show why that distinction is
necessary. No universal equality, density-one claim, or proof of
Chan–Pak is made.

## Reproduction

Python replays require only the standard library and the previous
free-seam packet at its checked-in relative path:

    python3 arithmetic/check.py
    python3 phase/verify_phase.py
    python3 paired/verify.py
    python3 collars/collar.py
    python3 bank/replay.py
    python3 bank/steer.py
    python3 extra/replay.py
    g++ -O3 -std=c++20 one_arc/replay.cpp -o /tmp/one_arc && /tmp/one_arc
    g++ -O3 -std=c++20 extra/apex_no_go.cpp -o /tmp/apex_no_go && /tmp/apex_no_go

The arithmetic and phase controls are exact integer evaluations.
The four-corner switch and small witnesses use the actual tournament
arcs and independent subset path DP. Separate finite tests are
controls for the proofs, not a substitute for all-order reasoning.
