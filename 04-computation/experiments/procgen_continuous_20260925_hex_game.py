#!/usr/bin/env python3
"""procgen_continuous_20260925_hex_game.py

Lane: natural boundaries, Mahler/Cobham and harmonic trees
(session collatz-procgen-20260922, 2026-09-25).

Part D: Hex/Brouwer (Gale 1979) versus Althofer's 3n+-1 game (Conway's Beans-Don't-Talk).

  HG1  No topological no-draw theorem is possible: from every odd n >= 3 there is an infinite
       play.  Of 3n+1 and 3n-1 exactly one is 2 mod 4, and that move leads to (3n+-1)/2 > n,
       an odd number; repeating it never reaches 1.  (PROVED; checked for odd n < 10^6, and a
       10^4-move ascending play is exhibited.)  In Hex every complete play has a winner; here
       'no draw' can only be a statement about OPTIMAL play.
  HG2  Strategy stealing has no purchase: the mover loses at P-positions, which exist (7 is the
       first).  Capped retrograde values for odd n < 200.

Every check raises on failure.  Runtime a few seconds.
"""
import sys
import time

T0 = time.time()


def check(cond, msg):
    if not cond:
        raise AssertionError("CHECK FAILED: " + msg)


def oddpart(m):
    while m % 2 == 0:
        m //= 2
    return m


def main():
    print("=" * 78)
    print("HG1  Every odd n >= 3 has an infinite (never-ending) play")
    print("=" * 78)
    for n in range(3, 10**6, 2):
        a, b = 3 * n + 1, 3 * n - 1
        check((a % 4 == 2) != (b % 4 == 2), "exactly one of 3n+-1 is 2 mod 4")
        up = a if a % 4 == 2 else b
        check((up // 2) % 2 == 1 and up // 2 > n, "the 2-mod-4 move gives an odd (3n+-1)/2 > n")
    x = 3
    for _ in range(10**4):
        a, b = 3 * x + 1, 3 * x - 1
        x = (a if a % 4 == 2 else b) // 2
        check(x != 1, "ascending play never reaches 1")
    print("  for all odd 3 <= n < 10^6: exactly one of 3n+1, 3n-1 is 2 mod 4 and it yields the odd")
    print("  successor (3n+-1)/2 > n.  A 10^4-move ascending play from 3 ends at a number with "
          f"{len(str(x))} digits.")
    print("  => infinite plays exist from every odd n >= 3; unlike Hex, a 'complete play' need not")
    print("     have a winner, so Gale's Hex <=> Brouwer mechanism (every filled board has exactly one")
    print("     winner) has no analogue.  The prize question is determinacy of optimal play.")

    print()
    print("=" * 78)
    print("HG2  P-positions exist (the mover can lose): capped retrograde for odd n < 200")
    print("=" * 78)
    C = 2**18
    val = {}
    children = {}
    for n in range(1, C + 1, 2):
        children[n] = (oddpart(3 * n + 1), oddpart(3 * n - 1)) if n > 1 else (1, 1)
    # least fixed point: N if some child is 1 or P; P if both children are N (children above C unknown)
    pending = list(range(1, C + 1, 2))
    rounds = 0
    while True:
        rounds += 1
        nxt = []
        for n in pending:
            c1, c2 = children[n]
            v1 = 'P' if c1 == 1 else val.get(c1)
            v2 = 'P' if c2 == 1 else val.get(c2)
            if v1 == 'P' or v2 == 'P':
                val[n] = 'N'
            elif v1 == 'N' and v2 == 'N':
                val[n] = 'P'
            else:
                nxt.append(n)
        if len(nxt) == len(pending):
            break
        pending = nxt
    small = [n for n in range(1, 200, 2)]
    Ps = [n for n in small if val.get(n) == 'P']
    Ns = [n for n in small if val.get(n) == 'N']
    U = [n for n in small if n not in val]
    print(f"  retrograde to cap {C} ({rounds} rounds): among odd n < 200: {len(Ns)} N, {len(Ps)} P, {len(U)} unresolved")
    print(f"  P-positions below 200: {Ps}")
    check(7 in Ps and 1 not in Ps, "7 is a P-position")
    check(val.get(7) == 'P' and val.get(5) == 'N' and val.get(11) == 'N', "7 -> {11, 5}, both N")
    print("  7 is P: its moves go to 11 and 5, and both are N (5 -> 1 wins; 11 -> 32 -> 1 wins).")
    print("  So the player to move loses at about half the positions (repo: P-density about 0.48 below")
    print("  2^32); no strategy-stealing argument gives a first-player win, and none addresses draws.")
    print()
    print(f"ALL CHECKS PASSED  ({time.time() - T0:.1f}s)")


if __name__ == "__main__":
    main()
