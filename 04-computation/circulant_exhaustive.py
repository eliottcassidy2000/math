"""
Exhaustive search over ALL circulant tournaments on Z_n to find max H.
Validates that QR_p is optimal for prime p ≡ 3 mod 4,
and finds the best circulant for non-Paley n.

For n vertices: 2^{(n-1)/2} connection sets to try.
Feasible for n ≤ 17 (2^8 = 256 options, fast DP).

Session: opus-2026-05-27-S7
"""

import subprocess
import sys
from itertools import combinations


def count_hp_bitmask(n, adj):
    """Bitmask DP for HP count. adj[v] = bitmask of v's out-neighbors."""
    full = (1 << n) - 1
    dp = [[0] * n for _ in range(full + 1)]
    for v in range(n):
        dp[1 << v][v] = 1
    total = 0
    for mask in range(1, full + 1):
        for v in range(n):
            if not (mask >> v & 1): continue
            d = dp[mask][v]
            if not d: continue
            if mask == full:
                total += d
                continue
            out = adj[v] & ~mask
            while out:
                lsb = out & -out
                w = lsb.bit_length() - 1
                dp[mask | lsb][w] += d
                out ^= lsb
    return total


def count_hp_circulant(n, S):
    """Circulant-reduced DP: O(2^{n-1}) memory, p passes."""
    sz = 1 << (n - 1)
    dp = [0] * sz
    dp[0] = n  # weight = n (all starting vertices equivalent)

    total = 0
    for length in range(1, n + 1):
        target_pc = length - 1
        for mask_hi in range(sz):
            if bin(mask_hi).count('1') != target_pc:
                continue
            cnt = dp[mask_hi]
            if not cnt:
                continue
            if length == n:
                total += cnt
                continue
            mask = (mask_hi << 1) | 1  # restore bit 0
            for d in S:
                if (mask >> d) & 1:
                    continue  # already visited
                # Shift by -d mod n
                new_mask = 1
                m = mask
                while m:
                    lsb = m & -m
                    v = lsb.bit_length() - 1
                    nv = (v - d + n) % n
                    if nv:
                        new_mask |= (1 << nv)
                    m ^= lsb
                dp[new_mask >> 1] += cnt
    return total


def is_prime(n):
    if n < 2: return False
    for i in range(2, int(n**0.5) + 1):
        if n % i == 0: return False
    return True


def paley_S(p):
    qr = set()
    for i in range(1, p):
        qr.add((i * i) % p)
    return sorted(i for i in range(1, p) if i in qr)


def exhaustive_circulant(n):
    """Try all 2^{(n-1)/2} circulant connection sets on Z_n."""
    half = (n - 1) // 2
    best_h = 0
    best_S = None

    # Generate all valid connection sets:
    # For each pair (i, n-i) for i=1..half, choose one of the two
    pairs = [(i, n - i) for i in range(1, half + 1)]
    count = 0
    for bits in range(1 << half):
        S = []
        for j, (a, b) in enumerate(pairs):
            if (bits >> j) & 1:
                S.append(a)  # choose lower element
            else:
                S.append(b)  # choose upper element (= n-i)

        h = count_hp_circulant(n, S)
        count += 1
        if h > best_h:
            best_h = h
            best_S = sorted(S)
            print(f"  New best: H={h}  S={S}", flush=True)

    print(f"  Tried {count} connection sets, best H={best_h}, S={best_S}")
    return best_h, best_S


def main():
    print("=== Exhaustive Circulant Tournament Search ===\n")

    # Known values from A038375
    known_a = {
        1: 1, 2: 1, 3: 3, 4: 5, 5: 15, 6: 45, 7: 189,
        8: 661, 9: 3357, 10: 15745, 11: 95095, 12: 531205, 13: 3719831
    }

    results = {}

    for n in range(3, 18, 2):  # odd n from 3 to 17
        print(f"--- n={n} ({2**((n-1)//2)} connection sets) ---")

        best_h, best_S = exhaustive_circulant(n)
        results[n] = (best_h, best_S)

        # Compare to Paley if applicable
        if is_prime(n) and n % 4 == 3:
            pS = paley_S(n)
            ph = count_hp_circulant(n, pS)
            status = "= Paley" if ph == best_h else f"< Paley ({ph})"
            print(f"  Paley QR_{n}: H={ph}  {status}")

        # Compare to known a(n)
        if n in known_a:
            a = known_a[n]
            print(f"  Known a({n}) = {a}", end="")
            if best_h == a:
                print("  ← BEST CIRCULANT = a(n)!")
            elif best_h < a:
                print(f"  ← circulant is SUBOPTIMAL by {a - best_h}")
            else:
                print(f"  ← NEW LOWER BOUND a({n}) ≥ {best_h}")

        print()

    print("\n=== Summary ===")
    print(f"{'n':4} {'Best circulant H':25} {'Best S':40} {'Notes'}")
    for n in sorted(results):
        h, S = results[n]
        notes = ""
        if n in known_a:
            if h == known_a[n]: notes = f"= a({n}) ✓"
            elif h < known_a[n]: notes = f"suboptimal (a={known_a[n]})"
        if is_prime(n) and n % 4 == 3:
            notes += " [Paley prime]"
        print(f"{n:4} {h:25} {str(S):40} {notes}")


if __name__ == "__main__":
    main()
