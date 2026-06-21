"""
lrc14_paley_p19_actual_max_angle6.py

SHOCK from the sample test: at p=19 a SAMPLED circulant beat Paley's H, and H was NOT
Schur-concave (121 violations). This script does the HONEST exhaustive check:

  - Compute H for ALL 2^9 = 512 circulant tournaments on Z_19 (2^19 DP each; ~512*19*2^19 ops).
  - Is Paley (QR_19) the global circulant H-max at p=19, or does it LOSE?
  - If Paley loses, the whole "Paley H-max" premise breaks at p=19 (independent of my mechanism).
  - Report the true argmax, its connection set, its spectrum spread, and whether the winner is a
    difference set / flat spectrum.
  - Recompute Schur-concavity exhaustively at p=19.

This is the falsification gate: if Paley is NOT even the H-max at p=19, then the conjecture
"Paley = circulant H-max for all p=3 mod 4" is FALSE and my majorization mechanism only governs
p where Paley genuinely wins (7,11). Honest.

Run: python3 04-computation/lrc14_paley_p19_actual_max_angle6.py
"""
import sys, cmath, math

sys.stdout.reconfigure(line_buffering=True)


def qr_set(p):
    return frozenset((a * a) % p for a in range(1, p))


def tournament_connection_sets(p):
    m = (p - 1) // 2
    pairs = [(t, (p - t) % p) for t in range(1, m + 1)]
    sets = []
    for bits in range(1 << m):
        S = set()
        for i, (a, b) in enumerate(pairs):
            S.add(a if (bits >> i & 1) else b)
        sets.append(frozenset(S))
    return sets


def ham_paths_circulant(S, p):
    n = p
    adj = [0] * n
    for i in range(n):
        mask = 0
        for j in range(n):
            if i != j and (j - i) % n in S:
                mask |= (1 << j)
        adj[i] = mask
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1, 1 << n):
        row = dp[mask]
        for v in range(n):
            if not (mask >> v & 1):
                continue
            cnt = row[v]
            if cnt == 0:
                continue
            avail = adj[v] & ~mask
            mm = avail
            while mm:
                u = (mm & -mm).bit_length() - 1
                dp[mask | (1 << u)][u] += cnt
                mm &= mm - 1
    return sum(dp[(1 << n) - 1])


def spectrum_spread(S, p):
    omega = cmath.exp(2j * math.pi / p)
    vec = [abs(sum(omega ** (k * s) for s in S)) ** 2 for k in range(1, p)]
    return max(vec) - min(vec)


def autocorr_flat(S, p):
    c = [0] * p
    for a in S:
        for b in S:
            c[(a - b) % p] += 1
    return len(set(c[t] for t in range(1, p))) == 1


def main():
    p = 19
    QR = qr_set(p)
    sets = tournament_connection_sets(p)
    print(f"p={p}: exhaustive H over {len(sets)} circulant tournaments ...")
    data = []
    for idx, S in enumerate(sets):
        H = ham_paths_circulant(S, p)
        data.append((H, S))
    data.sort(key=lambda x: -x[0])
    maxH = data[0][0]
    pal_H = ham_paths_circulant(QR, p)
    print(f"  Paley QR_19 H = {pal_H}")
    print(f"  global circulant max H = {maxH}")
    print(f"  Paley IS the circulant H-max at p=19: {pal_H == maxH}")
    print(f"  top 5 circulants:")
    for H, S in data[:5]:
        print(f"    H={H}  S={sorted(S)}  flat_spectrum={autocorr_flat(S,p)}  "
              f"spread={spectrum_spread(S,p):.3f}")
    # Paley rank
    rank = next(i for i, (H, S) in enumerate(data) if S == QR) + 1
    print(f"  Paley rank among 512 circulants: {rank}")

    # how far is the winner from Paley structurally?
    winner = data[0][1]
    print(f"  winner flat? {autocorr_flat(winner, p)}, "
          f"winner spread {spectrum_spread(winner,p):.4f} vs Paley spread {spectrum_spread(QR,p):.2e}")


if __name__ == "__main__":
    main()
