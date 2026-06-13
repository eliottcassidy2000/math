#!/usr/bin/env python3
r"""
lrc_tight_config_regular_orbit_chi_s581b.py    oracle-2026-06-03-S581o

Follow-up to s581. Core finding there: the rotational R_m (=AP) and Paley_m are DISTINCT
vertex-transitive regular tournaments with the SAME #3-cycles, separated only by the
DICHROMATIC number chi -- R_m has chi=2, Paley chi=3. So chi adds beyond vertex-transitivity.

Here: (A) verify chi(R_m)=2 (AP-rotational is the unique 'barely cyclic' regular orbit) for
m=5..13; (B) determine which regular orbit the LRC TIGHT configs realize at their tight time
-- does any tight config land on the chi=3 (Paley) orbit rather than the chi=2 (AP) one?
Tight family (S576o): n=6 -> {AP}; n=8 -> {AP, (1,2,3,4,5,7,12), (1,4,5,6,7,11,13)}.
"""
from itertools import permutations, combinations
from fractions import Fraction as Fr
from math import gcd

def circulant(m, conn):
    cs = set(c % m for c in conn)
    return [[1 if i != j and (j - i) % m in cs else 0 for j in range(m)] for i in range(m)]

def canon(adj, m):
    best = None
    for p in permutations(range(m)):
        key = tuple(adj[p[i]][p[j]] for i in range(m) for j in range(m) if i != j)
        if best is None or key < best: best = key
    return best

def is_acyclic(adj, S):
    for i, j, k in combinations(S, 3):
        f = adj[i][j] + adj[j][k] + adj[k][i]
        if f == 3 or f == 0: return False
    return True

def dichromatic(adj, m):
    for k in range(1, m + 1):
        color = [-1] * m
        def bt(v):
            if v == m: return True
            for c in range(k):
                color[v] = c
                if is_acyclic(adj, [u for u in range(v + 1) if color[u] == c]) and bt(v + 1):
                    return True
            color[v] = -1; return False
        if bt(0): return k
    return m

def regular(adj, m):
    d = (m - 1) // 2
    return all(sum(adj[i]) == d for i in range(m))

def num3(adj, m):
    return sum(1 for i, j, k in combinations(range(m), 3)
               if adj[i][j] + adj[j][k] + adj[k][i] in (0, 3))

# ---- LRC exact pinch + tight-time half-turn tournament ----
def circ(r, C): r %= C; return min(r, C - r)
def pinch_M(S):
    best = Fr(0); arg = None
    for C in set(a + b for i, a in enumerate(S) for b in S[i+1:]):
        for mm in range(1, C):
            v = Fr(min(circ(x * mm, C) for x in S), C)
            if v > best: best, arg = v, (mm, C)
    return best, arg

def runner_tournament(S, mm, C, sign):
    """half-turn tournament among runners at t=mm/C + sign*eps; ties (antipodal) broken by eps."""
    k = len(S)
    pos = [(x * mm) % C for x in S]
    adj = [[0] * k for _ in range(k)]
    half = Fr(C, 2)
    for a in range(k):
        for b in range(a + 1, k):
            diff = (pos[b] - pos[a]) % C
            if diff < half: ab = True
            elif diff > half: ab = False
            else: ab = (sign * S[a] < sign * S[b])   # antipodal tie -> eps direction
            if ab: adj[a][b] = 1
            else: adj[b][a] = 1
    return adj

def main():
    print("=" * 78)
    print("(A) chi(R_m)=2 ?  Is the AP-rotational the unique 'barely cyclic' regular orbit?")
    print("=" * 78)
    for m in (5, 7, 9, 11, 13):
        R = circulant(m, range(1, (m - 1) // 2 + 1))
        qr = sorted({(x * x) % m for x in range(1, m)})
        P = circulant(m, qr)  # Paley if m prime=3 mod4, else just the QR circulant
        iso = (canon(R, m) == canon(P, m)) if m <= 9 else "n/a(m>9)"
        print(f"  m={m:2d}: chi(R_m rotational/AP)={dichromatic(R, m)}  #3cyc={num3(R,m)};  "
              f"chi(QR/Paley circulant)={dichromatic(P, m)}  #3cyc={num3(P,m)};  iso? {iso}")

    print("\n" + "=" * 78)
    print("(B) which regular orbit do LRC TIGHT configs realize? (chi=2 AP vs chi=3 Paley)")
    print("=" * 78)
    tights = {
        "n=6 AP {1..5}": (6, (1, 2, 3, 4, 5)),
        "n=8 AP {1..7}": (8, (1, 2, 3, 4, 5, 6, 7)),
        "n=8 sporadic1 {1,2,3,4,5,7,12}": (8, (1, 2, 3, 4, 5, 7, 12)),
        "n=8 sporadic2 {1,4,5,6,7,11,13}": (8, (1, 4, 5, 6, 7, 11, 13)),
    }
    for name, (n, S) in tights.items():
        M, (mm, C) = pinch_M(S)
        k = len(S)
        info = []
        for sign in (+1, -1):
            adj = runner_tournament(S, mm, C, sign)
            reg = regular(adj, k)
            chi = dichromatic(adj, k)
            t3 = num3(adj, k)
            Rm = circulant(k, range(1, (k - 1) // 2 + 1))
            isR = (canon(adj, k) == canon(Rm, k))
            info.append((reg, chi, t3, isR))
        print(f"  {name}:  M={M}(={float(M):.4f}), tight time {mm}/{C}")
        for sign, (reg, chi, t3, isR) in zip(('+e', '-e'), info):
            print(f"      [{sign}] regular={reg}  chi={chi}  #3cyc={t3}  ==R_m(AP,chi2)? {isR}")

    print("\n" + "=" * 78)
    print("READING")
    print("=" * 78)
    print("""  If chi(R_m)=2 for all m while Paley/QR has chi=3, the AP-rotational is the UNIQUE
  minimally-cyclic (chi=2) regular orbit, and Paley is a genuinely-more-cyclic (chi=3)
  one that VT and 3-cycle-count cannot distinguish but chi can. If every LRC tight config
  realizes the chi=2 (R_m) orbit (not chi=3 Paley), then the LRC extremal is precisely the
  MINIMALLY-cyclic regular tournament -- not the maximally-symmetric Paley -- and chi=2 is
  a clean characterization of the tight orbit. Any tight config landing on chi=3 would be a
  regular-but-not-AP tight config with a DIFFERENT chi (the user's target).""")

if __name__ == "__main__":
    main()
