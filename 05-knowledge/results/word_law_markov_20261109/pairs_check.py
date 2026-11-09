"""PREREGISTRATION Q1-Q3: allowed nu-pairs from admissible run-window words, n = 13..40 (opus-2026-11-09-S88).

A word is admissible if it avoids the factor 11011 (no 0-junction; then R(w) has nu = (P(A), P(B)) by (L1), (L2)) or
avoids 00100 (no 1-junction; then the complement's run-window line does). Palindromic admissible words give blue
lines; non-palindromic ones give black lines.

All admissible words are counted by a DP whose states are:
  s_A, s_B   greedy distance states;
  cA, cB     greedy counts;
  last4      the last four letters;
  f0, f1     whether 11011 / 00100 has occurred.
Palindromic admissible words are enumerated through their first half.
For each n and each unordered pair {i, j} the script reports whether a non-palindromic admissible word (black) and a
palindromic one (blue) attain it. It then checks Q1 (black: 1 <= i <= j <= k, i + j >= k) and Q2 (blue: i + j > k).
usage: python pairs_check.py NMIN NMAX
"""
import sys
from collections import defaultdict


def greedy_step(sA, sB, cA, cB, col):
    if col == 1:
        if sA >= 3:
            cA, sA = cA + 1, 1
        else:
            sA = min(3, sA + 1)
        sB = min(3, sB + 1)
    else:
        if sB >= 3:
            cB, sB = cB + 1, 1
        else:
            sB = min(3, sB + 1)
        sA = min(3, sA + 1)
    return sA, sB, cA, cB


def all_admissible_pairs(L):
    """dict (cA, cB) -> number of admissible words of length L"""
    st = defaultdict(int)
    st[(3, 3, 0, 0, '', False, False)] = 1
    for _ in range(L):
        ns = defaultdict(int)
        for (sA, sB, cA, cB, last, f0, f1), c in st.items():
            for col in (0, 1):
                a, b, x, y = greedy_step(sA, sB, cA, cB, col)
                s = last + str(col)
                g0 = f0 or s[-5:] == '11011'
                g1 = f1 or s[-5:] == '00100'
                if g0 and g1:
                    continue
                ns[(a, b, x, y, s[-4:], g0, g1)] += c
        st = ns
    out = defaultdict(int)
    for (sA, sB, cA, cB, last, f0, f1), c in st.items():
        out[(cA, cB)] += c
    return out


def P(xs):
    c, last = 0, -10
    for x in xs:
        if x - last >= 3:
            c += 1
            last = x
    return c


def palindromic_pairs(L):
    out = defaultdict(int)
    h = L // 2
    mids = ('',) if L % 2 == 0 else ('0', '1')
    for code in range(1 << h):
        half = format(code, f'0{h}b') if h else ''
        for mid in mids:
            w = half + mid + half[::-1]
            if '11011' in w and '00100' in w:
                continue
            A = [i for i in range(L) if w[i] == '1']
            B = [i for i in range(L) if w[i] == '0']
            out[(P(A), P(B))] += 1
    return out


def main():
    nmin, nmax = int(sys.argv[1]), int(sys.argv[2])
    q1 = q2 = True
    for n in range(nmin, nmax + 1):
        L, k = n - 2, n // 3
        allp = all_admissible_pairs(L)
        palp = palindromic_pairs(L)
        tot = defaultdict(int)
        pal = defaultdict(int)
        for (a, b), c in allp.items():
            tot[(min(a, b), max(a, b))] += c
        for (a, b), c in palp.items():
            pal[(min(a, b), max(a, b))] += c
        miss_black = [(i, j) for i in range(1, k + 1) for j in range(i, k + 1)
                      if i + j >= k and tot[(i, j)] - pal[(i, j)] <= 0]
        miss_blue = [(i, j) for i in range(1, k + 1) for j in range(i, k + 1)
                     if i + j > k and pal[(i, j)] <= 0]
        if miss_black:
            q1 = False
        if miss_blue:
            q2 = False
        print(f"n = {n}, k = {k}: admissible words {sum(allp.values())}, palindromic {sum(palp.values())}; "
              f"black pairs missing {miss_black}; blue (non-tight) pairs missing {miss_blue}")
    print(f"Q1 (black, all allowed pairs realised): {'CONFIRMED' if q1 else 'REFUTED'}")
    print(f"Q2 (blue, all non-tight allowed pairs realised): {'CONFIRMED' if q2 else 'REFUTED'}")


if __name__ == "__main__":
    main()
