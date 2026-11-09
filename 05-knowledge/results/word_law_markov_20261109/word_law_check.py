"""Brute-force checks of THM-5020 (the word law for tight lines at every order), opus-2026-11-09-S88.

Convention (THM-4951): path 0 -> 1 -> ... -> n-1; tile (i, j), j >= i + 2, reversed in t iff j -> i; u = t + 1.
A word w of length L = n - 2 colours the positions a = 0..L-1; A = {a : w_a = 1}, B = {a : w_a = 0}.
P(X) = largest subset of X with pairwise differences >= 3 (greedy).
R(w) = run-window tiling: the union of the staircases on the windows [p, q + 2] of the 1-runs [p, q] of w.
A junction is an isolated letter whose two neighbouring runs (of the other letter) both have length >= 2;
J0 = number of isolated-0 junctions, J1 = number of isolated-1 junctions, J = J0 + J1.

Checked for every word of length L = n - 2, n = 4..N_ALL:
  (L2) nu(R(w) + 1) = P(B)                                  [every word]
  (L1) nu(R(w)) = P(A)                                       [words with J0 = 0]
  (J)  P(A) + P(B) >= ceil((L + 2J)/3)                       [every word]
  (W)  for every tight colouring (P(A) + P(B) = floor(n/3)): J <= 1, and the line of R(w) (if J0 = 0) or of
       R(complement w) (if J0 = 1) is tight and has an endpoint with triple word w; it is blue iff w is a palindrome.
Then (W) alone for every tight colouring up to n = N_TIGHT (tight colourings generated directly).
usage: python word_law_check.py N_ALL N_TIGHT
"""
import sys
import itertools


def tournament(n, rev):
    out = [0] * n
    for i in range(n - 1):
        out[i] |= 1 << (i + 1)
    for i in range(n):
        for j in range(i + 2, n):
            if (i, j) in rev:
                out[j] |= 1 << i
            else:
                out[i] |= 1 << j
    return out


def nu(out, n):
    """maximum number of vertex-disjoint 3-cycles (= of disjoint cycles)"""
    tri = [[] for _ in range(n)]
    for a in range(n):
        for b in range(a + 1, n):
            for c in range(a + 1, n):
                if b != c and (out[a] >> b & 1) and (out[b] >> c & 1) and (out[c] >> a & 1):
                    tri[a].append((1 << a) | (1 << b) | (1 << c))
    best = [0]

    def rec(av, cur):
        if cur > best[0]:
            best[0] = cur
        if not av or cur + bin(av).count("1") // 3 <= best[0]:
            return
        lo = (av & -av)
        a = lo.bit_length() - 1
        for T in tri[a]:
            if T & av == T:
                rec(av & ~T, cur + 1)
        rec(av & ~lo, cur)
    rec((1 << n) - 1, 0)
    return best[0]


def P(xs):
    c, last = 0, -10
    for x in xs:
        if x - last >= 3:
            c += 1
            last = x
    return c


def runs_of(w):
    out, i, L = [], 0, len(w)
    while i < L:
        s = i
        while i < L and w[i] == w[s]:
            i += 1
        out.append((w[s], s, i - 1))
    return out


def junctions(w, iso):
    R = runs_of(w)
    return sum(1 for j in range(1, len(R) - 1)
               if R[j][0] == iso and R[j][1] == R[j][2]
               and R[j - 1][2] - R[j - 1][1] >= 1 and R[j + 1][2] - R[j + 1][1] >= 1)


def run_window(w, n):
    return {(i, j) for (ch, p, q) in runs_of(w) if ch == '1'
            for i in range(p, q + 3) for j in range(i + 2, q + 3)}


def complement_tiles(n, rev):
    return {(i, j) for i in range(n) for j in range(i + 2, n)} - rev


def word_of(n, rev):
    return ''.join('1' if (a, a + 2) in rev else '0' for a in range(n - 2))


def sigma_fixed(n, rev):
    return all(((n - 1 - j, n - 1 - i) in rev) == ((i, j) in rev) for i in range(n) for j in range(i + 2, n))


def check_tight(w, n):
    """returns None if the word-law construction works for the tight colouring w, else a message"""
    k = n // 3
    J0, J1 = junctions(w, '0'), junctions(w, '1')
    if J0 + J1 > 1:
        return f"tight colouring with J = {J0 + J1}"
    wb = ''.join('1' if c == '0' else '0' for c in w)
    base = w if J0 == 0 else wb
    rev = run_window(base, n)
    t, u = tournament(n, rev), tournament(n, complement_tiles(n, rev))
    a, b = nu(t, n), nu(u, n)
    if a + b != k:
        return f"not tight: nu = ({a}, {b}), k = {k}"
    ends = {word_of(n, rev), word_of(n, complement_tiles(n, rev))}
    if w not in ends:
        return "word not realised"
    if sigma_fixed(n, rev) != (w == w[::-1]):
        return "colour mismatch"
    return None


def main():
    n_all, n_tight = int(sys.argv[1]), int(sys.argv[2])
    bad = {'L2': 0, 'L1': 0, 'J': 0, 'W': 0}
    words = 0
    for n in range(4, n_all + 1):
        L, k = n - 2, n // 3
        for bits in itertools.product('01', repeat=L):
            w = ''.join(bits)
            words += 1
            A = [i for i in range(L) if w[i] == '1']
            B = [i for i in range(L) if w[i] == '0']
            pa, pb = P(A), P(B)
            J0, J1 = junctions(w, '0'), junctions(w, '1')
            if pa + pb < -(-(L + 2 * (J0 + J1)) // 3):
                bad['J'] += 1
            rev = run_window(w, n)
            t, u = tournament(n, rev), tournament(n, complement_tiles(n, rev))
            if nu(u, n) != pb:
                bad['L2'] += 1
            if J0 == 0 and nu(t, n) != pa:
                bad['L1'] += 1
            if pa + pb == k and check_tight(w, n):
                bad['W'] += 1
    print(f"all words, n = 4..{n_all}: {words} words; failures {bad}")
    # tight colourings only, generated by a pruned search over prefixes (P is monotone in the prefix)
    for n in range(n_all + 1, n_tight + 1):
        L, k = n - 2, n // 3
        cnt = fails = blue = 0
        stack = ['']
        while stack:
            pre = stack.pop()
            A = [i for i in range(len(pre)) if pre[i] == '1']
            B = [i for i in range(len(pre)) if pre[i] == '0']
            if P(A) + P(B) > k:
                continue
            if len(pre) == L:
                if P(A) + P(B) == k:
                    cnt += 1
                    blue += pre == pre[::-1]
                    msg = check_tight(pre, n)
                    if msg:
                        fails += 1
                        if fails <= 3:
                            print("  FAIL", n, pre, msg)
                continue
            stack.append(pre + '0')
            stack.append(pre + '1')
        print(f"n = {n}: tight colourings {cnt} (palindromic {blue}); word-law construction failures {fails}")


if __name__ == "__main__":
    main()
