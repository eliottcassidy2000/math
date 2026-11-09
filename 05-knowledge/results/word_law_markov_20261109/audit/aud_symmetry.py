"""THM-5020 Remarks / README section 3: 'the construction commutes with sigma and with complementation'.
Checks, for every word of length L = n - 2 (n <= 16):
  sigma(R(w)) = R(reverse w)                       (commutes with sigma)
  R(complement w) = R(w) + 1 ?                     (would be 'commutes with complementation')
and lists the words where the second identity holds.
usage: python aud_symmetry.py
"""
import itertools


def runwindow(w, n):
    rev = set()
    a, L = 0, len(w)
    while a < L:
        if w[a] != '1':
            a += 1
            continue
        p = a
        while a < L and w[a] == '1':
            a += 1
        q = a - 1
        for i in range(p, q + 3):
            for j in range(i + 2, q + 3):
                rev.add((i, j))
    return rev


def main():
    for n in range(3, 17):
        L = n - 2
        tiles = {(i, j) for i in range(n) for j in range(i + 2, n)}
        sig_bad = 0
        comp_ok = []
        for bits in itertools.product('01', repeat=L):
            w = ''.join(bits)
            R = runwindow(w, n)
            sigR = {(n - 1 - j, n - 1 - i) for (i, j) in R}
            if sigR != runwindow(w[::-1], n):
                sig_bad += 1
            wb = ''.join('1' if c == '0' else '0' for c in w)
            if runwindow(wb, n) == tiles - R:
                comp_ok.append(w)
        print(f"n = {n}: sigma(R(w)) != R(reverse w) for {sig_bad} words; R(complement w) = R(w) + 1 for "
              f"{len(comp_ok)} of {2 ** L} words: {comp_ok[:6]}")


if __name__ == "__main__":
    main()
