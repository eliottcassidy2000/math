"""THM-5020 (5), central-pair lemma: if w is a palindrome of even length L and P(X) is odd (X = A or B), then the two
central positions L/2 - 1 and L/2 lie in X. So P(A) and P(B) are not both odd. Checked for every palindrome of even
length L <= 24, together with the central-position statement.
usage: python palindrome_parity_check.py LMAX
"""
import sys


def P(xs):
    c, last = 0, -10
    for x in xs:
        if x - last >= 3:
            c += 1
            last = x
    return c


def main():
    LMAX = int(sys.argv[1])
    for L in range(2, LMAX + 1, 2):
        h = L // 2
        bad = both = 0
        for code in range(1 << h):
            half = format(code, f'0{h}b')
            w = half + half[::-1]
            for col in '01':
                X = [i for i in range(L) if w[i] == col]
                if P(X) % 2 == 1 and not (w[h - 1] == col and w[h] == col):
                    bad += 1
            A = [i for i in range(L) if w[i] == '1']
            B = [i for i in range(L) if w[i] == '0']
            both += (P(A) % 2 == 1 and P(B) % 2 == 1)
        print(f"L = {L}: palindromes {1 << h}; odd P(X) without both central positions in X: {bad}; "
              f"P(A), P(B) both odd: {both}")


if __name__ == "__main__":
    main()
