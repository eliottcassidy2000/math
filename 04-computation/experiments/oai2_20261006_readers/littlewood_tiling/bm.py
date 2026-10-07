from fractions import Fraction
def bm(s):
    s = [Fraction(x) for x in s]
    C = [Fraction(1)]; B = [Fraction(1)]; L = 0; m = 1; b = Fraction(1)
    for n in range(len(s)):
        d = s[n] + sum(C[i] * s[n - i] for i in range(1, L + 1))
        if d == 0:
            m += 1
        elif 2 * L <= n:
            T = C[:]; coef = d / b
            C = C + [Fraction(0)] * max(0, len(B) + m - len(C))
            for i in range(len(B)): C[i + m] -= coef * B[i]
            L = n + 1 - L; B = T; b = d; m = 1
        else:
            coef = d / b
            C = C + [Fraction(0)] * max(0, len(B) + m - len(C))
            for i in range(len(B)): C[i + m] -= coef * B[i]
            m += 1
    return L, C[:L + 1]
if __name__ == "__main__":
    a = [2, 40, 544, 6912, 87552, 1094144, 13534208, 165978112, 2022215680, 24520433664, 296325611520, 3572784594944, 43009509490688, 517207712333824]
    for start in range(0, 4):
        L, C = bm(a[start:])
        print('start', start, 'order', L, 'from', len(a) - start, 'terms; coeffs', [str(c) for c in C])
