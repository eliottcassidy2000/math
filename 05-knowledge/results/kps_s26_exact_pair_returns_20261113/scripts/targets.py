"""Build the target files of the pre-registered searches S1, S2, S3 (THM-5151 PROOF section 7).

Run in the scratch directory that holds prof_n.txt (from fprof, step 2 of scripts/README.md) and c1p_12.txt
(from rprof on the C = 1 order-12 list, step 4).  Each output line is `c3 F_1 .. F_n`, where c3 is the triangle
count of the would-be tournament (half the coefficient of t^(n-2) in Q).

  S1: ytargets2.txt  - Q_X / t over the C = 1 order-12 classes X (would-be order-11 partners Y with Y -> TT1)
  S2: mtargets11.txt - Q(Y5 => Y7) / t over all odd strong Y5, Y7 (would-be order-11 odd merges Z)
  S3: targets10.txt  - all order-10 two-block sums TT1 + Y9, C3 => Y7, Y5 => Y5' (THM-5151 (1))

Searches: `ysearch 11 5 ytargets2.txt`, `ysearch 11 5 mtargets11.txt`, `yhash 10 5 targets10.txt`.
"""
from math import comb
from fractions import Fraction

def Qv(F, t):
    return sum(F[j - 1] * comb(t, j) for j in range(1, len(F) + 1))

def F_from_values(vals, n):
    """binomial-basis coefficients F_1..F_n from Q(0..n) (forward differences at 0)"""
    return tuple(sum((-1) ** (j - i) * comb(j, i) * vals[i] for i in range(j + 1)) for j in range(1, n + 1))

def c3_of(F):
    n = len(F)
    s = [[0] * (n + 2) for _ in range(n + 2)]; s[0][0] = 1          # signed Stirling numbers of the first kind
    for j in range(1, n + 1):
        for k in range(1, j + 1):
            s[j][k] = s[j - 1][k - 1] - (j - 1) * s[j - 1][k]
    fact = [1] * (n + 1)
    for i in range(1, n + 1):
        fact[i] = fact[i - 1] * i
    c = sum(Fraction(F[j - 1] * s[j][n - 2], fact[j]) for j in range(1, n + 1))
    assert c.denominator == 1
    return int(c) // 2

def strong_profiles(a):
    if a == 1:
        return [('TT1', (1,))]
    L = []
    for line in open('prof_%d.txt' % a):
        p = line.split()
        if int(p[2]) == 1:
            L.append((p[1], tuple(map(int, p[7:]))))
    return L

def write(name, profiles):
    with open(name, 'w') as f:
        for F in sorted(profiles):
            f.write(str(c3_of(F)) + ' ' + ' '.join(map(str, F)) + '\n')
    print(name, len(profiles))

def s1():
    T = set()
    for line in open('c1p_12.txt'):
        left, _ = line.split('|'); F = list(map(int, left.split()[2:]))
        n = len(F); vals = [0] + [Qv(F, t) // t for t in range(1, n + 1)]
        assert all(Qv(F, t) % t == 0 for t in range(1, n + 1))
        G = F_from_values(vals, n - 1)
        T.add(G)
    write('ytargets2.txt', T)

def s2():
    T = set(); n = 11
    for _, F1 in strong_profiles(5):
        for _, F2 in strong_profiles(7):
            vals = [0] + [Qv(F1, t) * Qv(F2, t) // t for t in range(1, n + 1)]
            T.add(F_from_values(vals, n))
    write('mtargets11.txt', T)

def s3():
    T = set(); n = 10
    for a, b in ((1, 9), (3, 7), (5, 5)):
        for _, F1 in strong_profiles(a):
            for _, F2 in strong_profiles(b):
                vals = [Qv(F1, t) * Qv(F2, t) for t in range(n + 1)]
                T.add(F_from_values(vals, n))
    write('targets10.txt', T)

if __name__ == '__main__':
    s1(); s2(); s3()          # expected counts: 1291, 655, 48908
