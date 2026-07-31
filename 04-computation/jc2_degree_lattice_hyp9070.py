"""JC(2) degree-pair geometry: the leading form is H^{n/g}, H^{m/g} with
deg H = g = gcd(n,m), and the Abhyankar-Moh reduction on (n,m) is the
EUCLIDEAN ALGORITHM.  So the continued fraction of n/m is the natural
complexity/alternation measure of a candidate counterexample, and the
metallic ratios [q;q,q,...] are its extremal rays.

This file establishes the elementary facts exactly (no literature needed):
 E1  For a Jacobian pair, Jac(P_n,Q_m)=0, so in TWO variables the leading
     forms are proportional powers of a common form H with deg H = g and
     P_n = c H^{n/g}, Q_m = c' H^{m/g}.   [g = gcd(n,m); a=n/g, b=m/g coprime]
 E2  n|m or m|n  =>  a=1 or b=1  =>  one leading form is (up to scalar) a
     power of the other, and the standard elementary reduction applies.
     So a counterexample needs a,b >= 2, coprime.
 E3  The Euclidean/continued-fraction length of a/b is the reduction depth.
     Among coprime pairs with max(a,b) <= N, depth is maximised exactly by
     consecutive Fibonacci pairs (golden), and among pairs whose partial
     quotients are all equal to q, by the q-metallic (Pell for q=2, ...).
"""
from math import gcd
from fractions import Fraction as Fr


def cf(a, b):
    """continued fraction of a/b (Euclidean quotients)."""
    out = []
    while b:
        out.append(a // b)
        a, b = b, a % b
    return out


def depth(a, b):
    return len(cf(a, b))


def metallic_convergents(q, N):
    """convergents of [q;q,q,...] = (q+sqrt(q^2+4))/2"""
    p0, p1 = 1, q
    q0, q1 = 0, 1
    out = []
    while p1 <= N:
        out.append((p1, q1))
        p0, p1 = p1, q * p1 + p0
        q0, q1 = q1, q * q1 + q0
    return out


def main():
    N = 60
    # E3: which coprime (a,b), a,b>=2, maximise Euclidean depth for given max?
    print("max Euclidean depth among coprime (a,b), 2<=a<b<=N, per N:")
    for NN in (8, 13, 21, 34, 55):
        best, arg = 0, None
        for b in range(3, NN + 1):
            for a in range(2, b):
                if gcd(a, b) != 1:
                    continue
                d = depth(a, b)
                if d > best:
                    best, arg = d, (a, b)
        print(f"   max b={NN:3d}: depth {best}  attained at {arg}  cf={cf(*arg)}")
    print("\nmetallic convergent pairs (a,b) up to 200:")
    for q in (1, 2, 3):
        name = {1: "golden/Fibonacci", 2: "silver/Pell", 3: "bronze"}[q]
        pairs = metallic_convergents(q, 200)
        cop = [(x, y) for (y, x) in zip([p for p, _ in pairs], [p for p, _ in pairs][1:])]
        seq = [p for p, _ in pairs]
        print(f"   q={q} ({name}): {seq}")
        print(f"        consecutive coprime ratios: "
              f"{[(seq[i], seq[i+1]) for i in range(len(seq)-1) if gcd(seq[i], seq[i+1]) == 1][:6]}")
    # the ratio n/m of a counterexample: a,b>=2 coprime.  Fibonacci pairs from
    # (2,3) on all qualify; Pell pairs (2,5),(5,12),(12,29),(29,70) qualify.
    print("\ncandidate degree ratios a/b with a,b>=2 coprime, by depth:")
    rows = []
    for b in range(3, 40):
        for a in range(2, b):
            if gcd(a, b) != 1:
                continue
            rows.append((depth(a, b), a, b, cf(a, b)))
    rows.sort(reverse=True)
    for d, a, b, c in rows[:12]:
        tag = ""
        if all(x == 1 for x in c[1:]) and c[0] == 1:
            tag = "  <- golden (Fibonacci)"
        if all(x == 2 for x in c[1:]):
            tag = "  <- silver (Pell)"
        print(f"   depth {d}: {a}/{b}  cf={c}{tag}")


if __name__ == "__main__":
    main()
