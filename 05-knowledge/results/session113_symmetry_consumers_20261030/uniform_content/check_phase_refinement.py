"""A proved finite-epoch upper bound using low multiplier precision."""
from fractions import Fraction as F
from check_obstruction import W, U, v5

def epoch_bound(m, h, L, r):
    if m.numerator * pow(m.denominator, -1, 5) % 5 not in (1, 4):
        return F(0)
    total = F(0)
    d = 0
    for e in range(1, min(h, r) + 1):
        q = 5 ** e
        T = 2 * 5 ** (e - 1)
        x = m.numerator * pow(m.denominator, -1, q) % q
        for d in range(T):
            if x == 1:
                break
            x = x * 3 * pow(2, -1, q) % q
        else:
            raise AssertionError('square multiplier not in generated subgroup')
        total += max(0, 1 + (L - d) // T)
    if h > r and L >= d:
        total += h - r + F(5, 8) * (F(1, 5 ** r) - F(1, 5 ** h)) * (L - d)
    return total

def check(word, r):
    P = Q = 1
    B = 0
    epochs = []
    for a in word:
        P, B, Q = (3 * P, 3 * B + Q, Q * 2 ** a)
        if a > 1:
            epochs.append([F(P, Q), v5(B + Q - P), 0])
        else:
            epochs[-1][2] += 1
    bounds = [epoch_bound(*e, r) for e in epochs]
    return {'r': r, 'bound': str(sum(bounds)), 'excess': str(sum(bounds) - F(5, 8) * sum(word)), 'epoch_bounds': list(map(str, bounds))}
if __name__ == '__main__':
    for r in (1, 2):
        result = check(U, r)
        if not F(result['excess']) < 0:
            raise AssertionError('verification failed')
        print('U5', result)
    V = [3, 1, 3] + [1] * 6 + [3, 1, 3]
    U3 = V + [1] * 231 + V
    if not (len(U3) == 255 and sum(U3) == 271):
        raise AssertionError('verification failed')
    for r in (1, 2):
        result = check(U3, r)
        if not F(result['excess']) < 0:
            raise AssertionError('verification failed')
        print('U3', result)
    from math import gcd
    checked = 0

    def sweep(word, remaining):
        global checked
        if word:
            P = Q = 1
            B = 0
            actual = 0
            for a in word:
                P, B, Q = (3 * P, 3 * B + Q, Q * 2 ** a)
                actual += v5(gcd(abs(P - Q), B))
            for r in (1, 2):
                if not actual <= F(check(word, r)['bound']):
                    raise AssertionError((word, r, actual))
            checked += 1
        for a in range(1 if word else 2, remaining + 1):
            sweep(word + [a], remaining - a)
    sweep([], 12)
    print('finite upper-envelope checks:', checked, 'words at each r=1,2')
