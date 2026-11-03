from math import gcd, isqrt
x, y = (9, 2)
for j in range(8):
    if not (x * x - 20 * y * y == 1 and x % 10 == 9):
        raise AssertionError()
    t = (x + 1) // 10
    a, b, c = (4 * t - 1, 2 * y, 6 * t - 1)
    if not (a * a + b * b == c * c and 3 * a + 1 == 2 * c and (gcd(a, c) == 1) and (a % 6 == 3) and (c % 6 == 5)):
        raise AssertionError()
    x, y = (161 * x + 720 * y, 36 * x + 161 * y)
x, y = (381, 144)
for j in range(8):
    if not (x * x - 7 * y * y == 9 and (x + 4) % 7 == 0 and (y % 3 == 0)):
        raise AssertionError()
    a, b = ((x + 4) // 7, y // 3)
    if not (4 * a - 1) % 3 == 0:
        raise AssertionError()
    c = (4 * a - 1) // 3
    if not (a * a + b * b == c * c and 3 * c + 1 == 4 * a and (gcd(a, c) == 1) and (a % 2 == c % 2 == 1)):
        raise AssertionError()
    x, y = (127 * x + 336 * y, 48 * x + 127 * y)
f = [0, 1]
for i in range(100):
    f.append(f[-1] + f[-2])
for n in range(2, 50):
    if n % 3 == 1:
        continue
    a, b, c = (f[n + 1] ** 2 - f[n] ** 2, 2 * f[n] * f[n + 1], f[2 * n + 1])
    if not (a * a + b * b == c * c and gcd(a, c) == 1 and (c % 6 == (1 if n % 4 in (0, 3) else 5))):
        raise AssertionError()
for a in range(1, 10000, 2):
    for q in range(1, 15):
        z = 3 * a + 1
        if z % (1 << q):
            continue
        c = z // (1 << q)
        if c <= a or c % 2 == 0:
            continue
        b = isqrt(c * c - a * a)
        if b * b == c * c - a * a:
            if not (a, b, c) in [(3, 4, 5), (1155, 1292, 1733)]:
                raise AssertionError()
print('All checks passed: 8 forward, 8 reverse, 32 Fibonacci examples; exhaustive forward a<10000.')
