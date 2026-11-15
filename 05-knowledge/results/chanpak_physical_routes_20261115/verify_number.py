"""Finite controls for the two-slope affine orbit theorem."""
from itertools import product
from math import gcd

def orbit(x, maps, m):
    seen = {x % m}
    todo = list(seen)
    while todo:
        y = todo.pop()
        for f in maps:
            z = f(y) % m
            if z not in seen:
                todo.append(z)
                seen.add(z)
    return seen

checks = 0
for m in range(1, 61):
    if gcd(m, 6) != 1:
        continue
    for a, b, x in product(range(-3, 4), repeat=3):
        d = gcd(m, 2*a-b)
        want = orbit(x+a, [lambda y: 2*y, lambda y: 3*y], d)
        expected = {y for y in range(m) if (y+a) % d in want}
        got = orbit(x, [lambda y: 2*y+a, lambda y: 3*y+b], m)
        if got != expected:
            raise RuntimeError((m, a, b, x, d, expected, got))
        checks += 1
print(f"verified {checks} affine orbit controls, including nonunit starts")
