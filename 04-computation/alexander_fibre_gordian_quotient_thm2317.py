#!/usr/bin/env python3
"""Exact companion checks for THM-2317.

This is deliberately small: it checks the load-bearing Fox polynomial
identities, the two mod-2 maximal-ideal controls, and the three-body word
metric hostile example.  It does not enumerate knots or assert an unknotting
number.
"""

from heapq import heappop, heappush


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def mul(a: tuple[int, ...], b: tuple[int, ...]) -> tuple[int, ...]:
    out = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i + j] += x * y
    return tuple(out)


def add3(x: tuple[int, int, int], y: tuple[int, int, int]) -> tuple[int, int, int]:
    return (x[0] + y[0], x[1] + y[1], x[2] + y[2])


def eval_mod2(poly: tuple[int, ...], value: int) -> int:
    return sum(coefficient * (value**degree) for degree, coefficient in enumerate(poly)) % 2


a = (1, 1, 1, 1, 1, 1, 1)
b = (1, 1)
delta = (1, -1, 1, -1, 1, -1, 1)
one_plus_t7 = (1, 0, 0, 0, 0, 0, 0, 1)
even_sum = (1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1)

require(mul(b, delta) == one_plus_t7, "b*Delta factorization failed")
require(mul(a, delta) == even_sum, "a*Delta factorization failed")
require(sum(((-1) ** j) for j in range(7)) == 1, "a(-1) is not one")

f = (1, 1, 0, 1)  # t^3+t+1, ascending coefficients
g = (1, 0, 1, 1)  # t^3+t^2+1
require(
    tuple(c % 2 for c in mul(f, g)) == tuple(c % 2 for c in delta),
    "mod-two factorization failed",
)
require(
    eval_mod2(f, 0) == eval_mod2(f, 1) == 1,
    "the cubic q1 factor has a root over F_2",
)
require(sum(delta) % 2 == 1, "Delta vanishes at the q0 control")

generators = (
    ((1, 0, 0), 1),
    ((-1, 0, 0), 1),
    ((0, 1, 0), 1),
    ((0, -1, 0), 1),
    ((0, 0, 1), 1),
    ((0, 0, -1), 1),
    ((1, 1, 1), 2),
    ((-1, -1, -1), 2),
)


def word_dist(target: tuple[int, int, int], radius: int = 4) -> int:
    origin = (0, 0, 0)
    distance = {origin: 0}
    queue = [(0, origin)]
    while queue:
        cost, state = heappop(queue)
        if distance[state] != cost:
            continue
        if state == target:
            return cost
        for step, weight in generators:
            nxt = add3(state, step)
            if max(map(abs, nxt)) > radius:
                continue
            new_cost = cost + weight
            if new_cost < distance.get(nxt, 10**9):
                distance[nxt] = new_cost
                heappush(queue, (new_cost, nxt))
    raise RuntimeError("target escaped the finite hostile-control box")


singletons = [word_dist((1, 0, 0)), word_dist((0, 1, 0)), word_dist((0, 0, 1))]
pairs = [word_dist((1, 1, 0)), word_dist((1, 0, 1)), word_dist((0, 1, 1))]
triple = word_dist((1, 1, 1))
require(singletons == [1, 1, 1], "wrong singleton word lengths")
require(pairs == [2, 2, 2], "wrong pair word lengths")
require(triple == 2, "wrong triple word length")

print("FOX b*Delta=1+t^7: exact")
print("FOX a*Delta=sum_(j=0)^6 t^(2j): exact")
print("MOD2 Delta=(t^3+t+1)(t^3+t^2+1): exact")
print("FIBRE controls: q1 sees Delta; q0=(2,t+1) does not")
print("WORD singleton lengths:", singletons)
print("WORD pair lengths:", pairs)
print("WORD triple length:", triple)
print("PAIR SHADOW: all zero defects; triple defect: 1")
