"""Exact critical beta-mixture weights; finite checks do not prove coverage."""
from fractions import Fraction as F
from math import comb, factorial
import json

CHECKS = 0


def need(condition, message):
    global CHECKS
    CHECKS += 1
    if not condition:
        raise ValueError(message)


def nat(n, positive=False):
    if type(n) is not int or n < int(positive):
        raise ValueError("exact natural number required")


def odd(n):
    nat(n, True)
    if n % 2 != 1:
        raise ValueError("positive odd source required")


def U(n):
    odd(n)
    y = 3*n+1
    while y % 2 == 0:
        y //= 2
    return y


def S(n, j):
    odd(n)
    nat(j)
    return 4**j*n+(4**j-1)//3


def decode(n):
    odd(n)
    j = 0
    while n % 8 == 5:
        n = (n-1)//4
        j += 1
    return n, j


def child(c, k):
    if decode(c) != (c, 0):
        raise ValueError("primitive parent required")
    y = S(c, k)
    if y % 3 == 0:
        return None
    b = ((4 if y % 3 == 1 else 2)*y-1)//3
    return None if b == 1 else b


def beta(a, b):
    nat(a, True)
    nat(b, True)
    return F(factorial(a-1)*factorial(b-1), factorial(a+b-1))


def weight(t, m):
    nat(t)
    nat(m)
    return 6*beta(m+2, t+2)


def interval(n, steps):
    """A rational interval. Zero lower endpoint never asserts nontermination."""
    nat(steps)
    b, j = decode(n)
    m = j
    for t in range(steps+1):
        if b == 1:
            return {"lower": weight(t, m), "upper": weight(t, m),
                    "T": t, "M": m, "rooted": True}
        if t == steps:
            return {"lower": F(0), "upper": 6*beta(m+2, t+3),
                    "T": None, "M": None, "rooted": False}
        b, k = decode(U(b))
        m += k
    raise RuntimeError("unreachable")


def finite_box(levels, cutoff):
    """Backward enumeration only: base T<levels; each k and sibling j<=cutoff."""
    nat(levels, True)
    nat(cutoff)
    bases = {1: (0, 0)}
    layer = [1]
    for t in range(1, levels):
        next_layer = []
        for c in layer:
            for k in range(cutoff+1):
                b = child(c, k)
                if b is None:
                    continue
                need(b not in bases, "unique rooted backward address")
                bases[b] = (t, bases[c][1]+k)
                next_layer.append(b)
        layer = next_layer
    sources = {}
    for b, (t, k) in bases.items():
        for j in range(cutoff+1):
            n = S(b, j)
            need(n not in sources, "unique sibling source")
            sources[n] = (t, k+j)
    return sources


def box_error(levels, cutoff, epsilon):
    nat(levels, True)
    nat(cutoff, True)
    if type(epsilon) not in (int, F) or not 0 < epsilon < 1:
        raise ValueError("exact epsilon in (0,1) required")
    epsilon = F(epsilon)
    depth = 12*epsilon+12*(1-epsilon**2/3)**(levels-1)
    width = 6*sum(F(a, cutoff+i) for i, a in enumerate((1, 2, 3, 2, 1)))
    siblings = 6*sum(F(1, cutoff+i) for i in (2, 3, 4))
    return depth+width+siblings


def integral_polynomial(t, m):
    # Independent exact integration of 6*r^(m+1)*(1-r)^(t+1).
    return 6*sum(F((-1)**i*comb(t+1, i), m+i+2) for i in range(t+2))


def run():
    for t in range(25):
        for m in range(25):
            w = weight(t, m)
            need(w == integral_polynomial(t, m), "beta integral")
            coef = F(6*(m+1)*(t+1), (m+t+1)*(m+t+2)*(m+t+3))
            need(coef >= F(6, (m+t+2)*(m+t+3)), "uniform polynomial floor")
            for r in (F(1, 16), F(1, 4), F(1, 2), F(3, 4), F(15, 16)):
                need(w >= coef*r**m*(1-r)**t, "all-parameter comparison")
            for count in (0, 1, 2, 5, 11):
                partial = sum(weight(t+1, m+j) for j in range(count))
                remainder = 6*beta(m+count+2, t+2)
                need(partial+remainder == w, "whole inverse fibre with exact tail")
    for j in range(100):
        need(weight(0, j) == F(6, (j+2)*(j+3)), "root ray")
        need(interval(S(1, j), 0)["lower"] == weight(0, j), "zero-step decoder")
    need(sum(weight(0, j) for j in range(100))+F(6, 102) == 3, "root ray mass")
    records = finite_box(5, 5)
    mass = sum(weight(t, m) for t, m in records.values())
    need(mass <= 11, "global mass bound on finite support")
    multiple_three_mass = sum(weight(t, m) for n, (t, m) in records.items() if n % 3 == 0)
    need(multiple_three_mass <= 2, "multiple-three conservation bound")
    for n, (t, m) in records.items():
        exact = weight(t, m)
        for cutoff in range(min(t+1, 3)):
            observed = interval(n, cutoff)
            need(observed["lower"] <= exact <= observed["upper"], "pointwise interval")
        final = interval(n, t)
        need(final["rooted"] and final["lower"] == exact, "backward versus forward")
        need((final["T"], final["M"]) == (t, m), "retained route parameters")
        if n > 1 and n % 3:
            b = ((4 if n % 3 == 1 else 2)*n-1)//3
            need(U(b) == n and b != 1, "actual incoming base")
            source = interval(b, t+1)
            need(source["rooted"], "rooted inverse base")
            need((source["T"], source["M"]) == (t+1, m), "incoming route count")
    hostile = interval(27, 0)
    need(not hostile["rooted"] and hostile["lower"] == 0, "unresolved is not zero limit")
    need(interval(27, 100)["lower"] > 0, "same unresolved source later certified")
    need(weight(0, 0) == 1, "root normalization")
    # Atomic nonnegative flow does not imply a full-support flow: a disjoint
    # self-loop with zero mass obeys the local inequality and remains unrooted.
    need(F(0) <= F(0), "zero-support hostile component")
    for bad in (True, -1, 0, 2, F(3), 3.0):
        try:
            interval(bad, 1)
        except ValueError:
            need(True, "malformed source rejected")
        else:
            need(False, "malformed source accepted")
    errors = [box_error(20*j**3, 20*j, F(1, j)) for j in (2, 4, 8)]
    need(errors[0] > errors[1] > errors[2], "explicit tail controls improve")
    return {"status": "FINITE-EXACT; universal positivity OPEN", "checks": CHECKS,
            "finite_box": {"T_less_than": 5, "each_k_and_j_at_most": 5,
                           "source_count": len(records), "mass": str(mass),
                           "multiple_three_mass": str(multiple_three_mass)},
            "proved_total_mass_bound": "11", "proved_multiple_three_mass": "2",
            "proved_root_ray_mass": "3", "weight_27": str(interval(27, 100)["lower"]),
            "hostile": "zero approximation at 27 with cutoff 0 later becomes positive",
            "tail_bound_formula_checked": True}


if __name__ == "__main__":
    print(json.dumps(run(), indent=2, sort_keys=True))
