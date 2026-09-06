#!/usr/bin/env python3
"""Independent small-universe audit of the trinomial classification.

Enumerates positive counts (y,z) and derives the negative count, rather than
using the producer's fixed-mass negative-count scan or modular inverse.
No producer imports, symbolic algebra, floating point, or Python asserts.
"""
from math import gcd
import sys

sys.stdout.reconfigure(newline="\n")
CHECKS = 0


def require(condition, label):
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError(label)


def channels(a, b, c, max_mass):
    rows = {}
    for y in range(max_mass+1):
        for z in range(max_mass-y+1):
            positive = b*y+c*z
            if positive == 0 or positive % a:
                continue
            x = positive//a
            mass = x+y+z
            if mass <= max_mass:
                rows.setdefault(mass, set()).add((x,y,z))
    return rows


def member(n, A, B):
    if n < 0:
        return False
    return any((n-B*z) % A == 0 for z in range(n//B+1))


def main():
    supports = collided = 0
    for a in range(1,13):
        for b in range(1,13):
            for c in range(b+1,13):
                if gcd(a,gcd(b,c)) != 1:
                    continue
                supports += 1
                g = gcd(a+b,a+c)
                A, B = (a+b)//g, (a+c)//g
                bound = min((a+b)//gcd(a,b),(a+c)//gcd(a,c))
                rows = channels(a,b,c,bound)
                m0 = min(rows)
                collision = len(rows[m0]) >= 2
                predicted = member(a-A*B,A,B)
                require(collision == predicted, ("classification",a,b,c))
                if not collision:
                    continue
                collided += 1
                require(m0 == g and g > B, "first mass and width")
                z0 = next(z for z in range(A) if (a-B*z) % A == 0)
                y0 = (a-B*z0)//A
                h, r = divmod(y0,B)
                higher = channels(a,b,c,3*g)
                first = higher[g]
                generated = {(0,0,0)}
                for k in range(1,4):
                    generated = {tuple(v[j]+w[j] for j in range(3))
                                 for v in generated for w in first}
                    predicted_count = k*h+(k*r)//B+(k*z0)//A+1
                    require(len(higher[k*g]) == predicted_count, "carry count")
                    require(generated <= higher[k*g], "addition preserves channel")
                    require(len(generated) == k*h+1, "generated first-slice count")
                    require(len(higher[k*g]-generated) == (k*r)//B+(k*z0)//A,
                            "exact missing channel count")
    relation_cases = 0
    for A,B,g in ((1,2,5),(2,3,7),(3,4,13),(3,5,16),(4,5,21),(5,7,36)):
        require(g > 2*B and gcd(A,B) == gcd(g,A*B) == 1, "family hypotheses")
        a,b,c = A*B,A*(g-B),B*(g-A)
        min_norm = None
        for R in range(-2*B,2*B+1):
            for S in range(-2*B,2*B+1):
                numerator = a*R-b*S
                if numerator % c:
                    continue
                T = numerator//c
                norm = abs(R)+abs(S)+abs(T)
                if norm:
                    min_norm = norm if min_norm is None else min(min_norm,norm)
        require(min_norm == 2*B, "exact bounded shortest signed relation")
        rows = channels(a,b,c,g)
        require(min(rows) == g and len(rows[g]) == 2, "family two first channels")
        relation_cases += 1
    print("PASS: independent positive-count enumeration of trinomial classification")
    print("primitive supports:", supports, "; collided:", collided)
    print("all-level carry tests: k=1,2,3 on every collided support")
    print("shortest-relation family cases:", relation_cases)
    print("exact_checks:", CHECKS)
    print("Scope: supports with 1<=a<=12,1<=b<c<=12; six relation controls.")


if __name__ == "__main__":
    main()
