"""Exact integer certificates after choosing test vectors by rounded cosines."""
from itertools import product
from math import cos, pi

MOD = 125
TURN = 57


def orbit(x, n):
    return tuple(sorted({x*pow(TURN, j, n) % n for j in range(4)}))


def test_transversals():
    representatives = sorted({orbit(i,25)[0] for i in range(1,25)})
    frequencies = sorted({orbit(i,MOD)[0] for i in range(1,MOD) if i%5})
    assert len(representatives)==6 and len(frequencies)==25
    orbit_options = [[orbit(i+25*t,MOD) for t in range(5)] for i in representatives]
    table, norms = [], []
    for k in frequencies:
        # Only the integer Rayleigh inequality below certifies a rejection.
        rounded = [round(10000*cos(2*pi*(k*i%MOD)/MOD)) for i in range(MOD)]
        s = sum(rounded)
        v = [MOD*x-s for x in rounded]
        assert sum(v)==0
        norms.append(sum(x*x for x in v))
        corr = [sum(v[j]*v[(j+d)%MOD] for j in range(MOD)) for d in range(MOD)]
        table.append([[sum(corr[d] for d in o) for o in opts] for opts in orbit_options])
    ncert = 0
    for choices in product(range(5), repeat=6):
        for n, offsets in zip(norms, table):
            dot = n + sum(offsets[i][j] for i,j in enumerate(choices))
            if 5*abs(dot) > 39*n:
                ncert += 1
                break
        else:
            raise AssertionError(choices)
    assert ncert == 15625
    return ncert


def mod101():
    q, b = 101, 16
    assert pow(b,25,q)==1 and pow(b,5,q)!=1
    hits = []
    for tail in product(range(5), repeat=4):
        exponents = tuple(r+5*j for r,j in enumerate((0,)+tail))
        if sum(pow(b,r,q) for r in exponents)%q==0:
            hits.append(sum(pow(b,7*r,q) for r in exponents)%q)
    assert hits == [62,83,56,74,19,73,95,15,80,31]
    return hits


if __name__ == "__main__":
    print("Affine quarter-turn transversal rejections at 39/5:", test_transversals())
    print("Surviving second evaluations modulo 101:", mod101())
