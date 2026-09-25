"""Guarded Collatz jumps on the inherited Berggren spine, with even states."""

from fractions import Fraction as F
from math import gcd, isqrt
from pathlib import Path
import importlib.util


def check(test, label):
    if not test:
        raise RuntimeError(label)


def v2(n):
    check(n > 0, "positive valuation input")
    return (n & -n).bit_length()-1


def shortcut(n, sigma=1):
    return (3*n+sigma)//2 if n % 2 else n//2


def accelerated(n, sigma=1):
    return (3*n+sigma)//2**v2(3*n+sigma)


def phi(n):
    check(n > 0 and n % 2, "odd spine domain")
    return n, (n*n-1)//2, (n*n+1)//2


def psi(n):
    d = 2 if n % 2 else 1
    return 2*n//d, (n*n-1)//d, (n*n+1)//d


def matvec(m, p):
    return tuple(sum(a*b for a,b in zip(row,p)) for row in m)


SPINE = ((1,-2,2),(2,-1,2),(2,-2,3))
SPINE_INV = ((1,2,-2),(-2,-1,2),(-2,-2,3))
MIDDLE = ((1,2,2),(2,1,2),(2,2,3))


def lorentz_lift(k, sigma):
    d = 2*4**k
    return tuple(tuple(F(v,d) for v in row) for row in (
        (3*2**(k+1), -sigma*2**(k+1), sigma*2**(k+1)),
        (6*sigma, 8+4**k, 10-4**k),
        (6*sigma, 8-4**k, 10+4**k)))


def main():
    for n in range(1,10002):
        a,b,c = psi(n)
        check(a*a+b*b == c*c and gcd(a,b) == 1, "all-integer primitive chart")
        check(F(a,c-b) == n, "marked chart inverse")
        check(c-b == (1 if n % 2 else 2), "parity gap")
    print('ALL-INTEGER MARKED CHART n1..10001 primitive/inverse/gap PASS')
    drops = rises = 0
    middle_checks = 0
    for n in range(1,10002,2):
        p = phi(n)
        check(matvec(SPINE,p) == phi(n+2), "inherited parabolic spine step")
        if n > 1:
            check(matvec(SPINE_INV,p) == phi(n-2), "spine inverse")
        child = matvec(MIDDLE,p)
        euclid_m, euclid_r = isqrt((child[2]+child[0])//2), isqrt((child[2]-child[0])//2)
        check((euclid_m,euclid_r) == (shortcut(n),(n+1)//2), "shortcut as marked middle-child parameter")
        check(euclid_m == 3*euclid_r-1 and gcd(euclid_m,euclid_r) == 1
              and (euclid_m+euclid_r) % 2 == 1, "exact marked-image condition")
        middle_checks += 1
        target = accelerated(n)
        j, jp = (n-1)//2, (target-1)//2
        if n % 4 == 1:
            check(4*jp <= 3*j, "guarded spine contraction")
            if n > 1:
                check(target < n and phi(target)[2] < p[2], "actual strict odd-return descent")
                drops += 1
        else:
            check(2*jp == 3*j+1 and target > n, "odd-depth hostile expansion")
            rises += 1
        k = v2(3*n+1)
        check(matvec(lorentz_lift(k,1),p) == phi(target), "guarded Lorentz similarity")
    print('SPINE/SHORTCUT controls all5001 odd states;',middle_checks,'middle-child projections;',drops,'drops;',rises,'rises PASS')

    metric = (-1,-1,1)
    for sigma in (-1,1):
        for k in range(1,13):
            m = lorentz_lift(k,sigma)
            for i in range(3):
                for j in range(3):
                    value = sum(metric[r]*m[r][i]*m[r][j] for r in range(3))
                    check(value == (F(9,4**k)*metric[i] if i == j else 0), "Lorentz similarity matrix identity")
    print('LORENTZ SIMILARITY both signs,k1..12: M^t J M=(9/4^k)J PASS')

    for a in range(1,2501):
        n,target = 8*a+1,6*a+1
        check(v2(3*n+1) == 2 and shortcut(shortcut(n)) == target,
              "exact two-shortcut macro")
        check((n-1)//2 == 4*a and (target-1)//2 == 3*a, "depth ratio3/4")
        check(phi(n)[2]-phi(target)[2] == 2*a*(7*a+1), "hypotenuse drop")
        check(n*n-target*target == 4*a*(7*a+1), "quadratic elliptic height drop")
    print('MACRO n=8a+1->6a+1, a1..2500: clock2/depth3/4/heightdrops PASS')

    route = tuple(psi(n) for n in (9,14,7))
    check(route == ((9,40,41),(28,195,197),(7,24,25)), "actual9-14-7 marked path")
    check(matvec(SPINE_INV,phi(9)) == phi(7), "one Berggren parent at9")
    check(matvec(MIDDLE,phi(9)) == (171,140,221), "marked middle child at9")
    check(psi(4)[2] > psi(5)[2], "primitive hypotenuse alone not global integer order")
    solutions = [y for y in range(1,10002,2) if 2*y+shortcut(y) == phi(y)[2]]
    check(solutions == [7], "14+11=25 unique positive odd state")
    check((abs(3*3-4*4),2*3*4,5*5) == phi(7), "inherited Gaussian square3-4-5")
    print('LOCAL RIGIDITY 2y+Tplus(y)=hypPhi(y) iff y7; Gaussian square(3,4,5)->(7,24,25) PASS')
    print('ACTUAL TRIANGLE ROUTE',route,'; accelerated edge is one inherited spine parent')
    tail = (7,11,17,13,5,1)
    check(all(accelerated(a) == b for a,b in zip(tail,tail[1:])), "verified7 root tail")
    for r in range(1,25):
        n = (7*4**r-1)//3
        check(3*n+1 == 7*4**r and v2(3*n+1) == 2*r, "exact common-target family")
        current = n
        for _ in range(2*r-1):
            current = shortcut(current)
        check(current == 14 and shortcut(current) == 7, "all family routes share14-7 junction")
        check(accelerated(n) == 7 and all(phi(x)[0] == x for x in (n,)+tail), "lifted root certificate")
    print('INFINITE CERTIFICATE SCHEMA n_r=(7*4^r-1)/3 ->7->11->17->13->5->1; exact controls r1..24 PASS')

    # Exact curve-group comparison with the separately proved coordinate decoder.
    source = Path(__file__).with_name('creation_elliptic_20260925.py')
    spec = importlib.util.spec_from_file_location('creation_elliptic_control',source)
    ec = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(ec)
    for a in range(1,5):
        q = ec.mul(a,ec.G)
        p = ec.add(ec.mul(8,q),ec.G)
        recovered = p
        digits = []
        for _ in range(3):
            recovered,r,t = ec.decode_step(recovered)
            digits.append((r,t))
        check(recovered == q and digits == [(1,0),(0,0),(0,0)], "creation prefix certifies macro cylinder")
        numerator = ec.add(ec.mul(3,p),ec.G)
        for _ in range(2):
            numerator = next(x for x in ec.halves(numerator) if abs(ec.alpha(x)) == 1)
        check(numerator == ec.add(ec.mul(6,q),ec.G), "elliptic actual two-step macro")
    print('ELLIPTIC GROUP macro a1..4: creation digits100 retain tailQ; actual endpoint6Q+G PASS')
    print('HOSTILE n3->5: same primitive spine but forward depth1->2 and hypotenuse5->13')
    print('SCOPE actual guarded local descent; no global Collatz/Berggren-tree equivalence')
    print('ALL CHECKS PASS')


if __name__ == '__main__':
    main()
