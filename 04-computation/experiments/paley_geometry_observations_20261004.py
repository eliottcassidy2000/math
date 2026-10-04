"""Exact Paley projective transport and certified finite-observation aliases.

Standard library; reuses the declared inverse-ray certificate implementation.
Every check is active under -O. No conjectural convergence is assumed.
"""
from itertools import permutations, product
from math import gcd
from functools import lru_cache

from inverse_ray_ternary_addresses_20261004 import (
    ROOT, extend, kappa, mod3, mod2, exponent, chain, ranks, bit_bounds, expand,
)

CHECKS = 0


def need(condition, label):
    global CHECKS
    CHECKS += 1
    if not condition:
        raise ArithmeticError(label)


def chi(x, p):
    x %= p
    if x == 0:
        return 0
    return 1 if pow(x, (p-1)//2, p) == 1 else -1


def vec(x, p):
    return (1, 0) if x == p else (x, 1)


def sign(x, y, p):
    a, b = vec(x, p)
    c, d = vec(y, p)
    return chi(a*d-b*c, p)


def determinant(M, p):
    a, b, c, d = M
    return (a*d-b*c) % p


def projective_action(M, x, p):
    a, b, c, d = M
    u, v = vec(x, p)
    r, s = (a*u+b*v) % p, (c*u+d*v) % p
    return ((r*pow(s, -1, p)) % p, s) if s else (p, r)


def projective_matrices(p):
    output = set()
    for M in product(range(p), repeat=4):
        if determinant(M, p):
            first = next(x for x in M if x)
            output.add(tuple(x*pow(first, -1, p) % p for x in M))
    return sorted(output)


def compose(a, b):
    return tuple(a[b[i]] for i in range(len(a)))


def generated_permutations(generators):
    identity = tuple(range(len(generators[0])))
    seen, todo = {identity}, [identity]
    while todo:
        x = todo.pop()
        for g in generators:
            y = compose(g, x)
            if y not in seen:
                seen.add(y)
                todo.append(y)
    return seen


def perm_order(p):
    x, identity = p, tuple(range(len(p)))
    n = 1
    while x != identity:
        x = compose(p, x)
        n += 1
    return n


def projective_controls():
    p = 7
    S = [[sign(x, y, p) for y in range(p+1)] for x in range(p+1)]
    need(all(S[x][y] == -S[y][x] for x in range(8) for y in range(8)),
         "antisymmetric determinant character")
    projective = projective_matrices(p)
    need(len(projective) == 336, "PGL2(7) size")
    square_permutations = set()
    for M in projective:
        delta = chi(determinant(M, p), p)
        values = [projective_action(M, x, p) for x in range(8)]
        perm = tuple(v[0] for v in values)
        eps = [chi(v[1], p) for v in values]
        for x in range(8):
            for y in range(x):
                need(S[perm[x]][perm[y]] == delta*eps[x]*eps[y]*S[x][y],
                     "projective switching cocycle")
        if delta == 1:
            square_permutations.add(perm)
    counts = {1: 0, -1: 0}
    strict = 0
    for perm in permutations(range(8)):
        if all(S[perm[x]][perm[y]] == S[x][y] for x in range(8) for y in range(x)):
            strict += 1
        for delta in (1, -1):
            eps = [1] + [delta*S[perm[0]][perm[x]]*S[0][x] for x in range(1, 8)]
            if all(S[perm[x]][perm[y]] == delta*eps[x]*eps[y]*S[x][y]
                   for x in range(8) for y in range(x)):
                counts[delta] += 1
                if delta == 1:
                    need(perm in square_permutations, "no extra switching permutations")
    need(strict == 21 and counts == {1: 168, -1: 168}, "symmetry census")
    matrices = [(0, -1, 1, 0), (0, -1, 1, 1)]
    generators = [tuple(projective_action(M, x, 7)[0] for x in range(8))
                  for M in matrices]
    need([perm_order(g) for g in generators] == [2, 3], "triangle generators")
    need(perm_order(compose(*generators)) == 7, "triangle product order7")
    need(generated_permutations(generators) == square_permutations,
         "triangle generators generate PSL2(7)")
    # The signed lift has genuine ties above each projective point.
    lifts = list(product(range(8), (-1, 1)))
    tie_pairs = sum(x == y for i, (x, e) in enumerate(lifts) for y, f in lifts[:i])
    need(tie_pairs == 8, "eight antipodal ties in signed cover")
    sl_count = 0
    for M in product(range(7), repeat=4):
        if determinant(M, 7) != 1:
            continue
        sl_count += 1
        images = [(projective_action(M, x, 7)[0],
                   e*chi(projective_action(M, x, 7)[1], 7)) for x, e in lifts]
        for i, (x, e) in enumerate(lifts):
            for j, (y, f) in enumerate(lifts[:i]):
                xx, ee = images[i]
                yy, ff = images[j]
                need(ee*ff*S[xx][yy] == e*f*S[x][y], "SL2 signed-cover observable")
    need(sl_count == 336, "SL2 signed-cover group size")
    need(all(projective_action((-1, 0, 0, -1), x, 7) == (x, 6)
             for x in range(7)) and projective_action((-1, 0, 0, -1), 7, 7) == (7, 6),
         "central minus-I changes every lift sheet")
    need(all(sign(x,y,5) == sign(y,x,5) for x in range(6) for y in range(6)),
         "characteristic5 is symmetric: not a tournament")
    print("Projective7: strict21, switching168, anti-switching168; SL2 signed cover336.")
    print("Signed cover:16 vertices,112 oriented pairs,8 antipodal ties.")


def word_data(word):
    P, Q, B = 1, 1, 0
    for a in word:
        P, Q, B = 3*P, Q*2**a, 3*B+Q
    return P, Q, B


def paley_word_controls():
    primes = (7, 11, 19, 23, 31, 43, 47)
    words = [w for length in range(1,5) for w in product(range(1,5), repeat=length)]
    for p in primes:
        for w in words:
            P, Q, B = word_data(w)
            delta = chi(3,p)**len(w)*chi(2,p)**sum(w)
            invQ = pow(Q,-1,p)
            images = [(P*x+B)*invQ % p for x in range(p)]
            for x in range(p):
                for y in range(x):
                    need(chi(images[x]-images[y],p) == delta*chi(x-y,p),
                         "word action on Paley arcs")
    for p in primes:
        P,Q,B = word_data((1,))
        R,T,C = word_data((3,))
        need(chi(P*Q,p) == chi(R*T,p) and P>Q and R<T,
             "same quadratic orientation, opposite real drift")
    need(word_data((1,2)) == (9,8,5) and word_data((2,1)) == (9,8,7),
         "equal clocks can lose carry and axis")
    print("Paley word census:",len(words),"words at",len(primes),"primes.")
    print("p mod24 -> orientation:7=(-1)^j;11=(-1)^A;19=(-1)^(j+A);23=+1.")
    print("Hostiles: words1 and3 have opposite drift and every quadratic signature equal.")


def phi(n):
    result, divisor = n, 2
    while divisor*divisor <= n:
        if n % divisor == 0:
            result -= result//divisor
            while n % divisor == 0:
                n //= divisor
        divisor += 1
    if n > 1:
        result -= result//n
    return result


def split_modulus(M):
    if type(M) is not int or M < 1:
        raise ValueError("M must be a positive integer")
    h = k = 0
    m = M
    while m % 2 == 0:
        h += 1
        m //= 2
    while m % 3 == 0:
        k += 1
        m //= 3
    return h, k, m


@lru_cache(None)
def coprime3_residue(cert, modulus):
    if cert.parent is None:
        return 1 % modulus
    return ((pow(2,exponent(cert),modulus)*coprime3_residue(cert.parent,modulus)-1)
            * pow(3,-1,modulus)) % modulus


def residue(cert, M):
    if M == 1:
        return 0
    h,k,m = split_modulus(M)
    pairs = []
    if h: pairs.append((mod2(cert,h),2**h))
    if k: pairs.append((mod3(cert,k),3**k))
    if m>1: pairs.append((coprime3_residue(cert,m),m))
    x, modulus = 0, 1
    for y, n in pairs:
        x += modulus*((y-x)*pow(modulus,-1,n) % n)
        modulus *= n
    need(modulus == M and 0<=x<M, "CRT reconstruction")
    return x


def inverse_exponent(parent, a):
    target = mod3(parent,2)
    numerator = (pow(2,a,9)*target-1) % 9
    need(numerator % 3 == 0, "chosen exponent integrality")
    row = numerator//3
    first = kappa(target,row)
    need(a>=first and (a-first)%6==0, "chosen exponent channel")
    return extend(parent,row,(a-first)//6)


def alias_pair(M,D):
    if type(D) is not int or D < 0:
        raise ValueError("D must be a nonnegative integer")
    h,k,m = split_modulus(M)
    j = D+h+1
    L = 4*j*3**(j+k-1)*phi(m)
    hub = inverse_exponent(ROOT,3**j+1)
    small = inverse_exponent(hub,1)
    large = inverse_exponent(hub,1+L)
    for _ in range(j-1):
        small = inverse_exponent(small,1)
        large = inverse_exponent(large,1)
    return small,large,hub,j,L


def alias_controls():
    expanded = 0
    for M in range(1,49):
        for D in range(5):
            small,large,hub,j,L = alias_pair(M,D)
            s_nodes,l_nodes = chain(small),chain(large)
            need(ranks(small) == (j+1,2*j+3**j+2), "small first-hit ranks")
            need(ranks(large) == (j+1,2*j+3**j+2+L), "large first-hit ranks")
            for i in range(D+1):
                need(residue(s_nodes[i],M)==residue(l_nodes[i],M), "aliased observed vertex")
            for i in range(D):
                need(exponent(s_nodes[i])==exponent(l_nodes[i])==1, "aliased observed edge")
            if bit_bounds(large)[1] <= 4096:
                n,np,u = expand(small,4096),expand(large,4096),expand(hub,4096)
                need(n<u<np, "actual opposite net drift")
                need(n==2**j*(u+1)//3**j-1, "small closed formula")
                need(np==2**j*(2**L*u+1)//3**j-1, "large closed formula")
                for source,word in ((n,[1]*j),(np,[1]*(j-1)+[1+L])):
                    value=source
                    for a in word:
                        v=3*value+1
                        actual=(v & -v).bit_length()-1
                        need(actual==a, "independent literal exact valuation")
                        value=v>>actual
                    need(value==u and (3*u+1)==2**(3**j+1), "literal completed suffix")
                for i in range(D+1):
                    need(residue(s_nodes[i],M)==expand(s_nodes[i],4096)%M,
                         "independent modular small source")
                    need(residue(l_nodes[i],M)==expand(l_nodes[i],4096)%M,
                         "independent modular large source")
                expanded+=1
    M,D = 2**64*3**12*5*7,10
    small,large,hub,j,L=alias_pair(M,D)
    s_nodes,l_nodes=chain(small),chain(large)
    for i in range(D+1):
        need(residue(s_nodes[i],M)==residue(l_nodes[i],M), "huge symbolic alias")
    for i in range(D):
        need(exponent(s_nodes[i])==exponent(l_nodes[i])==1, "huge observed edge")
    need(bit_bounds(small)[0]>10**35 and bit_bounds(large)[0]>10**35,
         "huge sources intentionally unexpanded")
    print("Certified aliases:240 modulus/depth cases;",expanded,"literal expanded pairs.")
    print("Huge control: M=2^64*3^12*35,D=10,j=75; both sources exceed10^35 binary digits.")
    small,large,hub,j,L=alias_pair(1,0)
    print("Smallest declared control:",expand(small),expand(large),"->",expand(hub),"->1.")
    print("Boundary: fixed finite observations cannot determine later net drift; adaptive tests remain available.")


if __name__ == '__main__':
    projective_controls()
    paley_word_controls()
    alias_controls()
    print("PASS:",CHECKS,"exact checks; all guards active under -O.")
