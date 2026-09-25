"""Independent intrinsic-triple and ternary-residue bridge audit.
See 05-knowledge/results/ternary_bridge_20260925.md. Standard library only.
"""
from math import gcd, isqrt, comb


def check(p, message="check failed"):
    if not p:
        raise RuntimeError(message)


def valuation(n, p):
    check(type(n) is int and n != 0 and p in (2,3))
    n, count = abs(n), 0
    while n % p == 0:
        n //= p
        count += 1
    return count


def validate(T):
    A,B,C = T
    check(all(type(v) is int for v in T), "integer triple required")
    check(A > 0 and C > 0 and A % 2 == 1 and B % 2 == 0,
          "odd leg/even signed leg convention")
    check(A*A+B*B == C*C and gcd(A,gcd(abs(B),C)) == 1,
          "primitive Pythagorean equation required")


def phi(s,t):
    check(type(s) is type(t) is int and s > 0 and t > 0 and s % 2 and t % 2)
    check(gcd(s,t) == 1, "primitive spinor required")
    T = (s*t, (s*s-t*t)//2, (s*s+t*t)//2)
    validate(T)
    return T


def decode(T):
    validate(T)
    A,B,C = T
    s,t = isqrt(C+B),isqrt(C-B)
    check(s*s == C+B and t*t == C-B)
    check(phi(s,t) == T)
    return s,t


def intrinsic_step(T, sign=1):
    validate(T)
    check(sign in (-1,1))
    A,B,C = T
    if sign == -1:
        check(4*C+5*B > 0, "minus map leaves positive rational domain")
    R = 10*C+8*B+6*sign*A
    k2 = valuation(R,2)
    check(k2 >= 2 and k2 % 2 == 0, "square carry valuation")
    k = k2//2
    K = 2**k
    g = 3 if (C-B) % 9 == 0 else 1
    nums = (3*A+sign*(C-B),
            6*sign*A+(8+K*K)*B+(10-K*K)*C,
            6*sign*A+(8-K*K)*B+(10+K*K)*C)
    dens = (K*g*g, 2*K*K*g*g, 2*K*K*g*g)
    check(all(n % d == 0 for n,d in zip(nums,dens)), "integer target guard")
    target = tuple(n//d for n,d in zip(nums,dens))
    validate(target)
    return target,k,g


def spinor_step(s,t,sign=1):
    z = 3*s+sign*t
    check(z > 0, "positive numerator required")
    k = 0
    while z % 2 == 0:
        z //= 2
        k += 1
    g = gcd(z,t)
    return (z//g,t//g),k,g


def repunit_mod(base,j,level):
    check(base > 1 and base % 3 == 1 and j >= 0 and level >= 0)
    q = 3**level
    return ((pow(base,j,q*(base-1))-1)//(base-1)) % q


def inverse_repunit_residue(base,h,level):
    j = 0
    for a in range(level):
        matches = [j+d*3**a for d in range(3)
                   if repunit_mod(base,j+d*3**a,a+1) == h % (3**(a+1))]
        check(len(matches) == 1, "unique ternary lift")
        j = matches[0]
    return j


def main():
    print("PROVED scoped bridge; FINITE-EXACT independent controls; Collatz remains OPEN")
    pair_checks = {1:0,-1:0}
    contractions = {1:0,-1:0}
    for s in range(1,256,2):
        for t in range(1,256,2):
            if gcd(s,t) != 1:
                continue
            T = phi(s,t)
            check(decode(T) == (s,t))
            for sign in (1,-1):
                if 3*s+sign*t <= 0:
                    continue
                out,k,g = intrinsic_step(T,sign)
                pair,kk,gg = spinor_step(s,t,sign)
                check((k,g) == (kk,gg) and out == phi(*pair))
                check(g == gcd(3,t))
                # The complete 3-free denominator is unchanged.
                q = t//3**valuation(t,3)
                tt = pair[1]
                check(tt//3**valuation(tt,3) == q)
                if t % 3 == 0:
                    check(out[2]*3 < T[2] if sign == 1 else out[2]*4 < T[2])
                    check(out[2]-out[1] == (T[2]-T[1])//9)
                    contractions[sign] += 1
                pair_checks[sign] += 1
    print("Intrinsic map vs independent spinor normalization:", pair_checks)
    print("Strict hypotenuse contraction checks off integer denominator level:", contractions)
    entry = 0
    for r in range(9):
        t = 3**r
        for s in range(1,1024,2):
            if gcd(s,t) != 1:
                continue
            T = phi(s,t)
            for j in range(r):
                T,k,g = intrinsic_step(T)
                check(g == 3 and T[2]-T[1] == 9**(r-j-1))
            check(decode(T)[1] == 1)
            entry += 1
    print("Exact r-step entry into unit-gap stratum: r0..8, odd s<=1023; cases=", entry)
    for seed in ((117,44,125),(39,80,89),(45,28,53)):
        x = intrinsic_step(seed)[0]
        y = intrinsic_step(x)[0]
        print("Marked triple trajectory:", seed, "->", x, "->", y)
    check(intrinsic_step((5,-12,13))[0] == (5,-12,13))
    check(intrinsic_step((5,12,13))[0] == (1,0,1))
    print("Orientation hostile: same unmarked5-12-13 triangle has fixed vs root-directed dynamics")
    for k in range(3,20):
        T = phi(1,2**k-3)
        check(intrinsic_step(T)[0] == T)
    print("Noninteger fixed point family1/(2^k-3) independently checked k3..19")
    minus = (783,56,785)
    minus1 = intrinsic_step(minus,-1)[0]
    minus2 = intrinsic_step(minus1,-1)[0]
    check(minus2 == (3,-4,5))
    try:
        intrinsic_step(minus2,-1)
        raise RuntimeError("minus zero boundary accepted")
    except RuntimeError as error:
        check(str(error) == "minus map leaves positive rational domain")
    print("Minus boundary:", minus,"->",minus1,"->",minus2,"-> rejected zero")
    # Ternary tree permutations and compatible cross-base conjugacy.
    for level in range(1,9):
        q = 3**level
        maps = {}
        for base in (4,10):
            values = [repunit_mod(base,j,level) for j in range(q)]
            check(len(set(values)) == q)
            for j,h in enumerate(values):
                check(inverse_repunit_residue(base,h,level) == j)
                check(h % (q//3) == repunit_mod(base,j % (q//3),level-1))
                check(values[(j+1)%q] == (base*h+1)%q)
            maps[base] = values
        psi = {maps[10][j]: maps[4][j] for j in range(q)}
        check(all(psi[(10*x+1)%q] == (4*psi[x]+1)%q for x in range(q)))
        print("Ternary repunit tree level",level,":",q,"vertices; both permutations, inverses, parent maps, conjugacy PASS")
    for base in (4,10):
        for j in range(41):
            exact = (base**j-1)//(base-1)
            expansion = sum(comb(j,h+1)*(base-1)**h for h in range(j))
            check(exact == expansion)
            check((exact-j) % (base-1) == 0)
    print("Exact binomial carry expansion for bases4,10 and lengths0..40 PASS")
    # Integral Berggren ray clock, checked in actual primitive triples.
    rays = 0
    for target in range(1,80,2):
        if target % 3 == 0:
            continue
        k0 = 2 if target % 3 == 1 else 3
        x0 = (2**k0*target-1)//3
        while x0 <= target:
            k0 += 2
            x0 = (2**k0*target-1)//3
        for j in range(7):
            h = 2**(k0-1)*(4**j-1)//3
            x = (2**(k0+2*j)*target-1)//3
            check(x == x0+2*h*target)
            check(phi(x,target)[0] == x*target)
            check(spinor_step(x,1)[0] == (target,1))
            rays += 1
    print("Collatz inverse fibre / literal Berggren ray sampling checks:",rays)
    print("Height hole: R4(j)=0,1,5,21,... omits2 despite bijection at every ternary level")
    print("No universal integer convergence or full Berggren-tree conjugacy inferred")


if __name__ == "__main__":
    main()
