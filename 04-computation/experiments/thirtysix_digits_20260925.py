"""Literal decimal phases and the complete signed odd inverse fibre.

Standard library only; all checks remain active under python -O.
"""

from fractions import Fraction
from math import gcd, isqrt, prod
from pathlib import Path


def check(ok, message="check failed"):
    if not ok:
        raise RuntimeError(message)


def v2(n):
    check(n > 0)
    return (n & -n).bit_length()-1


def prime_trial(n):
    if n < 2:
        return False
    if n % 2 == 0:
        return n == 2
    return all(n % d for d in range(3, isqrt(n)+1, 2))


def factor(n):
    out = {}
    d = 2
    while d*d <= n:
        while n % d == 0:
            out[d] = out.get(d, 0)+1
            n //= d
        d += 1
    if n > 1:
        out[n] = out.get(n, 0)+1
    return out


def order(a, p):
    check(prime_trial(p) and a % p)
    d = p-1
    for q in factor(d):
        while d % q == 0 and pow(a, d//q, p) == 1:
            d //= q
    check(pow(a, d, p) == 1)
    check(all(pow(a, d//q, p) != 1 for q in factor(d)))
    return d


def a(k):
    check(k >= 1)
    return (10**k-7)//3


def U(n, sigma):
    check(n > 0 and n % 2 and sigma in (-1, 1))
    z = 3*n+sigma
    k = v2(z)
    return z >> k, k


def F(n, sigma):
    return 2*n+sigma, -sigma


def source(y, k):
    check(y > 0 and y % 2 and y % 3 and k >= 1)
    r = pow(2, k, 3)*(y % 3) % 3
    sigma = 1 if r == 1 else -1
    n = (2**k*y-sigma)//3
    check(n > 0 and n % 2 and U(n, sigma) == (y, k))
    return n, sigma


def main():
    lines = []

    def report(s):
        lines.append(str(s))

    report("THIRTYSIX DIGIT AND SIGNED-FIBRE AUDIT 2026-09-25")
    report("Universe: decimal k1..120, exact prime return certificate k18; full signed source germs odd n<=10001; targets odd y<=501,3 not dividing y, exponents1..60; residue cycles depths1..7; literal parabolic heights y<=99,k1..8.")
    check(all(prime_trial(a(k)) for k in range(2, 9)))
    witnesses = dict(zip(range(9, 18), (17,673,307,19,523,607,181,199,31)))
    for k,p in witnesses.items():
        check(1 < p < a(k) and a(k) % p == 0)
    report(f"Seven prime terms k2..8; composite witnesses k9..17={list(witnesses.items())}")
    n = a(18)
    qs = (2,3,5,2071723,5363222357)
    check(prod(qs) == n-1 and all(prime_trial(q) for q in qs))
    check(pow(2,n-1,n) == 1)
    cert = [(q,pow(2,(n-1)//q,n),gcd(pow(2,(n-1)//q,n)-1,n)) for q in qs]
    check(all(g == 1 for _,_,g in cert))
    report(f"PROVED prime return a18={n}; n-1 factors={qs}; base2 Lucas order certificate(q,residue,gcd)={cert}")
    report("Completeness proof: for any prime divisor ell of n, these tests force ord_ell(2)=n-1; hence ell>=n, so n is prime.")
    for k in range(1, 121):
        check(a(k+1) == 10*a(k)+21)
        check(a(k+1)-a(k) == 3*10**k)
        check(Fraction(a(k+1),a(k)) == 10+Fraction(21,a(k)))
        check(a(k+8) % 17 == (1-a(k)) % 17)
        check((a(k)%17 == 0) == (k%16 == 9))
        check(all(a(k) % p for p in (2,3,5,7,11,13)))
    residues = [a(k)%17 for k in range(1,17)]
    check(set(residues) == set(range(17))-{9})
    check(order(10,17)==16 and pow(10,8,17)==16)
    check(order(2,17)==8 and pow(2,4,17)==16 and order(4,17)==4)
    report(f"a_k mod17,k1..16={residues}; excluded fixed class9, eight-step antipode a->1-a, sixteen-step return.")
    report("Literal gaps=3*10^k and ratios=10+21/a_k. Carry21 becomes42/84 only after multiplying the entire recurrence by2/4.")
    orders = [(k,a(k),order(10,a(k))) for k in range(2,9)]
    for k,p,d in orders:
        check(pow(10,k,p)==7 and pow(10,k+d,p)==7)
    report(f"Each initial prime seeds a later composite phase: (k,p,ord_p(10))={orders}; e.g.31 divides a_k exactly at k=2 mod15, with k2 the prime itself and k17 composite.")

    source_checks = 0
    for n in range(1,10002,2):
        for sigma in (-1,1):
            y,k = U(n,sigma)
            check(y%3 and source(y,k)==(n,sigma))
            nn,ss = F(n,sigma)
            check(U(nn,ss)==(y,k+1))
            inverse_legal = (n+sigma)%4 == 2
            check(inverse_legal == (k>=2))
            if inverse_legal:
                parent=((n+sigma)//2,-sigma)
                check(F(*parent)==(n,sigma) and U(*parent)==(y,k-1))
            source_checks += 1
    fibre_checks = 0
    for y in range(1,502,2):
        if y%3==0:
            continue
        for k in range(1,61):
            n,sigma=source(y,k)
            check(F(n,sigma)==source(y,k+1))
            for r in (0,1,2,3,6):
                state=(n,sigma)
                for _ in range(r): state=F(*state)
                expected=(2**r*n+(2**r-(-1)**r)*sigma//3,(-1)**r*sigma)
                check(state==expected and U(*state)==(y,k+r))
            fibre_checks += 1
    root=[source(1,k) for k in range(1,13)]
    check(root[:4]==[(1,-1),(1,1),(3,-1),(5,1)])
    check(all(n==(2**k-(-1)**k)//3 for k,(n,sigma) in enumerate(root,1)))
    report(f"All signed-germ decoder checks={source_checks}; inverse-ray exponent checks={fibre_checks}; root fibre k1..12={root}")
    report("F=(2n+sigma,-sigma) preserves odd TARGET VALUE and adds1 to exact halving count; F^2=(4n+sigma,sigma), F^3=(8n+3sigma,-sigma), F^6=(64n+21sigma,sigma).")
    for depth in range(1,8):
        modulus=2*3**depth
        state=(1,1)
        seen=set()
        for _ in range(2*3**depth):
            check(state not in seen)
            seen.add(state)
            x,sig=F(*state)
            state=(x%modulus,sig)
        check(state==(1,1) and len(seen)==2*3**depth)
        check(seen=={(n,sig) for n in range(1,modulus,2) for sig in (-1,1)})
    report("Joint odd residue/sheet clock: one complete cycle of length2*3^d modulo2*3^d, depths1..7.")
    parabolic=half_obstructions=0
    for y in range(1,100,2):
        if y%3==0:
            continue
        for k in range(1,9):
            n,sigma=source(y,k)
            if n<=y:
                continue
            nn,ss=F(*F(n,sigma))
            check(ss==sigma and nn==n+2*(2**(k-1))*y)
            check(gcd(n,y)==gcd(nn,y)==1)
            halfheight=Fraction(F(n,sigma)[0]-n,2*y)
            check((halfheight.denominator==1)==(y==1))
            parabolic+=1
            half_obstructions+=(y>1)
    report(f"Same-sheet literal Berggren parabolic increments=2^(k-1) checked at{parabolic} points; half-step has nonintegral same-ray height for all{half_obstructions} checked y>1 points.")
    check(U(5,-1)[0]==U(9,1)[0]==7)
    check(U(7,-1)[0]==5 and U(7,1)[0]==11)
    report("Forward-conjugacy hostile:5 on minus andF(5,-)=(9,+) share output7, but next fixed-sheet outputs are5 and11. The fibre decoder is not a time-map conjugacy.")

    for k in range(5,121):
        n=a(k); x,e1=U(n,1); y,e2=U(x,1)
        check((e1,e2)==(1,3) and y==(3*10**k)//16-1 and y<n)
        check(v2(y+1)==k-4)
        z=y
        for j in range(1,k-4):
            z,e=U(z,1)
            check(e==1 and z==3**(j+1)*2**(k-4-j)*5**k-1 and z>y)
        if k>=7:
            z=U(U(y,1)[0],1)[0]
            check(z>n)
    report("Decimal carry refinement: U+^2(a_k)+1=3*2^(k-4)*5^k for k>=5; the next k-5 odd steps all grow. For every k>=7, U+^4(a_k)>a_k again.")
    report("ALL CHECKS PASSED. Infinite statements have proofs in the note; decimal phase, fibre height, and forward orbit time are separate clocks.")
    output="\n".join(lines)+"\n"
    path=Path(__file__).resolve().parents[2]/"05-knowledge/results/thirtysix_digits_20260925.out"
    path.write_text(output,encoding="utf-8")
    print(output,end="")


if __name__=="__main__":
    main()
