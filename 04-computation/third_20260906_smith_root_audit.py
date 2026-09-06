#!/usr/bin/env python3
"""Full integer Hasse matrices, independent SymPy Smith forms; no producer import."""
from fractions import Fraction
from itertools import combinations
from math import comb, gcd
from sympy import Matrix, ZZ
from sympy.matrices.normalforms import smith_normal_form

GATES = 0


def need(ok, why):
    global GATES
    GATES += 1
    if not ok:
        raise RuntimeError(why)


def val(x, p):
    x = abs(int(x))
    if x == 0:
        return 10**9
    r = 0
    while x % p == 0:
        x //= p
        r += 1
    return r


def hasse_matrix(m, frames):
    degree = m+2
    rows = []
    for ((vx,vy),(wx,wy)), multiplicity in zip(frames,(m,2,1)):
        need(vx*wy-vy*wx == 1, 'unimodular tangent')
        for r in range(multiplicity):
            row = []
            for k in range(degree+1):
                row.append(sum(comb(degree-k,i)*comb(k,r-i)
                    *vx**(degree-k-i)*wx**i*vy**(k-r+i)*wy**(r-i)
                    for i in range(r+1) if i <= degree-k and 0 <= r-i <= k))
            rows.append(row)
    return Matrix(rows)


def affine_frames(a,b):
    return [((1,x),(0,1)) for x in (0,a,b)]


def smith(M):
    S = smith_normal_form(M,domain=ZZ)
    return tuple(abs(int(S[i,i])) for i in range(S.rows))


def global_formula(m,a,b):
    first = gcd(a**m,m*a**(m-1),b**m)
    second = gcd(a**(2*m),a**m*b**m*(b-a),
                 a**(m-1)*b**m*(m*b-(m+1)*a))
    determinant = abs(a)**(2*m)*abs(b)**m*(b-a)**2
    return (1,)*m+(first,second//first,determinant//second)


def classified(m,a,b,p):
    A,B,C,mu = val(a,p),val(b,p),val(b-a,p),val(m,p)
    if A <= B:
        h = min(A,mu)
        rest = ((m-1)*A+h,(m+1)*A-h,m*B+2*C)
    else:
        e,d = B,A-B
        if d != mu:
            rho = (m-1)*e+min(e,(m-1)*d+mu)
            sigma = 2*m*e+(m-1)*d+min(d,mu)
            total = 2*m*A+m*B+2*C
            rest = (rho,sigma-rho,total-sigma)
        else:
            K = min(e,m*d)
            u,v = a//p**(e+d),b//p**e
            kappa = min(K,val((m//p**d)*v-(m+1)*u,p))
            rest = ((m-1)*e+K,(m+1)*e+m*d-K+kappa,
                    (m+2)*e+m*d-kappa)
    return (0,)*m+rest


def main():
    matrices = partitions = projective = inverses = 0
    bank = {(m,a,b) for m in range(2,13) for a in range(-6,7)
            for b in range(-6,7) if a and b and a != b}
    for m,p in ((3,3),(4,2),(5,5),(6,2),(8,2),(9,3),(10,5),(12,3),(16,2)):
        d = val(m,p)
        for e in sorted({0,1,2,3,m*d,m*d+1}):
            K = min(e,m*d)
            for k in sorted({0,1,K,max(0,K-1)}):
                c = m//p**d
                candidates = [c+p**k,c+2*p**k,c]
                for u in candidates:
                    if u % p:
                        a,b = p**(e+d)*u,p**e*(m+1)
                        bank.add((m,a,b))
    for m,a,b in sorted(bank):
        M = hasse_matrix(m,affine_frames(a,b))
        result = smith(M)
        need(result == global_formula(m,a,b),'entire integer Smith list')
        need(all(y % x == 0 for x,y in zip(result,result[1:])),'ordered divisibility')
        for p in (2,3,5,7):
            need(tuple(val(x,p) for x in result) == classified(m,a,b,p),'full local case split')
            partitions += 1
        matrices += 1

    transforms = ((1,1,0,1),(1,0,1,1),(0,-1,1,0),(2,1,1,1))
    for m in (2,3,4,7,11):
        for a,b in ((8,4),(9,12),(16,-8),(81,108),(3,2)):
            expected = global_formula(m,a,b)
            for aa,bb,cc,dd in transforms:
                frames = []
                for i,((vx,vy),(wx,wy)) in enumerate(affine_frames(a,b)):
                    v = (aa*vx+bb*vy,cc*vx+dd*vy)
                    w = (aa*wx+bb*wy,cc*wx+dd*wy)
                    w = (w[0]+i*v[0],w[1]+i*v[1])
                    frames.append((v,w))
                need(smith(hasse_matrix(m,frames)) == expected,'full projective/tangent covariance')
                projective += 1
        frames = [((1,0),(0,1)),((0,1),(-1,0)),((1,1),(0,1))]
        need(smith(hasse_matrix(m,frames)) == (1,)*(m+3),'three dyadic residue classes')
        projective += 1

    for m,a,b in ((2,8,4),(3,81,108),(3,162,108),(4,32,12),(5,125,30),
                  (6,144,56),(8,256,36),(9,243,90)):
        M = hasse_matrix(m,affine_frames(a,b))
        inverse = M.inv()
        denominators = [1]*(m+3)
        for j in range(m+3):
            for i in range(m+3):
                q = int(inverse[i,j].q)
                denominators[j] = denominators[j]*q//gcd(denominators[j],q)
        for p in (2,3,5,7):
            loss = val(global_formula(m,a,b)[-1],p)
            need(max(val(denominators[j],p) for j in (m,m+2)) == loss,
                 'small-bank value owns full inverse loss')
            if loss:
                j = next(j for j in (m,m+2) if val(denominators[j],p) == loss)
                source = inverse[:,j]*denominators[j]
                need(all(x.q == 1 for x in source),'integral attained value perturbation')
                need(min(val(x,p) for x in source) == 0,'attained source precision')
                image = M*source
                need(all(image[i] == (denominators[j] if i == j else 0)
                         for i in range(m+3)),'only selected value changes')
        inverses += 1

    windows = 0
    for m,p in ((2,2),(3,3),(4,2),(5,5),(8,2),(9,3)):
        d = val(m,p)
        for e in range(1,m*d+3):
            K = min(e,m*d)
            for k in range(int(p==2),K):
                P,Q = (m+1)*e+m*d-K,(m+2)*e+m*d
                for N in range(1,Q+2):
                    diff = min(N,P+k+1)+min(N,Q-k-1)-min(N,P+k)-min(N,Q-k)
                    need(diff == int(P+k+1 <= N <= Q-k-1),'full adjacent kernel window')
                    windows += 1
    print('THIRD_SESSION_INDEPENDENT_MIXED_SMITH_AUDIT')
    print(f'full_integer_matrices={matrices} full_prime_partitions={partitions}')
    print(f'full_projective_controls={projective} full_inverse_owner_controls={inverses}')
    print(f'adjacent_kernel_precision_checks={windows}')
    print(f'PASS gates={GATES}')


if __name__ == '__main__':
    main()
