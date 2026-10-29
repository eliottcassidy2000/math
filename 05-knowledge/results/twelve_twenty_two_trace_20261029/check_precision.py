"""Exact 13-adic clock tower, rank-four mixed probe, and historical repair."""
from itertools import combinations
from math import prod
import json
from sympy import Matrix, divisors, isprime, kronecker_product

from check_clocks import HERE, ONE, U, mul, power, tr, need, rank_mod, smith


def layers(rows, p, precision):
    modulus = p**precision
    A = [[int(x) % modulus for x in row] for row in rows]
    nr, nc, rank = len(A), len(A[0]), 0
    while rank < min(nr,nc):
        pivot = next(((i,j) for i in range(rank,nr) for j in range(rank,nc) if A[i][j] % p),None)
        if pivot is None:
            break
        i,j = pivot
        A[rank],A[i] = A[i],A[rank]
        for row in A:
            row[rank],row[j] = row[j],row[rank]
        inverse = pow(A[rank][rank],-1,modulus)
        A[rank] = [inverse*x % modulus for x in A[rank]]
        for i in range(rank+1,nr):
            c = A[i][rank]
            A[i] = [(x-c*y) % modulus for x,y in zip(A[i],A[rank])]
        rank += 1
    if precision == 1:
        return [rank]
    if rank == min(nr,nc):
        return [rank]+[0]*(precision-1)
    tail = [[A[i][j] for j in range(rank,nc)] for i in range(rank,nr)]
    need(all(x % p == 0 for row in tail for x in row),'primary residual divisible')
    return [rank]+layers([[x//p for x in row] for row in tail],p,precision-1)


def run():
    u12 = power(U,12)
    need(all((a-b) % 13 == 0 for a,b in zip(u12,ONE)), 'scalar clock at split prime13')
    A = tuple((a-b)//13 for a,b in zip(u12,ONE))
    change = [ONE,A,mul(A,A)]
    determinant = int(Matrix(change).det())
    need(abs(determinant) == 29*43 and determinant % 13 != 0,
         'normalized collision has an integral local power basis')
    tower = []
    for r in range(9):
        s, exponent = r+1,12*13**r
        precision, modulus = 2*s+1,13**(2*s+1)
        O = [[tr(power(U,i+j*exponent,modulus)) % modulus for i in range(3)] for j in range(3)]
        actual = layers(O,13,precision)
        expected = [int(j in (0,s,2*s)) for j in range(precision)]
        need(actual == expected,'all certified primary precision layers')
        joined = [[tr(power(U,i+t,modulus)) % modulus for i in range(3)]
                  for t in (0,exponent,2*exponent,22,44)]
        need(rank_mod(joined,13) == 3,'fixed second clock restores the complete local observer')
        tower.append(dict(r=r,exponent=exponent,smith_valuations=[0,s,2*s],layers=actual))

    # Ramanujan130 factors: a=2+sqrt5, b=(3+sqrt13)/2.
    # Basis order is (1,b,a,ab).
    Ta = Matrix([[2,4],[4,18]])
    Tb = Matrix([[2,3],[3,11]])
    G = kronecker_product(Ta,Tb)
    z = Matrix([12,-8,-6,4])  # (2a-4)(2b-3)=2 sqrt65.
    values = list(map(int,G*z))
    need(values == [0,0,0,260],'all separate-clock moments vanish but a mixed one is visible')
    need(smith(G.tolist()) == [2,2,130,130], 'rank-four trace lattice')
    need(rank_mod(G.tolist()[:3],3)==3 and rank_mod(G.tolist(),3)==4,
         'one missing ternary amplitude at the first rank-four boundary')
    Ma = Matrix([[0,1],[1,4]])
    Mb = Matrix([[0,1],[1,3]])
    aop,bop = kronecker_product(Ma,Matrix.eye(2)),kronecker_product(Matrix.eye(2),Mb)
    trace_row = Matrix([G.tolist()[0]])
    for j in range(25):
        need((trace_row*aop**j*z)[0] == (trace_row*bop**j*z)[0] == 0,
             'separate-power controls; quadratic relations prove all degrees')
    need((trace_row*aop*bop*z)[0] == 260,'literal mixed clock response')

    def selector(n):
        return [d+1 for d in divisors(n) if isprime(d+1)]
    bernoulli = [dict(n=n,primes=selector(n),next=prod(selector(n))) for n in (6,8,12,22,42,1806)]
    need([r['next'] for r in bernoulli] == [42,30,2730,138,1806,1806], 'corrected divisor-prime map')
    need(not isprime(9) and prod(selector(8)) != 8,'composite successor does not imply a fixed point')
    need([p for p in selector(12) if p%7 in (1,6)] == [13], 'Bernoulli12 intersection with split cubic primes')
    return dict(normalized_clock=A,normalized_basis_determinant=determinant,
                precision_tower=tower,mixed_gram=G.tolist(),mixed_witness=list(z),mixed_responses=values,
                bernoulli=bernoulli)


if __name__ == '__main__':
    import sys
    result=json.loads(json.dumps(run(),default=int))
    path=HERE/'precision.json'
    if '--write' in sys.argv:
        path.write_text(json.dumps(result,separators=(',',':'))+'\n',encoding='utf-8')
    else:
        need(result==json.loads(path.read_text(encoding='utf-8')),'unchanged precision and archive-repair certificate')
    print('PASS 13-adic precision recurrence; mixed130 probe; March Bernoulli correction')
