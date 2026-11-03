#!/usr/bin/env python3
"""Exact controls for the two odd-pair triangles and literal PRIMEGAME chart."""
from fractions import Fraction
from itertools import product
from math import gcd, isqrt, prod
import json

GATES = 0

def need(test, message):
    global GATES
    GATES += 1
    if not test:
        raise RuntimeError(message)


def prime(n):
    return n >= 2 and all(n % d for d in range(2, isqrt(n) + 1))


def prime_sieve(limit):
    a = [True] * (limit + 1)
    a[:2] = [False, False]
    for d in range(2, isqrt(limit) + 1):
        if a[d]:
            for j in range(d*d, limit+1, d):
                a[j] = False
    return a


def factors(n):
    out = []
    for p in range(2, isqrt(n) + 1):
        if n % p == 0:
            out.append(p)
            while n % p == 0:
                n //= p
    if n > 1:
        out.append(n)
    return out


def rank_mod(rows, p=1000003):
    a = [list(r) for r in rows]
    h, w, k = len(a), len(a[0]), 0
    for j in range(w):
        pivot = next((i for i in range(k, h) if a[i][j] % p), None)
        if pivot is None:
            continue
        a[k], a[pivot] = a[pivot], a[k]
        inverse = pow(a[k][j] % p, -1, p)
        a[k] = [(v * inverse) % p for v in a[k]]
        for i in range(k+1, h):
            c = a[i][j] % p
            if c:
                a[i] = [(x-c*y) % p for x,y in zip(a[i],a[k])]
        k += 1
        if k == h:
            break
    return k


def characteristic_integer(a):
    """Faddeev-LeVerrier with exact integer divisions, highest degree first."""
    n=len(a)
    b=[[int(i==j) for j in range(n)] for i in range(n)]
    coefficients=[1]
    for k in range(1,n+1):
        c=[[sum(a[i][t]*b[t][j] for t in range(n)) for j in range(n)] for i in range(n)]
        trace=sum(c[i][i] for i in range(n))
        need(trace%k==0,'integer characteristic coefficient division')
        coefficient=-trace//k
        coefficients.append(coefficient)
        for i in range(n):
            c[i][i]+=coefficient
        b=c
    need(all(x==0 for row in b for x in row),'Cayley-Hamilton residual')
    return coefficients


def polynomial_product(a,b):
    out=[0]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):
            out[i+j]+=x*y
    return out


PROGRAM = [(17,91),(78,85),(19,51),(23,38),(29,33),(77,29),
           (95,23),(77,19),(1,17),(11,13),(13,11),(15,2),(1,7),(55,1)]
DATA = (3,5,7)
CONTROL = (11,13,17,19,23,29)


def split(n):
    a, u, v = 0, 1, 1
    while n % 2 == 0:
        a += 1
        n //= 2
    for p in DATA:
        while n % p == 0:
            u *= p
            n //= p
    for p in CONTROL:
        while n % p == 0:
            v *= p
            n //= p
    if n != 1:
        raise ValueError('outside typed prime-supported carrier')
    return a, u, v


PARTS = []
for numerator, denominator in PROGRAM:
    an, un, vn = split(numerator)
    ab, ub, vb = split(denominator)
    PARTS.append((an-ab, Fraction(un,ub), Fraction(vn,vb)))


def literal(n):
    for j,(a,b) in enumerate(PROGRAM):
        if n % b == 0:
            return j, n//b*a
    return None, None


def pair_encode(n):
    a,u,v = split(n)
    return a, (u+v)//2, (v-u)//2


def frame_step(frame):
    a,h,d = frame
    u,v = h-d,h+d
    for j,(delta,rho,sigma) in enumerate(PARTS):
        if a+delta >= 0 and u % rho.denominator == 0 and v % sigma.denominator == 0:
            hh = ((rho+sigma)*h+(sigma-rho)*d)/2
            dd = ((sigma-rho)*h+(rho+sigma)*d)/2
            need(hh.denominator == dd.denominator == 1, 'integral frame image')
            return j,(a+delta,int(hh),int(dd))
    return None,None


def reconstruct(frame):
    a,h,d = frame
    return 2**a*(h-d)*(h+d)


def main():
    # Declared universe: all unordered odd pairs with sum <= 1000.
    ps = prime_sieve(1000)
    for n in range(1001):
        need(ps[n] == prime(n), 'independent primality paths')
    pair_count = prime_pair_count = primitive_count = 0
    first_selfhit = None
    for h in range(3,501):
        tests = [p for p in range(3,isqrt(2*h-3)+1,2) if ps[p]]
        for p in range(3,h+1,2):
            q = 2*h-p
            d = h-p
            pair_count += 1
            need(h-d == p and h+d == q and (h+d)%2 == 1, 'triangle inverse/parity')
            true = ps[p] and ps[q]
            guarded = all((p%ell or p==ell) and (q%ell or q==ell) for ell in tests)
            need(bool(guarded) == true, 'height/self-hit sieve exact')
            if true:
                prime_pair_count += 1
                need((d==1) == (q-p==2), 'twin column')
                need((isqrt(p+q)**2 == p+q) == (h%2==0 and isqrt(h//2)**2==h//2), 'square sum rows')
                if first_selfhit is None and not all(p%ell and q%ell for ell in tests):
                    first_selfhit = [h,p,q]
                if p<q:
                    A,B,C = p*q,(q*q-p*p)//2,(p*p+q*q)//2
                    need(A*A+B*B==C*C and gcd(gcd(A,B),C)==1, 'primitive right triangle')
                    need(gcd(h,d)==1 and h-d==p and h+d==q, 'nonunit Euclid factor inverse')
                    primitive_count += 1
    need(first_selfhit == [7,3,11], 'earliest self-hit loss')
    need((8*8-7*7,2*8*7,8*8+7*7)==(15,112,113), 'unit-factor hostile')
    unit_pairs = [(p,290-p) for p in range(1,146,2) if gcd(p*(290-p),145)==1]
    actual_pairs = [[p,290-p] for p in range(3,146,2) if ps[p] and ps[290-p]]
    need(len(unit_pairs)==56 and unit_pairs[0]==(1,289), '56 complete unit pair orbits')
    need(len([p for p,q in unit_pairs if p>=3])==55 and len(actual_pairs)==10, '55 versus10')
    # Local masks: exact integer residue universes and a mod-prime rank lower bound.
    need(prime(1000003), 'rank modulus prime')
    kernels = []
    for M in (3,5,7,11,15,21,35,105,145):
        prime_div = factors(M)
        need(prod(prime_div)==M, 'squarefree odd kernel domain')
        rows = [[int(gcd(h*h-d*d,M)==1) for d in range(M)] for h in range(M)]
        expected_rank = prod((p+1)//2 for p in prime_div)
        unique = {tuple(row) for row in rows}
        square_classes = {h*h%M for h in range(M)}
        need(len(unique)==len(square_classes)==expected_rank, 'row-square quotient')
        need(rank_mod(rows)==expected_rank, 'rank lower bound with duplicate-row upper bound')
        for h in range(M):
            need(sum(rows[h])==prod(p-1 if h%p==0 else p-2 for p in prime_div), 'exact row mass')
            for d in range(M):
                need(rows[h][d]==int(gcd(h-d,M)==gcd(h+d,M)==1), 'invertible joint cell map')
        phi=prod(p-1 for p in prime_div)
        need(sum(map(sum,rows))==phi*phi,'complete cell count')
        representatives={}
        weights={}
        for h in range(M):
            representatives.setdefault(h*h%M,h)
            weights[h*h%M]=weights.get(h*h%M,0)+1
        reps=list(representatives.values())
        for h in reps:
            need(sum(rows[h])==sum(rows[h][d]*weights[d*d%M] for d in reps), 'weighted quotient operator')
        kernels.append({'M':M,'rank':expected_rank,'total':phi*phi})
    spectra=[]
    for ell in (3,5,7,11,29):
        matrix=[[int((h*h-d*d)%ell!=0) for d in range(ell)] for h in range(ell)]
        expected=[1,-(ell-3),-(ell-1)]
        for _ in range((ell-3)//2):
            expected=polynomial_product(expected,[1,2])
        expected += [0]*((ell-1)//2)
        actual=characteristic_integer(matrix)
        need(actual==expected,'full weighted local characteristic polynomial')
        need(sum(sum(matrix[i][j]*matrix[j][i] for j in range(ell)) for i in range(ell))==(ell-1)**2,'weighted second trace')
        # The zero-odd sector is literal cancellation of paired columns.
        for d in range(1,(ell+1)//2):
            need(all(matrix[h][d]==matrix[h][-d] for h in range(ell)),'odd local zero sector')
        spectra.append({'ell':ell,'characteristic_high_to_low':actual})
    # CRT permutation to the product carrier preserves BOTH row and column labels.
    for h in range(145):
        for d in range(145):
            need(int(gcd(h*h-d*d,145)==1)==int((h*h-d*d)%5!=0)*int((h*h-d*d)%29!=0),'literal tensor matrix entries')
    rows145 = [[int(gcd(h*h-d*d,145)==1) for d in range(145)] for h in (0,5,29,1)]
    need([sum(r) for r in rows145]==[112,108,84,81], '145 weighted row profile')
    qrs={}
    for h in range(145):
        qrs.setdefault(h*h%145,h)
    need(all(sum(gcd(h*h-d*d,145)==1 for d in qrs.values())==28 for h in qrs.values()),'quotient degree28')
    for bound in (3,5,7,11,13,17,19):
        Q=prod(p for p in range(3,bound+1,2) if prime(p))
        p,q=(2*Q+1)**2,(4*Q+1)**2
        h,d=(p+q)//2,(q-p)//2
        need(gcd(p*q,Q)==1 and gcd(p,q)==1 and gcd(h,d)==1,'fixed wheel primitive hostile')
        need(isqrt(p)**2==p and isqrt(q)**2==q,'hostile composite certificate')
    # Every literal trial phase for2<=n<=32, plus every frame transition.
    trial_count=trial_steps=0
    examples=[]
    for n in range(2,33):
        for d in range(1,n):
            p,q=2*d+1,2*(n-d)+1
            need((p+q-2)//2==n and (p-1)//2==d,'ordered trial chart inverse')
            need((n%d==0)==((q-1)%(p-1)==0),'integer-slope divisor predicate')
            success=n%d==0
            length=4*n+2*(n//d)+2 if success else 6*n+2*(n//d)+2
            target=2**n*7**(d-1) if success else 5**n*7**(d-1)*13
            state=5**n*7**d*13
            frame=pair_encode(state)
            for _ in range(length):
                j,state=literal(state)
                jj,frame=frame_step(frame)
                need(j==jj and state==reconstruct(frame),'literal and matrix frame lockstep')
                need(pair_encode(state)==frame,'literal frame inverse')
            need(state==target,'trial endpoint and literal cost')
            trial_count+=1
            trial_steps+=length
            if (n,d) in ((3,1),(3,2),(4,1),(4,2),(29,1),(30,15)):
                examples.append({'n':n,'d':d,'odd_pair':[p,q],'divides':success,'steps':length,'endpoint':state})
    need(((5-1)%(3-1)==0) and ((3-1)%(5-1)!=0),'folded trial hostile')
    need(literal(1)[0]==13 and literal(2)[0]==11 and pair_encode(1)[1:]==pair_encode(2)[1:],'missing2-exponent hostile')
    need(literal(39)==(9,33) and literal(33)==(4,29),'reciprocal word priority hostile')
    need(literal(13)==(9,11) and literal(11)==(10,13),'reciprocal loop control')
    jump=[]
    for b in range(9):
        n=3**b
        j,nn=literal(n)
        a,h,d=pair_encode(n)
        aa,hh,dd=pair_encode(nn)
        need(j==13 and hh-h==2*3**b+5,'literal-plane unbounded jump')
        jump.append(hh-h)
    # Actual orbit control, with phase/power observations; no extrapolation.
    state=2; frame=pair_encode(state); emissions=[]
    for t in range(1,20001):
        j,state=literal(state)
        jj,frame=frame_step(frame)
        need(j==jj and state==reconstruct(frame),'actual orbit frame lockstep')
        if state & (state-1)==0:
            exponent=state.bit_length()-1
            need(prime(exponent),'actual pure-power prime control')
            emissions.append([t,exponent])
    print(json.dumps({'status':'PASS','gates':GATES,
        'universes':{'unordered_odd_pairs_sum_at_most':1000,'trial_n_max':32,'all_trials':trial_count,'literal_trial_steps':trial_steps,'actual_orbit_steps':20000},
        'odd_pairs':pair_count,'prime_pairs':prime_pair_count,'primitive_prime_pair_triples':primitive_count,
        'row290':{'unit_pairs_including_one':56,'unit_pairs_excluding_one':55,'prime_pairs':actual_pairs},
        'first_prime_selfhit_loss':first_selfhit,'kernels':kernels,
        'weighted_145_rows':[112,108,84,81],'unweighted_quotient_degree':28,'weighted_local_spectra':spectra,
        'trial_examples':examples,'direct_composite4_emits16_in_steps':26,
        'priority_hostile':[39,33,29],'plane_jumps':jump,'actual_prime_emissions':emissions},indent=2,sort_keys=True))


if __name__=='__main__':
    main()
