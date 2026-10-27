#!/usr/bin/env python3
"""Standalone exact parity and reduced-observability tests; standard library, no external inputs."""
from collections import defaultdict
from itertools import combinations, permutations
from math import factorial
import json
from fractions import Fraction
from math import comb

def require(condition, message):
    if not condition:
        raise ValueError(message)


def tournament(n, code):
    require(n >= 1 and 0 <= code < (1 << (n*(n-1)//2)), 'bad tournament code')
    out = [0]*n
    for k, (i,j) in enumerate(combinations(range(n), 2)):
        if (code >> k) & 1:
            out[i] |= 1 << j
        else:
            out[j] |= 1 << i
    return out


def path_cover_profile(out):
    """Subset Hamiltonian DP, then root-block set-partition recurrence."""
    n = len(out)
    require(all(not (out[v] >> v & 1) for v in range(n)), 'loop')
    require(all(((out[v] >> w)&1) + ((out[w] >> v)&1) == 1
                for v in range(n) for w in range(v)), 'not a tournament')
    dp = [[0]*n for _ in range(1 << n)]
    h = [0]*(1 << n)
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1, 1 << n):
        h[mask] = sum(dp[mask])
        for v in range(n):
            rem = out[v] & ~mask
            while rem:
                bit = rem & -rem
                rem -= bit
                dp[mask | bit][bit.bit_length()-1] += dp[mask][v]
    pc = [[0]*(n+1) for _ in range(1 << n)]
    pc[0][0] = 1
    for mask in range(1, 1 << n):
        first = mask & -mask
        sub = mask
        while sub:
            if sub & first:
                for d in range(1, n+1):
                    pc[mask][d] += h[sub]*pc[mask ^ sub][d-1]
            sub = (sub-1) & mask
    return pc[-1][1:]








def phi(out):
    """T -> a -> b -> T, retaining every original named vertex."""
    n = len(out)
    require(n > 0, 'empty tournament excluded')
    return [bits | (1 << n) for bits in out] + [1 << (n+1), (1 << n)-1]


def hamiltonian_dp(out):
    n = len(out)
    dp = [[0]*n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1, 1 << n):
        for v, count in enumerate(dp[mask]):
            if not count:
                continue
            rem = out[v] & ~mask
            while rem:
                bit = rem & -rem
                rem -= bit
                dp[mask | bit][bit.bit_length()-1] += count
    return sum(dp[-1])


def matrix_entries(d):
    require(d >= 1, 'positive path index required')
    return ((d-2, 1), (d-1, 2*d-1), (d, d*(d+2)), (d+1, d*(d+1)))


def transform(p):
    require(len(p) >= 1, 'empty profile excluded')
    n = len(p)
    return [sum(coef*p[r-1] for r,coef in matrix_entries(d) if 1 <= r <= n)
            for d in range(1,n+3)]








def must_reject(call, message):
    try:
        call()
    except ValueError:
        return
    raise ValueError('hostile input accepted: '+message)


def add(a,b):
    result = [0]*max(len(a),len(b))
    for i,v in enumerate(a): result[i] += v
    for i,v in enumerate(b): result[i] += v
    while len(result)>1 and result[-1]==0: result.pop()
    return result


def mul(a,b):
    c=[0]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b): c[i+j]+=x*y
    return c


def psi(p):
    f=[0]; falling=[1]
    for d,value in enumerate(p,1):
        falling=mul(falling,[-(d-1),1])
        f=add(f,[value*x for x in falling])
    return f


def shift(f,a):
    g=[0]*len(f)
    for d,x in enumerate(f):
        for j in range(d+1): g[j]+=x*comb(d,j)*a**(d-j)
    return g


def operator(f):
    return add([0,0]+f,[0]+add(shift(f,1),[-x for x in shift(f,-1)]))


def evaluate(f,q):
    return sum(v*q**i for i,v in enumerate(f))


def solve(a,b):
    n=len(b)
    m=[[Fraction(x) for x in row]+[Fraction(y)] for row,y in zip(a,b)]
    for j in range(n):
        pivot=next((i for i in range(j,n) if m[i][j]),None)
        require(pivot is not None,'singular interpolation')
        m[j],m[pivot]=m[pivot],m[j]
        div=m[j][j]; m[j]=[x/div for x in m[j]]
        for i in range(n):
            if i!=j:
                scale=m[i][j]; m[i]=[x-scale*y for x,y in zip(m[i],m[j])]
    return [row[-1] for row in m]


def recover_polynomial(n,observations):
    require(n>=1,'nonempty order required')
    r=(n-1)//2
    require(len(observations)==r,'wrong number of observations')
    weights={1:1}; values={0:0}
    for k,h in enumerate(observations):
        require(weights.get(k+1)==factorial(k),'evaluation pivot')
        remainder=h-sum(w*values[q] for q,w in weights.items() if q<=k)
        require(remainder%factorial(k)==0,'nonintegral recovered evaluation')
        values[k+1]=remainder//factorial(k)
        nxt=defaultdict(int)
        for q,w in weights.items():
            nxt[q]+=w*q*q
            nxt[q+1]+=w*q
            if q>0: nxt[q-1]-=w*q
        weights=dict(nxt)
    powers=list(range(1 if n%2 else 2,n,2))
    require(len(powers)==r,'parity dimension')
    coefficients=solve([[q**j for j in powers] for q in range(1,r+1)],
                       [values[q]-q**n for q in range(1,r+1)])
    result=[0]*(n+1); result[n]=1
    for j,c in zip(powers,coefficients):
        require(c.denominator==1,'nonintegral polynomial coefficient')
        result[j]=int(c)
    return result


def from_psi(f):
    n=len(f)-1; p=[0]*n; residual=f[:]
    falling=[[1]]
    for d in range(1,n+1): falling.append(mul(falling[-1],[-(d-1),1]))
    for d in range(n,0,-1):
        p[d-1]=residual[d]
        for j,v in enumerate(falling[d]): residual[j]-=p[d-1]*v
    require(not any(residual),'falling factorial residual')
    return p


def long_reversal(s,n=None):
    require(s>=3,'long reversal needs at least three vertices')
    if n is None:n=s
    require(n>=s,'join order')
    out=tournament(n,(1 << (n*(n-1)//2))-1)
    out[0] &= ~(1 << (s-1));out[s-1] |= 1
    return out


def main():
    checked=0
    for n in range(1,6):
        codes=range(1 << (n*(n-1)//2)) if n<=4 else (69,76)
        for code in codes:
            out=tournament(n,code);p=path_cover_profile(out);f=psi(p)
            require(len(f)==n+1 and f[-1]==1 and f[0]==0,'monic/constant')
            require(all(v==0 for j,v in enumerate(f) if j%2!=n%2),'parity')
            require(psi(transform(p))==operator(f),'operator identity')
            # Independent permutation/descent evaluation, with ordinary binomial for q>=1.
            for q in range(1,5):
                count=0
                for word in permutations(range(n)):
                    desc=sum(not(out[word[j]] >> word[j+1]&1) for j in range(n-1))
                    top=q+n-1-desc
                    count+=comb(top,n) if top>=n else 0
                require(count==evaluate(f,q),'descent formula')
            obs=[];cur=out
            for k in range((n-1)//2):
                obs.append(hamiltonian_dp(cur));cur=phi(cur)
            recovered=recover_polynomial(n,obs)
            require(recovered==f and from_psi(recovered)==p,'reduced-probe reconstruction')
            checked+=1
    # Basis vectors tested by actual profiles, not the asserted polynomial formula.
    basis_cases=0
    for n in range(3,9):
        rows=[];r=(n-1)//2;powers=list(range(1 if n%2 else 2,n,2))
        for s in range(3,n+1,2):
            actual=psi(path_cover_profile(long_reversal(s,n)))
            inner=add(shift([0]*(s-2)+[1],1),[-x for x in shift([0]*(s-2)+[1],-1)])
            expected=add([0]*n+[1],[0]*(n-s+1)+inner)
            require(actual==expected,'long reversal/join formula')
            rows.append([actual[j] for j in powers]);basis_cases+=1
        require(len(rows)==r,'basis size')
        solve(rows,[0]*r)  # Throws if basis is singular.
    # Even s also checked directly.
    for s in (4,6,8):
        f=psi(path_cover_profile(long_reversal(s)))
        expected=add([0]*s+[1],[0]+add(shift([0]*(s-2)+[1],1),[-x for x in shift([0]*(s-2)+[1],-1)]))
        require(f==expected,'even long reversal formula');basis_cases+=1
    # Same-profile/nonisomorphic hostile pair: source-C3 and its converse.
    source_triangle=[14,4,8,2]
    opposite=[sum(1 << j for j in range(4) if source_triangle[j] >> i & 1)
              for i in range(4)]
    require(path_cover_profile(source_triangle)==path_cover_profile(opposite)==[3,9,6,1],
            'converse profile invariant')
    require(psi([3,9,6,1])==[0,0,2,0,1],'source triangle polynomial')
    require(sorted(x.bit_count() for x in source_triangle)==[1,1,1,3] and
            sorted(x.bit_count() for x in opposite)==[0,2,2,2],
            'nonisomorphism by degree sequences')
    # Operator identity on each basis vector, not just realizable profiles.
    for a in range(1,11):
        basis=[0]*(a-1)+[1]
        require(psi(transform(basis))==operator(psi(basis)), 'falling basis operator')
    # The all-integer ideal argument is proved in the supplement; iterations illustrate it.
    f_a=psi([15,39,33,10,1]);f_b=psi([15,45,35,10,1])
    delta=add(f_b,[-x for x in f_a])
    require(delta==[0,-2,0,2], 'symbolic modular obstruction')
    require(all(evaluate(delta,q)%12==0 for q in range(-20,21)), 'base integer divisibility')
    for k in range(12):
        require(evaluate(delta,1)%12==0, 'modular observation obstruction')
        delta=operator(delta)
    must_reject(lambda:recover_polynomial(5,[15]),'insufficient probes')
    must_reject(lambda:recover_polynomial(0,[]),'empty order')
    must_reject(lambda:recover_polynomial(5,[15,124]),'fractional interpolation')
    print(json.dumps({'status':'PASS','tournaments':checked,'basis_cases':basis_cases,
                      'hostile_rejections':3,'operator_basis_cases':10,'modular_illustration_iterations':12,'converse_boundary_check':True,'five_vertex_probe_vectors':[[15,123],[15,135]]},sort_keys=True,indent=2))


if __name__=='__main__': main()
