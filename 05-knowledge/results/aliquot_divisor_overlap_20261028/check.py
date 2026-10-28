"""Exact divisor arity, factor-pair compression, and inherited 1001 controls."""
from array import array
from collections import Counter
from math import isqrt,prod
import json
import hashlib,sys

def need(b,s):
    if not b:raise RuntimeError(s)

def factor(n):
    out={};p=2
    while p*p<=n:
        while n%p==0:out[p]=out.get(p,0)+1;n//=p
        p+=1
    if n>1:out[n]=out.get(n,0)+1
    return out

def divisors(n):
    out=[]
    for d in range(1,isqrt(n)+1):
        if n%d==0:
            out.append(d)
            if d*d!=n:out.append(n//d)
    return sorted(out)

def sigma_factor(n):return prod(sum(p**j for j in range(e+1)) for p,e in factor(n).items())
def signs(n):
    es=list(factor(n).values());lam=(-1)**sum(es)
    return (0 if any(e>1 for e in es) else lam),lam

def build(N):
    sig=array('I',[0])*(N+1);tau=array('I',[0])*(N+1)
    for d in range(1,N+1):
        for n in range(d,N+1,d):sig[n]+=d;tau[n]+=1
    return sig,tau

def digest(a):
    need(a.itemsize==4,'uint32 array encoding')
    b=array('I',a)
    if sys.byteorder!='little':b.byteswap()
    return hashlib.sha256(b.tobytes()).hexdigest()

def roots(n,s):
    delta=(s-1)**2-4*n
    if delta<0:return None
    r=isqrt(delta)
    if r*r!=delta:return None
    need((s-1-r)%2==0,'Vieta parity')
    a,b=(s-1-r)//2,(s-1+r)//2
    need(0<a<=b and a*b==n and 1+a+b==s,'positive compression roots')
    return a,b

def main():
    N=1000001;sig,tau=build(N);s=[sig[n]-n for n in range(N+1)]
    powers=set()
    for p in range(2,N+1):
        if tau[p]!=2:continue
        n=p
        while n<=N:powers.add(n);n*=p
    compressed=[];unit=[];overlap=[];classes=Counter();arities=Counter()
    for n in range(2,N+1):
        r=roots(n,s[n]);c=r is not None
        need(c==(tau[n]==4 or s[n]==n+2),'all-source compression characterization')
        if c:compressed.append(n)
        if s[n]==n+2:unit.append(n);need(r==(1,n),'unit pair')
        if tau[n]==6 and c:overlap.append(n)
        label='prime_power' if n in powers else ('factor_pair_compression' if c else 'remaining')
        classes[label]+=1
        if n not in powers and not c:need(tau[n]-1>=5,'residual divisor arity')
        arities[tau[n]-1]+=1
    need(overlap==[20],'unique five-proper-divisor overlap in finite universe')
    arity_overlaps={t:[n for n in compressed if tau[n]-1==t] for t in [5,7,9,11]}
    need(arity_overlaps[7]==[104] and arity_overlaps[9]==[464],'seven/nine-term overlap controls')
    for n in range(1,5001):
        ds=divisors(n);need(sum(ds)==sig[n]==sigma_factor(n) and len(ds)==tau[n],'independent divisors and prime product')
        rr=[(a,n//a) for a in ds if a<=n//a and 1+a+n//a==s[n]]
        need(rr==([] if roots(n,s[n]) is None else [roots(n,s[n])]),'independent pair search')
    need(roots(16,s[16]) is None and s[16]==15,'explicit sixteen boundary')
    need(all(s[n]==sum(divisors(n)[:-1]) for n in [12,18,20,24,28]),'include proper divisor one')
    examples=[]
    for n in [4,9,25,8,49,6,10,15,14,16,21,35,22,12,18,20,28,24,104,464,650,1952,1001,343,57]:
        examples.append({'n':n,'s':s[n],'factorization':factor(n),'proper_divisors':divisors(n)[:-1],
                         'proper_divisor_count':tau[n]-1,'factor_pair':roots(n,s[n]),
                         'excess':s[n]-n,'mu_lambda':signs(n)})
    row=[]
    for p in [2,3,5,7,11]:
        for a in range(1,6):
            for q in [2,3,5,7,11,13,17,19,23,29,31]:
                if p==q:continue
                n=p**a*q;r=(p**a-1)//(p-1)
                need(sigma_factor(n)-n==1+r*(p+q),'prime-power times prime affine chart')
    for a in range(2,16):
        q=2**(a+1)-3
        if len(factor(q))!=1 or next(iter(factor(q).values()))!=1:continue
        n=2**a*q;need(sigma_factor(n)-2*n==2,'redundant-two classical family')
        row.append({'a':a,'prime':q,'n':n,'proper_divisors':2*a+1})
    need(s[650]==652 and factor(650)=={2:1,5:2,13:1},'family is not exhaustive')
    need([n for n in compressed if s[n]==n]==[6],'perfect overlap finite control')
    inverse343=[n for n in range(1,342**2+1) if s[n]==343]
    fibres={v:[n for n in range(1,(v-1)**2+1) if s[n]==v] for v in [6,8,15,16,21,22,28,36]}
    need(fibres[15]==[16,33] and fibres[22]==[20,38],'same aliquot row different divisor representation')
    triples=[n for n in inverse343 if len(factor(n))==3 and all(e==1 for e in factor(n).values())]
    need(len(inverse343)==19 and triples==[1001],'complete bounded predecessor recovery')
    pair_cases={}
    for p in [3,5,7]:
        v=343+p*(p+1)
        pair_cases[p]=[[d,v//d] for d in range(1,isqrt(v)+1) if v%d==0 and d*d<v]
    orbits=[]
    for seed in [1001,1000001]:
        n=seed;word=[];seen=set()
        while n:
            need(n not in seen,'terminating certified short orbit');seen.add(n);word.append(n)
            v=sigma_factor(n)-n;need(v==s[n],'orbit independent formulas');n=v
        word.append(0);orbits.append(word)
    need(orbits==[[1001,343,57,23,1,0],[1000001,10003,1437,483,285,195,141,51,21,11,1,0]],'inherited orbit controls')
    square_pp=[n for n in sorted(powers) if n>1 and isqrt(sig[n])**2==sig[n]]
    need(square_pp==[3,81,343],'finite prime-power square-sigma control')
    need(1001+1001000==1001**2 and 343+s[343]==20**2,'two different square-sum edges')
    need(isqrt(sig[1001])**2!=sig[1001],'1001 to343 is not square-sum')
    need(signs(1001)==(-1,-1) and signs(343)==(0,-1),'Mobius detects exponent overlap')
    output={'universe_max_source':N,'independent_divisor_range':5000,'class_counts':dict(classes),
            'sigma_sha256_uint32le':digest(sig),'tau_sha256_uint32le':digest(tau),
            'compressed_count':len(compressed),'unit_compression_sources':unit,
            'five_proper_divisor_overlap':overlap,'examples':examples,'classical_redundant_two_family':row,
            'finite_small_arity_overlaps':arity_overlaps,
            'complete_F343':inverse343,'three_distinct_prime_F343':triples,'factor_pair_cases':pair_cases,
            'complete_small_fibres':fibres,
            'orbits':orbits,'finite_prime_power_square_sigma_sources':square_pp,
            'scope':'General classification and unique20 use proofs in THEORY.md; prime-cube square uniqueness imports Ljunggren. No universal aliquot orbit theorem.'}
    print(json.dumps(output,sort_keys=True,indent=2))

if __name__=='__main__':main()
