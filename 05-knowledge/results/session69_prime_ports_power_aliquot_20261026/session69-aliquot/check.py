#!/usr/bin/env python3
"""Session69: exact aliquot-fiber and radical interfaces; explicit gates survive -O."""
import argparse, hashlib, json, sys
from math import gcd,isqrt
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent/'source'))
from aliquot_fibers import fiber

def demand(ok,msg):
    if not ok: raise RuntimeError(msg)
def primesieve(N):
    p=bytearray(b'\x01')*(N+1); p[0:2]=b'\x00\x00'
    for d in range(2,isqrt(N)+1):
        if p[d]: p[d*d:N+1:d]=b'\x00'*(((N-d*d)//d)+1)
    return p

def sigma_trial(n):
    s=0
    for d in range(1,isqrt(n)+1):
        if n%d==0:
            s+=d
            if d*d!=n:s+=n//d
    return s

def radical(n):
    r=1; p=2
    while p*p<=n:
        if n%p==0:
            r*=p
            while n%p==0:n//=p
        p+=1
    return r*n

def normalize(A,B):
    g=gcd(A,B); return (A//g,B//g,(A+B)//g,g)

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--negative-control',action='store_true');args=ap.parse_args()
    M=1000;B=(M-1)**2
    # Independent source-sieve path. Complete by n <= (m-1)^2 for m>1.
    sig=[0]*(B+1)
    for d in range(1,B//2+1):
        for n in range(2*d,B+1,d):sig[n]+=d
    expected={m:[] for m in range(2,M+1)}
    for n,m in enumerate(sig):
        if 2<=m<=M:expected[m].append(n)
    actual={m:fiber(m) for m in range(2,M+1)}
    demand(actual==expected,'inherited inverse algorithm disagrees with divisor sieve')
    demand(fiber(0)==[1],'zero fiber')
    try:fiber(1)
    except ValueError:pass
    else:raise RuntimeError('infinite prime fiber was represented as finite')
    target6=[6] if args.negative_control else [6,25]
    demand(actual[6]==target6,'negative-control target6 completeness')
    demand(actual[2]==[] and actual[5]==[],'small untouchables')
    demand(actual[3]==[4] and actual[4]==[9],'transcription control')
    demand(next(m for m in range(2,M+1) if len(actual[m])>1)==6,'first finite collision')
    demand(next(n for n in range(1,100) if sigma_trial(n)==2*n)==6,'first perfect')
    demand(actual[7]==[8] and actual[21]==[18,51,91],'H-gap touchability')
    # Independent witnesses for every prime target <= 10^6 except 2,5.
    X=1000000; isp=primesieve(X);primes=[p for p in range(2,X+1) if isp[p]]
    witnesses=[];special={3:4,7:8};maxfirst=0
    for m in primes:
        if m in (2,5):continue
        if m in special:
            n=special[m];demand(sigma_trial(n)-n==m,'special prime witness');witnesses.append((m,n));continue
        pair=None
        for p in primes:
            q=m-1-p
            if p>=q:break
            if isp[q]:pair=(p,q);break
        demand(pair is not None,f'no distinct-Goldbach witness for {m}')
        p,q=pair;demand(isp[p] and isp[q] and p!=q and 1+p+q==m,'bad semiprime witness')
        witnesses.append((m,p*q));maxfirst=max(maxfirst,p)
    # Affine prime-extension ABC normalization, including exceptional primes dividing T/d.
    branch_cases=0;generic=0;exceptional=0;already_supported=0
    for a in range(2,201):
        T=sigma_trial(a);S=T-a;d=gcd(S,T);U=S//d;V=T//d
        for p in primes:
            if p>199:break
            if a%p==0:continue
            branch_cases+=1;m=S*p+T
            demand(sigma_trial(a*p)-a*p==m,'guarded prime-extension identity')
            h=gcd(p,V);A=U*p//h;BB=V//h;C=(U*p+V)//h
            demand((A,BB,C,d*h)==normalize(S*p,T),'affine gcd normal form')
            demand(gcd(A,BB)==1 and A+BB==C,'primitive ABC branch')
            R=radical(A*BB*C)
            demand(C<=(U+V)*R,'all-branch fixed-core radical bound')
            if h==1:
                generic+=1
                demand(R==p*radical(U*V*C)//gcd(p,U),'generic radical factor with support overlap')
                if U%p==0:already_supported+=1
                demand(C<=(U+V)*R,'unconditional fixed-core exponent-one bound')
            else:
                exceptional+=1;demand(h==p and V%p==0,'exceptional branch prime')
    # Twin-prime collisions: exact witnesses, not an asserted entire fiber.
    twins=[]
    for p in primes:
        if p>997:break
        if p>2 and isp[p+2]:
            n1=2*p;n2=(p+2)**2;m=p+3
            demand(n1!=n2 and sigma_trial(n1)-n1==m==sigma_trial(n2)-n2,'twin-prime collision')
            demand((n1==m)==(p==3),'self-preimage unique within twin family')
            demand(n1 in actual[m] and n2 in actual[m],'complete-fiber twin witnesses')
            twins.append([p,m,n1,n2])
    demand(sigma_trial(27)-27==13 and sigma_trial(64)-64==63!=21,'composite prime-power-base hostile')
    for q in range(2,31):
        demand(1+q*q+q**4==(1+q+q*q)*(1-q+q*q),'cyclotomic squaring identity')
    # Genuine/false bridges and exact branch boundaries.
    demand(sigma_trial(4)-4==3 and (sigma_trial(2)-2)*2+sigma_trial(2)==5,'noncoprime hostile')
    demand(sigma_trial(9)-9==4 and 1+3+3==7,'equal-prime Goldbach hostile')
    demand(sigma_trial(2)-2==sigma_trial(3)-3==1 and sigma_trial(4)-4!=sigma_trial(6)-6,'fiber quotient not multiplicative')
    demand(normalize(6,6)[:3]==normalize(28,28)[:3]==(1,1,2),'perfect scale-loss hostile')
    demand(normalize(6,6)[:3]!=normalize(25,6)[:3],'same aliquot target does not fix primitive triple')
    demand(sigma_trial(8)-8==1+2+4==7,'1-2-4 divisor bridge')
    demand(sigma_trial(32)-32==31 and 32==2**2*2**3,'squarefull gcd hostile')
    # Small labelled tournament positive control H=5, independently brute counted.
    import itertools
    verts=range(4);edges=list(itertools.combinations(verts,2));five=[]
    for mask in range(1<<len(edges)):
        arcs=set((u,v) if mask>>i&1 else (v,u) for i,(u,v) in enumerate(edges))
        words=[p for p in itertools.permutations(verts) if all((p[i],p[i+1]) in arcs for i in range(3))]
        if len(words)==5:five.append(mask)
    demand(bool(five),'untouchable5 is nevertheless raw H-attainable')
    payload={'status':'PASS','source_bound':B,'complete_targets':[2,M],
      'small_complete_fibers':{'0':[1],'1':'all primes, infinite',**{str(m):actual[m] for m in range(2,65)}},
      'prime_target_witness_cutoff':X,'prime_targets':len(primes),'prime_targets_certified_touchable':len(witnesses),
      'excluded_prime_targets':[2,5],'largest_first_goldbach_prime':maxfirst,
      'witness_list_sha256':hashlib.sha256(json.dumps(witnesses,separators=(',',':')).encode()).hexdigest(),
      'abc_branch_core_range':[2,200],'abc_branch_prime_bound':199,'abc_branch_cases':branch_cases,
      'generic_noncancelled_prime_branches':generic,'prime_already_in_core_radical':already_supported,'exceptional_divisor_primes':exceptional,
      'twin_prime_pairs_tested':len(twins),'first_twin_fiber_witnesses':twins[:8],'labelled_order4_tournaments_with_H5':len(five),'first_mask_H5':five[0],
      'controls':['F(1) infinite','F(3)/F(4) typo','first finite collision F(6)',
      '7 and 21 touchable;5 H-attainable','p divides core','equal Goldbach primes',
      'squarefull may have gcd(u,v)>1','fiber quotient nonmultiplicative',
      'ABC primitive normalization loses perfect scale','fiber equality does not identify ABC triples','composite base4 destroys prime-cube divisor list','radical overlap when p divides S/d']}
    print(json.dumps(payload,indent=2,sort_keys=True))
if __name__=='__main__':main()
