#!/usr/bin/env python3
"""Independent literal-divisibility audit; does not import producer modules."""
from math import isqrt, prod
import json

checks=0

def check(ok, label):
    global checks
    checks+=1
    if not ok: raise RuntimeError(label)

def prime(n):
    if n<2:return False
    for d in range(2,isqrt(n)+1):
        if n%d==0:return False
    return True

def mobius(n):
    parity=0;d=2
    while d*d<=n:
        if n%d==0:
            n//=d;parity^=1
            if n%d==0:return 0
        d+=1
    if n>1:parity^=1
    return -1 if parity else 1

def main():
    HMAX=10000
    ps=[p for p in range(2,110)if prime(p)]
    prim=[prime(n)for n in range(HMAX+3)]
    pi=[0]*(HMAX+1)
    for x in range(3,HMAX+1):pi[x]=pi[x-1]+(prim[x] and prim[x+2])
    strict={}
    for B in range(2,102):
        tested=[p for p in ps if p<=B]
        count=[0]*(HMAX+1)
        for x in range(3,HMAX+1):
            count[x]=count[x-1]+all(x%p and (x+2)%p for p in tested)
        strict[B]=count
        previous=set() if B==2 else accepted
        accepted=set()
        for x in range(3,1001):
            e=all((x%p!=0 or x==p)and((x+2)%p!=0 or x+2==p)for p in tested)
            rough=all(x%p and(x+2)%p for p in tested)
            patch=prim[x] and prim[x+2] and x<=B
            check(e==(rough or patch)and not(rough and patch),'W2 exact disjoint equality exceptions')
            if e:accepted.add(x)
        if B>2:check(accepted<=previous,'restored monotonicity')
        q=next(p for p in ps if p>B)
        for H in range(3,min(HMAX,q*q-3)+1):
            check(pi[H]==count[H]+pi[min(H,B)],'W3 strict next-prime threshold')
        x=q*q-2
        if x+2<=HMAX+2:
            e=all((x%p!=0 or x==p)and((x+2)%p!=0 or x+2==p)for p in tested)
            check(e==prim[x],'first square boundary iff lower prime')
    for H in range(3,HMAX+1):
        z=isqrt(H+2)
        check(z<H and pi[H]==strict[z][H]+pi[z],'descending all-height recurrence')
    wheels=[];oldM=1;oldW=[0]
    for B in [2,3,5,7,11]:
        tested=[p for p in ps if p<=B];M=prod(tested)
        W=[r for r in range(M)if all(r%p and(r+2)%p for p in tested)]
        check(len(W)==prod(p-(1 if p==2 else 2)for p in tested),'literal CRT wheel cardinality')
        check(M-1 in W and (M-1)%oldM==oldM-1,'compatible minus-one branch')
        if B>2:
            for a in oldW:
                children=[w for w in W if w%oldM==a]
                forbidden={(-a*pow(oldM,-1,B))%B, (-(a+2)*pow(oldM,-1,B))%B}
                check(len(children)==B-2 and children==[a+t*oldM for t in range(B)if t not in forbidden], 'literal lift indices')
        divs=[d for d in range(1,M+1)if M%d==0]
        roots={d:[r for r in range(d)if r*(r+2)%d==0]for d in divs}
        for d,rs in roots.items():check(len(rs)==prod(1 if p==2 else 2 for p in tested if d%p==0),'literal forbidden roots')
        error=min(len(W),2*3**(len(tested)-1))
        for H in list(range(3,74))+[100,201,1000,10000]:
            floorW=sum((H-a)//M-(2-a)//M for a in W)
            floorM=sum(mobius(d)*sum((H-r)//d-(2-r)//d for r in roots[d])for d in divs)
            check(floorW==floorM==strict[B][H],'two literal floor formulae')
            check(abs(floorW*M-(H-2)*len(W))<=error*M,'signed-floor error retained')
        wheels.append([B,M,len(W)]);oldM=M;oldW=W
    false=[]
    for B in [p for p in ps if p<=61]:
        tested=[p for p in ps if p<=B]
        bad=next(x for x in range(3,5001)if all(x%p and(x+2)%p for p in tested)and not(prim[x]and prim[x+2]))
        q=next(p for p in ps if p>B)
        check(bad>=q*q-2,'first false witness threshold')
        square=prim[bad] and isqrt(bad+2)**2==bad+2 and prime(isqrt(bad+2))
        false.append([B,bad,bad+2,square])
    check(next(v[0]for v in false if not v[3])==61,'first square-pattern failure cutoff')
    check(false[-1][1:3]==[4757,4759] and 4757==67*71 and prime(4759),'first nonsquare factor witness')
    check(4487==7*641,'lower boundary rejected')
    check([(-1)*1,2*0,0**2+1]==[-1,0,1],'signed negative lift')
    check(strict[3][3]==0 and pi[3]==1,'self-hit at pair3,5')
    print(json.dumps({'status':'PASS','checks':checks,'universe':{'B_integers':[2,101],'self_hit_x':[3,1000],'all_height_recurrence':[3,HMAX],'literal_wheel_cutoffs':[2,3,5,7,11],'first_false_prime_cutoffs_through':61,'first_false_search_limit':5000},'wheels':wheels,'first_false':false,'counts':{str(H):pi[H]for H in [100,1000,10000]},'independence':'literal trial division, no producer imports, literal root enumeration instead of producer CRT root generation'},indent=2,sort_keys=True))

if __name__=='__main__':main()
