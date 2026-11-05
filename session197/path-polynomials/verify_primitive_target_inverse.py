#!/usr/bin/env python3
from verify_primitive_cancellation import *
import random
D0=-C0;H0=3*A0-D0

def count(t,s,j,m):
    return A0*4**t*22352**(2*s)*352**(2*j)*3**m-D0*1147**(2*s)*37**(2*j)

def inverse(N):
    if type(N) is not int or N<H0:return []
    recipes=[];s=0;xs=1;hs=1
    while H0*xs<=N:
        j=0;ys=1;hy=1
        while H0*xs*ys<=N:
            a=A0*xs*ys;d=D0*hs*hy
            if (N+d)%a==0:
                q=(N+d)//a;e2=e3=0
                while q%2==0:q//=2;e2+=1
                while q%3==0:q//=3;e3+=1
                if q==1 and e2%2==0 and e3>=1:
                    t=e2//2;m=e3
                    require(count(t,s,j,m)==N,'inverse algebraic witness')
                    recipes.append((t,s,j,m))
            j+=1;ys*=352**2;hy*=37**2
        s+=1;xs*=22352**2;hs*=1147**2
    return recipes

def bounded_table(cap):
    out={};s=0;xs=1
    while H0*xs<=cap:
        j=0;ys=1
        while H0*xs*ys<=cap:
            t=0;kt=1
            while H0*xs*ys*kt<=cap:
                m=1;km=1
                while H0*xs*ys*kt*km<=cap:
                    N=count(t,s,j,m)
                    if N<=cap:out.setdefault(N,[]).append((t,s,j,m))
                    m+=1;km*=3
                t+=1;kt*=4
            j+=1;ys*=352**2
        s+=1;xs*=22352**2
    return out

def run():
    positive=[5990037973,161868952021]
    require(inverse(positive[0])==[(0,0,0,1)],'exact prime recipe')
    require(inverse(positive[1])==[(0,0,0,4)],'exact7p recipe')
    missing=5990037979;require(inverse(missing)==[],'all-order missing-prime certificate')
    cert=prime_certificate(missing);recoveries=0
    for t,s,j,m in product(range(3),range(3),range(3),range(1,9)):
        N=count(t,s,j,m)
        require((t,s,j,m) in inverse(N),'every actual finite-box recipe recovered')
        recoveries+=1
    cap=10**13;table=bounded_table(cap);rng=random.Random(1971150)
    targets=set(table)
    for N in table:
        targets.update(n for n in (N-2,N+2,N-6,N+6) if 0<n<=cap)
    targets.update(rng.randrange(1,cap+1) for _ in range(500))
    targets.update((0,-1,1,H0-1))
    for N in targets:
        got=inverse(N)
        require(sorted(got)==sorted(table.get(N,[])),'complete independent bounded-table comparison')
    for bad in (None,True,1.5,'5990037973'):
        require(inverse(bad)==[],'noninteger target rejected')
    print(json.dumps({'status':'PASS','positive_targets':[{'N':N,'recipes':inverse(N)} for N in positive],'absent_prime':cert,'finite_box_recipes_recovered':recoveries,'complete_table_cap':cap,'complete_table_members':len(table),'tested_targets':len(targets)},sort_keys=True,indent=2))
if __name__=='__main__':run()
