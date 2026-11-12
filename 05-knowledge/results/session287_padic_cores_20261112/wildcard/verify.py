#!/usr/bin/env python3
"""Self-contained exact checks; no external dependencies, no assert-based gates."""
import itertools as it, math, json, random

def check(ok,msg):
    if not ok: raise RuntimeError(msg)

def tournament(n,code):
    a=[0]*n
    for k,(i,j) in enumerate(it.combinations(range(n),2)):
        if code>>k&1:a[i]|=1<<j
        else:a[j]|=1<<i
    return a

def edge(a,i,j):return a[i]>>j&1

def hp(a):
    n=len(a);dp=[{} for _ in range(1<<n)]
    for v in range(n):dp[1<<v][v]=1
    for mask in range(1,1<<n):
        for v,c in dp[mask].items():
            bits=a[v]&~mask
            while bits:
                b=bits&-bits;bits-=b;w=b.bit_length()-1
                d=dp[mask|b];d[w]=d.get(w,0)+c
    return sum(dp[-1].values())

def inflate(a,v,m):
    old=[x for x in range(len(a)) if x!=v]; reps=old+[v]*m;n=len(reps);b=[0]*n
    for i,j in it.combinations(range(n),2):
        z=(i<j) if reps[i]==reps[j] else edge(a,reps[i],reps[j])
        if z:b[i]|=1<<j
        else:b[j]|=1<<i
    return b

def skeleton(a,v):
    out=[x for x in range(len(a)) if x!=v];r=len(out)
    c=[0]*(r//2+2)
    for q in range(1,len(c)):
        for perm in it.permutations(out):
            for gaps in it.combinations(range(r+1),q):
                gs=set(gaps);w=[]
                for k in range(r+1):
                    if k in gs:w.append(v)
                    if k<r:w.append(perm[k])
                if all(edge(a,i,j) for i,j in zip(w,w[1:])):c[q]+=1
    while len(c)>2 and c[-1]==0:c.pop()
    return c

def modes(c):
    R=len(c)-1
    return [0]+[sum((-1)**(q-j)*math.comb(q,j)*c[q] for q in range(j,R+1)) for j in range(1,R+1)]

def val(b,m,mod=None,units=None):
    return sum(b[j]*(j**m if mod is None else pow(j,m,mod)) for j in range(1,len(b)) if units is None or j%units)%mod if mod else sum(b[j]*j**m for j in range(1,len(b)))

def stirling(m,q):
    return sum((-1)**(q-j)*math.comb(q,j)*j**m for j in range(q+1))

def windmill(s):
    # root0; ordered pairs (a_i,b_i) with root->a_i->b_i->root.
    n=2*s+1;a=[0]*n
    for i,j in it.combinations(range(n),2):
        z=(j%2==1) if i==0 else True
        if z:a[i]|=1<<j
        else:a[j]|=1<<i
    return a

def isprime(p):return p>=2 and all(p%d for d in range(2,math.isqrt(p)+1))
def vp(n,p):
    if n==0:return 10**9
    e=0
    while n%p==0:n//=p;e+=1
    return e

def parameters(b,p,m0):
    L,delta=(2,3) if p==2 else (p-1,1)
    J=[j for j in range(2,len(b)) if j%p and b[j]]
    if not J:return L,delta,0
    e=min(vp(b[j],p) for j in J)
    D=sum((b[j]//p**e)*pow(j,m0,p)*((j**L-1)//p**delta) for j in J)%p
    return L,delta+e,D

def shadow(b,p,m0,t,mod):
    L=2 if p==2 else p-1
    return sum(b[j]*pow(j,m0+L*t,mod) for j in range(1,len(b)) if j%p)%mod

def inverse(b,p,m0,target,a):
    L,d,D=parameters(b,p,m0)
    check(D!=0,'singular inverse requested');check(a>=d,'bad depth')
    check((target-shadow(b,p,m0,0,p**d))%p**d==0,'outside ball')
    t=0
    for ell in range(a-d):
        mod=p**(d+ell+1)
        choices=[t+x*p**ell for x in range(p) if shadow(b,p,m0,t+x*p**ell,mod)==target%mod]
        check(len(choices)==1,'unique Hensel child');t=choices[0]
    period=p**(a-d);m=m0+L*t
    while m<a:m+=L*period
    check(val(b,m,p**a)==target%p**a,'physical shadow realization')
    return m

def logmod(alpha,p,precision):
    mod=p**precision;z=alpha-1;total=0
    for k in range(1,2*precision+11):
        e=vp(k,p);q=k//p**e
        num=pow(z,k,mod*p**e)//p**e
        total=(total+(-1)**(k+1)*num*pow(q,-1,mod))%mod
    return total

def main():
    results={};count=0;profiles=set();examples={}
    # Exhaustive rooted cores through order4; sample order5 separately.
    for n in range(1,5):
        for code in range(1<<(n*(n-1)//2)):
            a=tournament(n,code)
            for v in range(n):
                c=skeleton(a,v);b=modes(c);profiles.add(tuple(b));count+=1
                for m in range(1,5):
                    h=hp(inflate(a,v,m))
                    check(h==val(b,m)==sum(c[q]*stirling(m,q) for q in range(1,len(c))),'skeleton direct')
                triangle=any(edge(a,v,x) and edge(a,x,y) and edge(a,y,v) for x in range(n) for y in range(n) if len({v,x,y})==3)
                check((len(c)>2)==triangle,'cycle criterion')
    results['rooted_core_controls_through_order4']=count
    rng=random.Random(287);sample=0
    for _ in range(40):
        a=tournament(5,rng.randrange(1024));v=rng.randrange(5);c=skeleton(a,v);b=modes(c);profiles.add(tuple(b))
        for m in (1,2,4,5):check(hp(inflate(a,v,m))==val(b,m),'order5 direct');sample+=1
    results['sampled_order5_direct_checks']=sample
    for s in range(1,4):
        a=windmill(s);c=skeleton(a,0);b=modes(c)
        check(len(c)-1==s+1,'sharp rank');examples[str(s)]={'rows':a,'skeleton':c,'modes':b}
        for m in range(1,5):check(hp(inflate(a,0,m))==val(b,m),'windmill direct')
        profiles.add(tuple(b))
    results['sharp_rank_examples']=examples
    # Complete p-adic finite chart permutations, also p<=R and dyadic.
    charts=[];checks=0
    for b in sorted(profiles):
        for p in (2,3,5,7,11):
            L,_=(2,3) if p==2 else (p-1,1)
            for m0 in range(1,L+1):
                _,d,D=parameters(b,p,m0)
                if not D:continue
                depth=3 if p<=5 else 2;a=d+depth;mod=p**a;base=shadow(b,p,m0,0,p**d)
                image={shadow(b,p,m0,t,mod) for t in range(p**depth)}
                expected={base+p**d*x for x in range(p**depth)}
                check(image==expected,'chart onto');checks+=1
                for target in sorted(expected)[::max(1,len(expected)//12)]:inverse(b,p,m0,target,a);checks+=1
                charts.append({'modes':b,'p':p,'m0':m0,'delta':d,'D':D})
    results['chart_and_physical_inverse_checks']=checks
    results['representative_charts']=charts[:10]
    # No false freezing through cancellation: Vandermonde finite criterion.
    tests=0
    for b in profiles:
        R=len(b)-1
        for p in range(R+1,60):
            if not isprime(p):continue
            live=[j for j in range(1,R+1) if b[j]*((j**(p-1)-1)//p)%p]
            Ds=[sum(b[j]*pow(j,r,p)*((j**(p-1)-1)//p) for j in range(1,R+1))%p for r in range(1,R+1)]
            check((not any(Ds))==(not live),'Vandermonde freezing')
            if live:
                s=len(live)
                for start in range(1,2*p):check(any(sum(b[j]*pow(j,r,p)*((j**(p-1)-1)//p) for j in range(1,R+1))%p for r in range(start,start+s)),'no s dead phases')
            tests+=1
    results['tangent_observability_controls']=tests
    interior_controls=0
    for b in profiles:
        S=[j for j in range(2,len(b)) if b[j]]
        for p in (2,3,5,7,11):
            J=[j for j in S if j%p];L=2 if p==2 else p-1;mod=p**8
            ders=[sum(b[j]*j*pow(j**L,t,mod)*logmod(j**L,p,8) for j in J)%mod for t in range(len(J))]
            check(bool(any(ders))==bool(J),'interior finite derivative witness')
            interior_controls+=1
    check(vp(pow(2,1092)-1,1093)==2,'Wieferich deeper positive chart')
    results['interior_derivative_controls']=interior_controls
    results['Wieferich1093_triangle_chart_depth']=2

    # Dyadic shadow != actual interpolation: triangle has only base2 + constant.
    b=(0,1,1)
    check(shadow(b,2,1,0,256)==1 and val(b,1,256)==3,'nonunit transient hostile')
    check(all(parameters(b,1093,r)[2]==0 for r in (1,2,3)),'base2 Wieferich hostile')
    c5=[sum(1<<((i+d)%5) for d in (1,2)) for i in range(5)]
    bc5=modes(skeleton(c5,0));check(bc5==[0,-5,7,2],'C5 mode formula')
    for m in range(1,8):check(hp(inflate(c5,0,m))==val(bc5,m),'C5 literal DP')
    for a in range(4,11):
        mod=2**a
        shadowset={shadow(bc5,2,r,t,mod) for r in (1,2) for t in range(2**(a-4))}
        check(shadowset=={x for x in range(mod) if x%16 in (1,13)},'C5 closure balls')
        for y in shadowset:inverse(bc5,2,1 if y%16==1 else 2,y,a)
    check([val(bc5,m) for m in range(1,4)]==[15,41,105],'C5 exceptions')
    results['C5_exact_dyadic_closure']={'isolated_counts':[15,41,105],'balls_mod16':[1,13],'modes':bc5}
    results['hostiles']=['triangle dyadic shadow is1, not its actual value3 at m1','base2 tangent freezes at Wieferich prime1093']
    results['status']='PASS'
    print(json.dumps(results,sort_keys=True,indent=2))
if __name__=='__main__':main()
