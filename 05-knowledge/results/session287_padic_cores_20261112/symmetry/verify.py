#!/usr/bin/env python3
"""Exact deterministic controls for symmetry/PROOF.md; no external packages."""
import itertools,json,math
from fractions import Fraction as F
from collections import Counter

def check(x,msg):
    if not x:raise RuntimeError(msg)

def mul(a,b):
    c=[0]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        if x:
            for j,y in enumerate(b):
                if y:c[i+j]+=x*y
    return c

def step(a):
    c=mul(mul(a,a),a)
    c += [0]*(3*len(a)-1-len(c))
    for i,x in enumerate(a):c[3*i]-=x;c[3*i+1]+=3*x
    check(all(x%3==0 and x>=0 for x in c),'orbit polynomial integers')
    while c[-1]==0:c.pop()
    return [x//3 for x in c]

def moments(a,e,q):
    I=(3**e-1)//2; N=q**(3**e)
    lab=[x*3**(I-i) for i,x in enumerate(a)]
    check(sum(lab)==N,'labelled mass')
    mean=F(sum(i*x for i,x in enumerate(lab)),N)
    var=F(sum(i*i*x for i,x in enumerate(lab)),N)-mean**2
    return lab,mean,var

def collision(a,e,q):
    I=(3**e-1)//2;N=q**(3**e)
    alpha=sum(F(x*3**(3*(I-i)),N**3) for i,x in enumerate(a))
    beta=sum(F(i*x*3**(3*(I-i)),N**3) for i,x in enumerate(a))
    return alpha,beta

def bounds(q,r,mean,var):
    rho=F(3,q*q); t=rho**(3**r)
    dm=t/(3**(r+2)*(1-t*t/3));dv=t/(3*(1-t*t))
    return (mean/3**r,mean/3**r+dm),(var/3**r-dv,var/3**r+dv)

def interval_string(lo,hi,d=18):
    s=10**d
    low=(lo*s).numerator//(lo*s).denominator
    high=-((-hi*s).numerator//(-hi*s).denominator)
    def fmt(x):return f'{x//s}.{x%s:0{d}d}'
    return [fmt(low),fmt(high)]

def tower(e):
    n=3**e; digs=[tuple((i//3**j)%3 for j in reversed(range(e))) for i in range(n)]
    return [sum(1<<j for j in range(n) if i!=j and next((b-a)%3 for a,b in zip(digs[i],digs[j]) if a!=b)==1) for i in range(n)]

def group(e):
    if not e:return [(0,)]
    old=group(e-1);m=3**(e-1);out=[]
    for r in range(3):
        for trip in itertools.product(old,repeat=3):
            out.append(tuple(((b+r)%3)*m+trip[b][i] for b in range(3) for i in range(m)))
    return out

def H(rows):
    n=len(rows);dp=[[0]*n for _ in range(1<<n)]
    for i in range(n):dp[1<<i][i]=1
    for mask in range(1,1<<n):
        for v,c in enumerate(dp[mask]):
            if c:
                rem=rows[v]&~mask
                while rem:
                    bit=rem&-rem;rem-=bit;w=bit.bit_length()-1
                    dp[mask|bit][w]+=c
    return sum(dp[-1])

def valuation(x,p=3):
    check(x>0,'positive valuation input');s=0
    while x%p==0:x//=p;s+=1
    return s

def main():
    out={'status':'PASS','fixed_q':{},'finite_scopes':{}}
    polys={}
    for q in (2,3,4,8):
        a=[q];prev_mean=F(0);prev_var=F(0);records=[]
        for e in range(1,6):
            alpha,beta=collision(a,e-1,q);old=a;a=step(a)
            lab,mean,var=moments(a,e,q);rho=F(3,q*q)
            check(alpha<=rho**(3**(e-1))/3,'collision upper bound')
            check(mean==3*prev_mean+alpha,'mean recurrence')
            delta=alpha*(1-alpha)+6*(beta-prev_mean*alpha)
            check(var==3*prev_var+delta,'variance recurrence with covariance')
            check(abs(delta)<3**e*alpha,'variance increment bound')
            # Independent labelled-pgf identity with exact rational coefficients.
            oldlab,_,_=moments(old,e-1,q)
            pred=mul(mul(oldlab,oldlab),oldlab)+[0]
            I=(3**(e-1)-1)//2
            for j,x in enumerate(oldlab):
                term=F(3**(2*I)*x,9**j)
                check(term.denominator==1,'integral labelled correction')
                pred[3*j]-=term.numerator;pred[3*j+1]+=term.numerator
            while pred[-1]==0:pred.pop()
            check(pred==lab,'independent labelled recurrence')
            if q==2:check(a[0]==2,'inherited two binary free orbits')
            records.append({'depth':e,'n':3**e,'mean_density':float(mean/3**e),'variance_density':float(var/3**e),'free_orbits':str(a[0])})
            if q==2 and e==3:
                check(var/27>F(56,1000),'depth3 variance positive control')
                check(rho**27/(3*(1-rho**54))<F(1,1000),'q2 tail positivity control')
            prev_mean,prev_var=mean,var
            if e==2:polys[q]=lab
        mb,vb=bounds(q,5,mean,var)
        check(vb[0]>0,'positive limiting variance certificate')
        out['fixed_q'][q]={'certified_mu_interval':interval_string(*mb),'certified_sigma2_interval':interval_string(*vb),'depths':records}
    # Symbolic q>=3 variance positivity reduces to monotone factors; boundary exact.
    check(F(3**4)*(1-F(1,9))*(1-F(729,3**12))>27,'q>=3 positivity boundary')
    # Explicit group and all binary colorings of the nine-tower.
    rows=tower(2);g=group(2);check(len(g)==81 and len(set(g))==81,'explicit wreath group')
    for p in g:
        check(all(bool(rows[i]>>j&1)==bool(rows[p[i]]>>p[j]&1) for i in range(9) for j in range(9)),'literal automorphism')
    hist=Counter();hv=Counter();different=0;rootfree=0
    for bits in itertools.product(range(2),repeat=9):
        stab=sum(all(bits[i]==bits[p[i]] for i in range(9)) for p in g)
        s=valuation(stab);hist[s]+=1;rootfree+=s==0
        bottom=sum(bits[j]==bits[j+1]==bits[j+2] for j in (0,3,6))
        check(s>=bottom,'nonnegative ancestor toll');different+=s!=bottom
        # New root beats vertex i when bits[i]=1.
        ext=[row | ((1<<9) if not bits[i] else 0) for i,row in enumerate(rows)]
        ext.append(sum(1<<i for i,b in enumerate(bits) if b))
        h=H(ext);check(h%3==0,'inherited endpoint divisor for every literal extension')
        check(h%stab==0,'root stabilizer free path action');hv[valuation(h)]+=1
    check([hist[i] for i in range(len(polys[2]))]==polys[2],'literal group enumeration matches distribution')
    bound=F(9,2**6)/(1-F(3,4)**6/3)
    check(F(different,512)<=bound,'finite coupling inequality')
    check(rootfree==2*81,'two free orbits literal count')
    out['finite_scopes']['all_binary_attachments_n9']={'colorings':512,'root_stabilizer_valuation_histogram':dict(sorted(hist.items())),'path_count_valuation_histogram':dict(sorted(hv.items())),'nonbottom_toll_probability':str(F(different,512)),'coupling_upper_bound':str(bound),'rootfree_labelled':rootfree}
    out['poisson_controls']=[]
    for e,k in [(8,7),(12,10),(16,13),(20,16),(20,20)]:
        n=3**e;q=2**k;rho=F(3,q*q)
        lam=F(n,3*q*q);err=F(n,q**6)/(1-rho**6/3)+F(2*n,3*q**4)
        E=sum(max(3**d-2*k,0) for d in range(e))
        falling=math.prod(range(n+1,n+k+1));floorlog=0
        while 3**(floorlog+1)<=falling:floorlog+=1
        out['poisson_controls'].append({'e':e,'newborns':k,'poisson_mean':float(lam),'tv_upper_bound':float(err),'universal_H_v3_floor':E,'quotient_v3_floor_if_pointwise_rigid':max(0,E-floorlog)})
    print(json.dumps(out,indent=2,sort_keys=True))
if __name__=='__main__':main()
