#!/usr/bin/env python3
"""Weighted martingale energy and legal Collatz word pressure controls."""
from fractions import Fraction as F
from itertools import product
import json
from weighted_ports import moments,berggren

CHECKS=0
def ck(ok,msg):
    global CHECKS
    CHECKS+=1
    if not ok:raise RuntimeError(msg)
def wordpoly(n,t):
    a,b=2,1
    for _ in range(n):a,b=b,b+t*a
    return a
def main():
    energy=[]
    for name,prob in [('uniform',(F(1,3),)*3),('biased',(F(1,2),F(1,3),F(1,6))),('atomic',(F(1),F(0),F(0)))]:
        for k in range(7):
            mu=[F(1)]
            for _ in range(k):mu=[x*p for x in mu for p in prob]
            mm=moments('A',k,[3**k*x for x in mu],2)
            collision=sum(p*p for p in prob)
            expected=1+2*sum((3*collision)**j for j in range(k))
            ck(mm[1]==1 and mm[2]==expected,'weighted martingale second moment')
            ck(expected>=1+2*k,'uniform source minimizes quadratic variation')
            if k in [3,6]:energy.append(dict(law=name,depth=k,second_moment=str(expected)))
    # A genuine arithmetic source flow: branch probabilities proportional to reciprocal hypotenuse.
    levels=[[F(1)]]
    for k in range(1,4):
        cs=[F(1,t[2]) for t in berggren(k)];nextlevel=[]
        for parent,mass in enumerate(levels[-1]):
            z=sum(cs[3*parent:3*parent+3]);nextlevel += [mass*x/z for x in cs[3*parent:3*parent+3]]
        levels.append(nextlevel)
        expected=1+2*sum(3**j*sum(x*x for x in levels[j]) for j in range(k))
        actual=moments('A',k,[3**k*x for x in nextlevel],2)[2]
        ck(actual==expected,'arithmetic source flow energy')
    arithmetic_energy=expected
    # Incoming S94's exact valuation clock, independently expanded at the root.
    z,d,z2,d2,zd=F(1),F(0),F(1),F(0),F(0)
    clock_regression=[]
    for k in range(1,13):
        z,d,z2,d2,zd=z,d+2,z2+2*z*z,d2+2*d*d+F(16,3)*d+6,zd+2*z*d+F(8,3)*z
        vz,vd,cov=z2-z*z,d2-d*d,zd-z*d
        ck(vz==2*k,'critical population variance')
        ck(vd==F(8,3)*k*k*(k-1)+2*k,'division-clock variance')
        ck(cov==F(2,3)*k*(3*k-2),'population-clock covariance')
        c=cov/vd;res=vz-cov*cov/vd
        ck(c==F(3*k-2,4*k*k-4*k+3),'optimal clock coefficient')
        ck(res==F(2*k*(3*k*k+5),3*(4*k*k-4*k+3)) and res>F(k,2),'clock-only residual lower bound')
        if k in [1,3,6,12]:clock_regression.append(dict(depth=k,coefficient=str(c),residual_variance=str(res)))
    def compositions(n):
        if not n:yield ()
        for a in range(1,n+1):
            for rest in compositions(n-a):yield (a,)+rest
    for total in range(1,13):
        for word in compositions(total):
            queue=1+sum(3*(a%2==0)-1 for a in word)
            ck((queue-total-1)%2==0,'clock carries exact queue parity')
            ck((queue-(1+3*total-4*len(word)))%6==0,'clock plus odd-step count carries mod-six queue')
    clock_hostiles=[]
    for word,frontier in [((4,1,1,1),0),((2,2,2,1),6)]:
        P,B,A=1,0,0
        for a in word:P,B,A=3*P,3*B+2**A,A+a
        modulus=2**(A+1);source=((2**A-B)*pow(P,-1,modulus))%modulus
        x=source;path=[x]
        for a in word:
            numerator=3*x+1
            ck(numerator&-numerator==2**a,'native hostile valuation')
            x=numerator//2**a;path.append(x)
        ck(3*sum(a%2==0 for a in word[1:])==frontier and A==7,'same terminal clock different frontier')
        clock_hostiles.append(dict(word=word,source=source,modulus=modulus,carry=B,
                                   division_clock=A,frontier=frontier,path=path))
    pressures=[]
    for n in range(3,19):
        ws=[v for v in product([0,1],repeat=n) if all(not(v[i] and v[(i+1)%n]) for i in range(n))]
        for t in [F(1,2),1,2,6]:
            ck(sum(t**sum(v) for v in ws)==wordpoly(n,t),'weighted cyclic word trace')
        ck(wordpoly(n,6)==3**n+(-2)**n,'Collatz drift pressure roots')
        for t in [1,2,6]:ck(wordpoly(3*n,t)==wordpoly(n,t)**3-3*(-t)**n*wordpoly(n,t),'tripling trace')
        expanding=sum(3**sum(v)>=2**(n-sum(v)) for v in ws)
        drift=sum(F(3**sum(v),2**(n-sum(v))) for v in ws)
        ck(expanding<=drift,'explicit first moment tail bound')
        if n in [5,9,15,18]:pressures.append(dict(n=n,words=len(ws),noncontracting=expanding,mean_drift=str(drift/len(ws))))
    # Proper-clock defect, exact prime-field witnesses and source odd-charge discrepancy.
    defects=[]
    for m,p in [(5,11),(7,29),(11,199),(13,521)]:
        a=next(x for x in range(p) if (x*x-x-1)%p==0 and pow(x,m,p)==1)
        n=3*m
        ck(pow(a,2,p)==(1+pow(a,m+1,p))%p,'same proper-clock input')
        diff=(pow(a,-2,p)-1-pow(a,-(m+1),p))%p
        ck(diff==(-2*pow(a,-1,p))%p and diff!=0,'mirror defect is twice a unit')
        ck(F(9,2**(n-2))==6*F(3,2**(n-1)),'hidden odd-charge drift ratio six')
        ck(wordpoly(n,1)==p*(p*p+3),'prime-base tripling fibre')
        defects.append(dict(base_period=m,prime=p,extended_period=n,hidden_fibre=p*p+3,mirror_defect=diff))
    print(json.dumps(dict(energy=energy,arithmetic_flow_second_moment=str(arithmetic_energy),pressure=pressures,
                         clock_regression=clock_regression,
                         clock_hostiles=clock_hostiles,
                         proper_clock_defects=defects,active_checks=CHECKS),indent=2,sort_keys=True))
if __name__=='__main__':main()
