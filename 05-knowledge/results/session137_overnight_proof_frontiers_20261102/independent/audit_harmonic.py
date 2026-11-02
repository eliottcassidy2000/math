"""Independent exact interval-intersection audit, no producer code imports."""
from fractions import Fraction as F
from math import gcd
import json

def require(c,m):
    if not c:raise RuntimeError(m)

def merge(xs):
    out=[]
    for a,b in sorted(xs):
        if out and a<=out[-1][1]:out[-1]=(out[-1][0],max(b,out[-1][1]))
        else:out.append((a,b))
    return out

def safe(speeds,r):
    old=[(F(0),F(1))]
    for v in speeds:
        new=[((k+r)/v,(k+1-r)/v) for k in range(v)]
        out=[];i=j=0
        while i<len(old) and j<len(new):
            a,b=old[i];c,d=new[j]
            if max(a,c)<=min(b,d):out.append((max(a,c),min(b,d)))
            if b<d:i+=1
            else:j+=1
        old=merge(out)
    return old

def claimed(q,A,r,boundary=False):
    eps=F(1,q)-r
    e=lambda p:p+q*(p in A)
    W=[]
    for a in range(1,q):
        if gcd(a,q)==1:
            p=pow(a,-1,q)
            W.append((F(a,q)-eps/e(p),F(a,q)+eps/e(q-p)))
    if boundary:
        d=next((j for j in range(1,q) if j not in A),q)
        W.extend((F(b,q+d),F(b,q+d)) for b in range(1,q+d) if gcd(b,q+d)==1)
    return merge(W)

def run():
    cases=0;boundary_cases=0;max_intervals=0;summary=[]
    for q in range(3,11):
        count=0
        for mask in range(1<<(q-1)):
            A={i for i in range(1,q) if mask>>(i-1)&1}
            C=list(range(1,q))+[q+i for i in sorted(A)]
            d=next((j for j in range(1,q) if j not in A),q)
            low,high=F(1,q+d),F(1,q)
            for r in [low,(3*low+high)/4,(low+high)/2,(low+3*high)/4,high]:
                actual=safe(C,r);pred=claimed(q,A,r,r==low)
                require(actual==pred,f'formula q={q},A={A},r={r},actual={actual},pred={pred}')
                if r==low:
                    positive=sum(a<b for a,b in actual)
                    isolated=sum(a==b for a,b in actual)
                    phi=lambda n:sum(gcd(a,n)==1 for a in range(1,n))
                    require(positive==phi(q),'positive count')
                    require(isolated==phi(q+d)-2*(gcd(q,d)==1),'isolated count')
                    boundary_cases+=1
                max_intervals=max(max_intervals,len(actual));cases+=1;count+=1
        summary.append({'q':q,'promotion_sets':1<<(q-1),'rational_cases':count})
    r=F(1,7);t=F(3,7);C=[1,2,3,4,6,9]
    owners=[v for v in C if min((v*t)%1,(-v*t)%1)==r]
    require(owners==[2,9],'tied owner hostile')
    iso=[str(a) for a,b in safe([1,2,4],F(1,5)) if a==b]
    require(iso==['1/5','4/5'],'isolated boundary hostile')
    print(json.dumps({'scope':'FINITE-EXACT; q=3..10, every promotion subset, five rational thresholds each',
                      'summary':summary,'cases':cases,'boundary_cases':boundary_cases,
                      'max_components':max_intervals,'controls':{'q5_A14_t3over7_owners':owners,
                      'q3_A1_isolated':iso}},sort_keys=True,indent=2))
if __name__=='__main__':run()
