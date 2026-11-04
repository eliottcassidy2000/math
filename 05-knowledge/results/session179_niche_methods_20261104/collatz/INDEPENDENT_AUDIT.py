"""Independent audit: symbolic Fractions, exact precision cylinders, no producer imports."""
from fractions import Fraction
from itertools import product
import json

def need(ok,what):
    if not ok: raise RuntimeError(what)
def valuation(x):
    if not x: raise ValueError('zero')
    k=0
    while x%2==0: x//=2;k+=1
    return k

def parameters(w):
    L=len(w);A=sum(w);D=2**A-3**L; Bs=[]
    for j in range(L):
        slope=Fraction(1);intercept=Fraction(0)
        for a in w[j:]+w[:j]:
            slope=3*slope/2**a
            intercept=(3*intercept+1)/2**a
        need(slope==Fraction(3**L,2**A),'slope')
        B=intercept*2**A;need(B.denominator==1,'carry integer');Bs.append(int(B))
    exponents=[L*sum(w[:i])-A*i for i in range(L)]
    return D,Bs,[3**(e-min(exponents)) for e in exponents]

def measure(w,pars,i,n):
    D,B,c=pars;z=D*n-B[i];k=valuation(z)
    E=c[i]*3**(len(w)*k)*abs(z//2**k)**sum(w)
    need(E>=abs(z),'properness bound')
    return E,k

def audit_case(w,K):
    pars=parameters(w);D,B,c=pars
    # Force exact K, independently of affine cylinder generation in producer.
    mod=2**(K+1)
    n=((B[0]+2**K)*pow(D,-1,mod))%mod+mod
    need(n>0 and n%2==1 and valuation(D*n-B[0])==K,'exact cylinder')
    cumulative=0; predicted=0
    while cumulative+w[predicted%len(w)]<K:
        cumulative+=w[predicted%len(w)];predicted+=1
    x=n;j=0
    while True:
        actual=valuation(3*x+1);a=w[j%len(w)]
        if actual!=a:
            need(j==predicted,'mismatch')
            remaining=valuation(D*x-B[j%len(w)])
            need((remaining<a and actual==remaining) or (remaining==a and actual>a),'mismatch direction')
            break
        old=measure(w,pars,j%len(w),x)
        x=(3*x+1)//2**actual;j+=1
        new=measure(w,pars,j%len(w),x)
        need(old[0]==new[0] and new[1]<old[1],'rank step')
    for r in [2,3]:
        rp=parameters(w*r);A=sum(w);L=len(w)
        H=sum(2**(A*(r-1-h))*3**(L*h) for h in range(r))
        need(rp[0]==H*D and rp[1][0]==H*B[0] and H%2==1,'repeat carry')
        E,k=measure(w,pars,0,n);Er,kr=measure(w*r,rp,0,n)
        need(k==kr and Er==H**(r*A)*E**r,'repeat energy exact constant')
    return predicted

def main():
    cases=0;maxrun=0;positiveD=negativeD=0
    for L in range(1,5):
        for tup in product(range(1,5),repeat=L):
            w=list(tup);D,B,c=parameters(w)
            positiveD+=int(D>0);negativeD+=int(D<0)
            for K in [1,2,3,4,5,7,11,19,31]:
                maxrun=max(maxrun,audit_case(w,K));cases+=1
    long=[]
    for w in [[1],[2],[1,1,3],[1,2,1,4],[4,3,2,1]]:
        K=sum(w)*205+1
        matched=audit_case(w,K);need(matched==len(w)*205,'long exact boundary');long.append(matched)
    # Fixed-center and negative signed-cofactor controls.
    for L in range(1,7):
        D,B,c=parameters([2]*L);need(all(b==D for b in B),'fixed center1')
    w=[1,1,3];D,B,c=parameters(w);need(D*1-B[0]<0,'negative cofactor exists')
    E,k=measure(w,(D,B,c),0,1);need(k==1 and valuation(4)==2,'negative cofactor mismatch')
    # Neutral tower arithmetic and old claimed lower-energy family.
    towers=0
    for r in range(1,12):
        modulus=3**r;t=(-232*pow(65536,-1,modulus))%modulus
        n=223+62208*t;m=233+65536*t
        mr=1+4**r*((m-1)//modulus)
        old=measure([2],parameters([2]),0,m)[0]
        need(old<measure([2],parameters([2]),0,n)[0],'old family energy decrease')
        need(old==measure([2],parameters([2]),0,mr)[0],'tower neutrality')
        x=mr
        for h in range(r): need(valuation(3*x+1)==2,'tower native');x=(3*x+1)//4
        need(x==m,'tower returns');towers+=1
    print(json.dumps({'status':'PASS','precision_cylinder_cases':cases,'positive_D_words':positiveD,'negative_D_words':negativeD,'largest_short_match':maxrun,'long_exact_matches':long,'neutral_towers':towers,'repeated_powers':[2,3]},indent=2))
if __name__=='__main__':main()
