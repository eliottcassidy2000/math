"""Exact periodic-template parser and phase-corrected integer rank; stdlib only."""
from itertools import product
from pathlib import Path
import json

def check(ok, message):
    if not ok: raise RuntimeError(message)
def v2(x):
    if not x: raise ValueError('valuation of zero')
    x=abs(x)
    return (x & -x).bit_length()-1

def atlas(word):
    L=len(word); A=sum(word); D=2**A-3**L
    B=[]
    for i in range(L):
        rot=word[i:]+word[:i]; b=0; s=0
        for a in rot:
            b=3*b+2**s; s+=a
        B.append(b)
    S=[sum(word[:i]) for i in range(L)]
    e=[L*s-A*i for i,s in enumerate(S)]
    c=[3**(x-min(e)) for x in e]
    return L,A,D,B,c

def predicted_match(word,n):
    L,A,D,B,c=atlas(word); z=D*n-B[0]
    if z==0:return None
    k=v2(z); q,rem=divmod(k-1,A)
    j=q*L; s=0
    for a in word:
        if s+a>rem:break
        s+=a;j+=1
    return j

def energy(data,i,n):
    L,A,D,B,c=data;z=D*n-B[i]
    if z==0:return None
    k=v2(z);t=abs(z)//2**k
    return c[i]*3**(L*k)*t**A,k

def replay(word,n):
    j=0;x=n
    while True:
        a=v2(3*x+1)
        if a!=word[j%len(word)]:return j,x,a
        x=(3*x+1)//2**a;j+=1
        if j>10000:raise RuntimeError('unexpected replay limit')

def main():
    checks=0; words=0; sources=0; fixed=0; maxmatch=0
    for L in range(1,6):
      for tup in product(range(1,5),repeat=L):
        w=list(tup);data=atlas(w);_,A,D,B,c=data;words+=1
        for i,a in enumerate(w):
            check(2**a*B[(i+1)%L]==3*B[i]+D,'cyclic carry');checks+=1
        mod=2**(A+1);r=((2**A-B[0])*pow(3**L,-1,mod))%mod
        candidates={r+mod*t for t in range(3)}|set(range(1,70,2))
        for n in candidates:
            if n<=0:continue
            expected=predicted_match(w,n)
            if expected is None:
                x=n
                for a in w:
                    check(v2(3*x+1)==a,'fixed point native');x=(3*x+1)//2**a
                check(x==n,'fixed point closes');fixed+=1;continue
            j,x,actual=replay(w,n);check(j==expected,'exact mismatch index');checks+=1;sources+=1;maxmatch=max(maxmatch,j)
            x=n
            for h in range(j):
                i=h%L;a=w[i];old=energy(data,i,x)
                x=(3*x+1)//2**a;new=energy(data,(i+1)%L,x)
                check(new[0]==old[0] and new[1]==old[1]-a,'rank transport');checks+=1
            k=energy(data,j%L,x)[1];a=w[j%L]
            check((k<a and actual==k) or (k==a and actual>a),'mismatch direction');checks+=1
        # Repeated words do not manufacture a new rank on the starting phase.
        repeated=w+w;dat2=atlas(repeated);H=2**A+3**L
        check(dat2[2]==D*H and dat2[3][0]==B[0]*H,'power redundancy');checks+=1
    long_controls=[]
    for w in ([1],[2],[1,1,3],[1,2,1,4],[4,3,2,1]):
      for repeats in (10,50,200):
        big=w*repeats;L=len(big);A=sum(big);b=0;s=0
        for a in big:b=3*b+2**s;s+=a
        mod=2**(A+1);n=((2**A-b)*pow(3**L,-1,mod))%mod+mod
        expected=predicted_match(w,n);actual=replay(w,n)[0]
        check(expected==actual and actual>=L,'long shadow');checks+=1
        long_controls.append({'word':w,'repeats':repeats,'match':actual,'source_bits':n.bit_length()})
    result={'long_controls':long_controls,'status' :'PASS','words':words,'sources':sources,'checks':checks,'fixed_point_controls':fixed,'maximum_matched_letters':maxmatch,
      'example':{'word':[1,1,3],'source':55,'predicted_matched_letters':predicted_match([1,1,3],55),'replay':replay([1,1,3],55)}}
    print(json.dumps(result,indent=2))
if __name__=='__main__':main()
