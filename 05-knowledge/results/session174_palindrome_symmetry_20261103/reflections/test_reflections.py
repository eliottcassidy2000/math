#!/usr/bin/env python3
from fractions import Fraction as F
from itertools import product
from math import isqrt
import json

def check(p,msg):
    if not p:raise RuntimeError(msg)
def trim(p):
    p=list(p)
    while len(p)>1 and p[-1]==0:p.pop()
    return p
def add(a,b):return trim([(a[i] if i<len(a) else 0)+(b[i] if i<len(b) else 0) for i in range(max(len(a),len(b)))])
def neg(a):return [-x for x in a]
def ev(p,x):
    y=0
    for c in reversed(p):y=y*x+c
    return y
def rem(a,b):
    a=list(map(F,a));b=list(map(F,b))
    while len(a)>=len(b) and any(a):
        c=a[-1]/b[-1];k=len(a)-len(b)
        for j in range(len(b)):a[j+k]-=c*b[j]
        a=trim(a)
    return a
def sturm(p):
    s=[list(map(F,p)),[F(i*p[i]) for i in range(1,len(p))]]
    while True:
        r=neg(rem(s[-2],s[-1]))
        if not any(r):break
        s.append(r)
    return s
def variations(s,x):
    v=[ev(p,x) for p in s];v=[z for z in v if z]
    return sum(a*b<0 for a,b in zip(v,v[1:]))
def main():
    checks=0
    for base in range(2,8):
        for length in range(1,6):
            for ds in product(range(base),repeat=length):
                p=ev(ds,base);q=ev(ds[::-1],base)
                check((p-q)%(base-1)==0,'reversal residue')
                if ds==ds[::-1] and length%2==0:check(p%(base+1)==0,'palindrome b+1')
                checks+=1
    reflection_cases=0
    for length in range(1,8):
        for ss in product(range(4),repeat=length):
            A=sum((-1)**(length-j)*s for j,s in enumerate(ss,1))
            check(A==(-1)**(length-1)*ev(ss,-1),'evaluation holonomy')
            for x in [-2,0,7]:
                y=x
                for s in ss:y=s-y
                check(y==(-1)**length*x+A,'affine composition')
            if ss==ss[::-1] and length%2==0:check(A==0,'even palindrome')
            reflection_cases+=1
    spectra=0;Ss=[[2],[0,1]]
    for k in range(2,25):Ss.append(add([0]+Ss[-1],neg(Ss[-2])))
    for m in range(1,25):
        D=[1]
        for k in range(1,m+1):D=add(D,Ss[k])
        for delta in [F(-8,9),F(0),F(1,9),F(8,9)]:
            Q=add(D,[-delta]);s=sturm(Q)
            check(len(s[-1])==1,'simple roots')
            check(ev(Q,-2)!=0 and ev(Q,2)!=0,'endpoints')
            check(variations(s,-2)-variations(s,2)==m,'root count')
            reconstruction={};power={0:F(1)}
            for c in Q:
                for k,v in power.items():reconstruction[k+m]=reconstruction.get(k+m,0)+c*v
                nxt={}
                for k,v in power.items():
                    for e in [k-1,k+1]:nxt[e]=nxt.get(e,0)+v
                power=nxt
            check(reconstruction=={i:F(1)-(delta if i==m else 0) for i in range(2*m+1)},'reconstruction')
            spectra+=1
    BQ=add(Ss[15],[6,6]);bs=sturm(BQ);bv=[variations(bs,x) for x in [-2,2]]
    check(bv==[8,7] and len(bs[-1])==1,'Belphegor count/gcd')
    check(ev(BQ,-2)!=0 and ev(BQ,2)!=0,'Belphegor endpoints')
    check(len(sturm([1,2,1])[-1])==2,'boundary repeated root')
    path=[8,1,15,10,6,3,13,12,4,5,11,14,2,7,9]
    labels=[a+b for a,b in zip(path,path[1:])]
    check(all(isqrt(s)**2==s for s in labels),'square labels')
    check(len(set(path))==15 and sorted(path)==list(range(1,16)),'path universe')
    A=sum((-1)**(len(labels)-j)*s for j,s in enumerate(labels,1))
    check(A==1 and labels!=labels[::-1],'Q15 holonomy/nonpalindromy')
    odd_path=[1,15,10,6]
    check([a+b for a,b in zip(odd_path,odd_path[1:])]==[16,25,16] and len(set(odd_path))==4,'odd hostile')
    check(7417-int(str(7417)[::-1])==270,'7417 defect')
    check(240-int('042')==198,'leading-zero hostile')
    check(all(7417%d for d in range(2,isqrt(7417)+1)),'7417 prime')
    order=next(k for k in range(1,7417) if pow(10,k,7417)==1)
    check(order==2472,'order')
    nmods={str(q):(pow(10,1888529,q)-pow(10,944264,q)-1)%q for q in [3,11,240,7417]}
    check(nmods=={'3':2,'11':8,'240':239,'7417':6663},'giant residues')
    for m in range(1,101):check(int('9'*m+'8'+'9'*m)==10**(2*m+1)-10**m-1,'word family')
    print(json.dumps({'status':'PASS','digit_words':checks,'reflection_words':reflection_cases,'exact_sturm_spectra':spectra,'Q15_labels':labels,'Q15_translation':A,'order_7417_10':order,'giant_residues':nmods,'Belphegor_quotient_coefficients_ascending':BQ,'Belphegor_Sturm_variations_at_minus2_plus2':bv,'Belphegor_unit_circle_roots':2,'Belphegor_off_unit_circle_roots':28,'huge_integer_primality':'NOT TESTED'},indent=2))
if __name__=='__main__':main()
