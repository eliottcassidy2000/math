#!/usr/bin/env python3
"""Finite exhaustive palette checks, independent of inherited checker code."""
import argparse, hashlib, itertools, json, math
from fractions import Fraction
from pathlib import Path

def need(v, why):
    if not v: raise RuntimeError(why)

def perms(n, mults):
    return [tuple((a*x+b)%n for x in range(n)) for a in mults for b in range(n)]

def moved(w,p,alpha=None):
    # inverse convention is immaterial for this full permutation group.
    return tuple(w[i] if alpha is None else alpha[w[i]] for i in p)

def check_seed(n, diffs, mults, maxq):
    g=perms(n,mults); ident=tuple(range(n)); neg=tuple((-x)%n for x in range(n))
    allp=list(itertools.permutations(range(n)))
    aut=[p for p in allp if all(((y-x)%n in diffs)==((p[y]-p[x])%n in diffs) for x in range(n) for y in range(n) if x!=y)]
    need(set(aut)==set(g),f'full automorphism group n={n}')
    anti=[p for p in allp if all(((y-x)%n in diffs)!=((p[y]-p[x])%n in diffs) for x in range(n) for y in range(n) if x!=y)]
    invol=[p for p in anti if all(p[p[i]]==i for i in range(n))]
    need(len(anti)==len(g) and len(invol)==n,'anti group counts')
    rows=[]
    for q in range(maxq+1):
        # Form ordinary orbits directly; regularity checked by all group elements.
        regular=[]; reps=[]
        for w in itertools.product(range(q),repeat=n):
            if any(moved(w,p)==w for p in g if p!=ident): continue
            regular.append(w)
            if w==min(moved(w,p) for p in g): reps.append(w)
        expectedQ={3:(q**3-q)//3,5:(q**5-q)//5,7:(q**7-7*q**3+6*q)//21}[n]
        need(len(reps)==expectedQ and len(regular)==len(g)*expectedQ,'Q')
        for f in range(q%2,q+1,2):
            alpha=list(range(q))
            for i in range(f,q,2): alpha[i],alpha[i+1]=i+1,i
            F=sum(min(moved(moved(w,neg,alpha),p) for p in g)==w for w in reps)
            expectedF={3:f*(q-1),5:f*(q*q-1),7:f*q*(q*q-1)//3}[n]
            need(F==expectedF,'twisted F')
            # Independent anti-stabilizer incidence census, all anti maps for q<=4.
            fixed0=sum(moved(w,neg,alpha)==w for w in regular)
            need(fixed0*n==F*len(g),'reflection incidence')
            if q<=4:
                counts=[sum(moved(w,p,alpha)==w for p in anti) for w in regular]
                need(max(counts,default=0)<=1,'unique regular anti-stabilizer')
                need(counts.count(0)==len(g)*(expectedQ-expectedF),'signed regular raw count')
            rows.append(dict(q=q,f=f,Q=expectedQ,F=F,regular_words=len(regular),fixed_reflection=fixed0))
    return dict(n=n,group_order=len(g),anti_involutions=len(invol),max_q=maxq,rows=rows)

def dynamics():
    rows=[]
    for c in range(2,7):
        q=f=c; a=1; R=c
        for k in range(5):
            n=3**k
            need(a==3**((n-1)//2),'tower aut')
            need(R==a*q,'raw colors')
            row=dict(c=c,k=k,n=n,q=str(q),f=str(f),aut=str(a))
            if k<3: row['ordinary_probability']=str(Fraction(R,c**n))
            rows.append(row)
            nq=(q**3-q)//3; nf=f*(q-1); na=3*a**3
            need(Fraction(nf,nq)==Fraction(f,q)*Fraction(3,q+1),'ratio')
            nR=na*nq
            need(nR==R**3*Fraction(q*q-1,q*q),'raw recurrence')
            if c==2: need(nq==nf==2,'binary fixed point')
            else: need(nq>q and nf<=nq,'growth')
            q,f,a,R=nq,nf,na,nR
    return rows

def labelled():
    rows=[]; coeff={}
    for d,ax in [(3,3),(5,5),(7,21)]:
        a=ax
        for k in range(6):
            n=d*3**k
            need(a==ax*3**((n-d)//2),'exception aut')
            need(n not in coeff,'family orders disjoint')
            coeff[n]=Fraction(1,a)
            if n<=63: rows.append(dict(d=d,k=k,n=n,aut=str(a),labelled=str(math.factorial(n)//a)))
            a*=3**n
    for n,v in coeff.items():
        expected=Fraction(1,{3:3,5:5,7:21}[n]) if n in [3,5,7] else coeff[n//3]/3**(n//3)
        need(v==expected,'labelled EGF functional equation')
    return rows

def main():
    p=argparse.ArgumentParser(); p.add_argument('--output',required=True); args=p.parse_args()
    obj=dict(status='FINITE-EXACT, supplements all-order proofs',universe='All q-color words on C3,R5,P7 for q=0..6; every palette involution conjugacy type f=0..q with f=q mod 2; full anti-stabilizer incidence checks q<=4. No imported checker or receipt.',seeds=[check_seed(3,{1},[1],6),check_seed(5,{1,2},[1],6),check_seed(7,{1,2,4},[1,2,4],6)],dynamics=dynamics(),labelled=labelled())
    obj['source_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    Path(args.output).write_text(json.dumps(obj,indent=2)+'\n'); print('PASS: '+args.output)
if __name__=='__main__': main()
