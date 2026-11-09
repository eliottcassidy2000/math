#!/usr/bin/env python3
"""Three-block escape verifier. Binary compiler inherited from session228 verify.py.
Standalone, standard library, explicit checks active under -O.
"""
from itertools import combinations, permutations
from math import comb, factorial
from fractions import Fraction
from hashlib import sha256
import json

X=4095
N=2+(X-1).bit_length()-1

def need(ok,msg):
    if not ok: raise RuntimeError(msg)
def val(p,t):return sum(c*t**i for i,c in enumerate(p))
def mul(a,b):
    out=[0]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):out[i+j]+=x*y
    return tuple(out)
def shift(p,s):
    return tuple(sum(p[j]*comb(j,i)*s**(j-i) for j in range(i,len(p))) for i in range(len(p)))
def sub(a,b):return tuple(x-y for x,y in zip(a,b))
def cyc(a,b):
    u=mul(sub(shift(a,1),a),sub(shift(b,1),b))
    v=mul(sub(a,shift(a,-1)),sub(b,shift(b,-1)))
    return (0,)+tuple(x+y+z for x,y,z in zip(mul(a,b),u,v))
def collar(p):return cyc(p,(0,1))
def tournaments(n):
    pairs=list(combinations(range(n),2))
    for bits in range(1<<len(pairs)):
        r=[0]*n
        for k,(i,j) in enumerate(pairs):
            a,b=(i,j) if bits>>k&1 else (j,i);r[a]|=1<<b
        yield tuple(r)
def join(a,b):
    n=len(a);bm=((1<<len(b))-1)<<n
    return tuple(x|bm for x in a)+tuple(x<<n for x in b)
def cycle(a,b):
    # A -> B -> z -> A.
    n,m=len(a),len(b);z=n+m;bm=((1<<m)-1)<<n
    return tuple(x|bm for x in a)+tuple((x<<n)|(1<<z) for x in b)+((1<<n)-1,)
def hsub(r):
    n=len(r);size=1<<n;dp=[[0]*n for _ in range(size)]
    for i in range(n):dp[1<<i][i]=1
    for mask in range(1,size):
        for i,z in enumerate(dp[mask]):
            if not z:continue
            todo=r[i]&~mask
            while todo:
                bit=todo&-todo;todo-=bit
                dp[mask|bit][bit.bit_length()-1]+=z
    h=list(map(sum,dp));h[0]=1
    return h
def qvalues(r,k):
    h=hsub(r);size=len(h);q=[0]*size;q[0]=1;out=[0]
    for _ in range(k):
        p=[0]*size
        for mask in range(size):
            s=mask
            while True:
                p[mask]+=h[s]*q[mask^s]
                if not s:break
                s=(s-1)&mask
        q=p;out.append(q[-1])
    return out
def permutation_poly(r):
    n=len(r);F=[0]*(n+1)
    for p in permutations(range(n)):
        b=sum(not(r[x]>>y&1) for x,y in zip(p,p[1:]))
        for k in range(b+1,n+1):F[k]+=comb(n-1-b,k-1-b)
    out=[0]*(n+1);fall=(1,)
    for d in range(1,n+1):
        need(F[d]%factorial(d)==0,'unordered cover divisibility')
        fall=mul(fall,(-(d-1),1))
        for i,c in enumerate(fall):out[i]+=c*(F[d]//factorial(d))
    return tuple(out)
def rank(A):
    A=[[Fraction(x) for x in r] for r in A];k=0
    for j in range(len(A[0])):
        at=next((i for i in range(k,len(A)) if A[i][j]),None)
        if at is None:continue
        A[k],A[at]=A[at],A[k];c=A[k][j];A[k]=[v/c for v in A[k]]
        for i in range(len(A)):
            if i!=k:
                c=A[i][j];A[i]=[a-c*b for a,b in zip(A[i],A[k])]
        k+=1
    return k

def compile_profiles(seeds,limit):
    # Values are (minimum cyclic depth, canonical expression).
    bank=[{} for _ in range(limit+1)]
    def put(n,p,d,e):
        if val(p,1)>X:return
        if p not in bank[n] or (d,e)<bank[n][p]:bank[n][p]=(d,e)
    for name,rows,p in seeds:
        if len(rows)<=limit:put(len(rows),p,0,name)
    for n in range(2,limit+1):
        for i in range(1,n//2+1):
            for a,(da,ea) in sorted(bank[i].items()):
                for b,(db,eb) in sorted(bank[n-i].items()):put(n,mul(a,b),max(da,db),'J('+ea+','+eb+')')
        for i in range(1,(n-1)//2+1):
            for a,(da,ea) in sorted(bank[i].items()):
                for b,(db,eb) in sorted(bank[n-1-i].items()):put(n,cyc(a,b),1+max(da,db),'C('+ea+','+eb+')')
    return bank

def graph_expr(expr,seedrows):
    if expr[0]=='s':return seedrows[expr]
    op=expr[0];body=expr[2:-1];depth=0
    for i,c in enumerate(body):
        if c=='(':depth+=1
        if c==')':depth-=1
        if c==',' and depth==0:break
    a=graph_expr(body[:i],seedrows);b=graph_expr(body[i+1:],seedrows)
    return join(a,b) if op=='J' else cycle(a,b)

def cover(p):
    vals=[val(p,k) for k in range(len(p))];out=[]
    while vals:
        out.append(vals[0]); vals=[b-a for a,b in zip(vals,vals[1:])]
    return out+[0]
def cscalar(A,B,C):
    return sum(3*A[r]*B[r]*C[r]+A[r+1]*B[r]*C[r]+A[r]*B[r+1]*C[r]+A[r]*B[r]*C[r+1]+A[r+1]*B[r+1]*C[r]+A[r+1]*B[r]*C[r+1]+A[r]*B[r+1]*C[r+1] for r in range(1,min(len(A),len(B),len(C))-1))
def tt(n):return tuple(((1<<n)-1)^((1<<(i+1))-1) for i in range(n))
def subst(R,blocks):
    offsets=[];n=0
    for b in blocks:offsets.append(n);n+=len(b)
    masks=[((1<<len(b))-1)<<o for b,o in zip(blocks,offsets)]
    out=[]
    for i,b in enumerate(blocks):
        cross=sum(masks[j] for j in range(len(R)) if R[i]>>j&1)
        out.extend((r<<offsets[i])|cross for r in b)
    return tuple(out)
def cycle_certificate(R):
    return next((p for p in permutations(range(len(R))) if p[0]==0 and all(R[p[i]]>>p[(i+1)%len(R)]&1 for i in range(len(R)))),None)
def balance(q,n):
    a,b=divmod(n,q);return (q-b)*2**a+b*2**(a+1)-q

def main():
    from itertools import product,combinations_with_replacement
    C3=(2,4,1)
    distinct={};seed_checks=0
    for n in range(1,6):
        for rows in tournaments(n):
            p=permutation_poly(rows)
            need(qvalues(rows,2)==[val(p,j) for j in range(3)],'seed independent Q1/Q2')
            distinct.setdefault((n,p),rows);seed_checks+=1
    seeds=[('s'+str(i),r,p) for i,((n,p),r) in enumerate(sorted(distinct.items()))]
    seedrows={name:r for name,r,p in seeds}
    bank=compile_profiles(seeds,N)
    generators={val(p,1) for row in bank for p in row}
    reachable={1}
    for h in range(3,X+1,2):
        if h in generators or any(h%d==0 and h//d in reachable for d in generators if 1<d<h):reachable.add(h)
    need(len(reachable)==577,'inherited complete scalar count')
    need(19 not in reachable and 137 not in reachable,'old grammar hostiles')
    digest=sha256(json.dumps(sorted(reachable),separators=(',',':')).encode()).hexdigest()
    profiles=[(n,p,cover(p),e) for n,row in enumerate(bank) for p,(d,e) in row.items() if n<=5]
    found={};candidate_checks=0
    for A,B,C in combinations_with_replacement(profiles,3):
        h=cscalar(A[2],B[2],C[2]);candidate_checks+=1
        if h<=X and h not in reachable:found.setdefault(h,(A,B,C))
    witnesses=[]
    for h,parts in sorted(found.items()):
        blocks=[graph_expr(a[3],seedrows) for a in parts]
        rows=subst(C3,blocks)
        need(hsub(rows)[-1]==h,'literal new scalar')
        witnesses.append({'H':h,'order':len(rows),'children':[a[3] for a in parts],'rows':rows})
    need(len(found)==17 and 137 in found,'escape finite image')
    exact_controls=0;literal_controls=0
    for ns in product(range(1,7),repeat=3):
        Fs=[cover((0,)*n+(1,)) for n in ns]
        h=cscalar(*Fs);aa=[2**n-2 for n in ns]
        lower=3+sum(aa)+sum(aa[i]*aa[j] for i in range(3) for j in range(i+1,3))+3*aa[0]*aa[1]*aa[2]
        need(h>=lower>=2**(sum(ns)-2)+1,'sharp ternary size bound')
        if sum(n>1 for n in ns)<=1:need(h==2**(sum(ns)-2)+1,'sharp witness')
        exact_controls+=1
        if max(ns)<=3:
            need(hsub(subst(C3,[tt(n) for n in ns]))[-1]==h,'literal ternary control');literal_controls+=1
    formula_checks=0
    for _,rows,p in seeds:
        h=4*val(p,3)+5*val(p,2)-11*val(p,1)
        need(hsub(subst(C3,[tt(2),tt(2),rows]))[-1]==h,'three-moment escape formula');formula_checks+=1
    general_checks=0;template_checks=0
    for q in (3,4):
        for R in tournaments(q):
            cycle=cycle_certificate(R)
            if cycle is None:continue
            template_checks+=1
            for ns in product((1,2),repeat=q):
                h=hsub(subst(R,[tt(n) for n in ns]))[-1]
                need(h>=sum(2**n-1 for n in ns)>=balance(q,sum(ns)),'general template bound')
                general_checks+=1
    for q in range(3,7):
        for ns in product(range(1,5),repeat=q):
            need(sum(2**n-1 for n in ns)>=balance(q,sum(ns)),'balancing lower bound')
    hostile=join(C3,tt(12));need(hsub(hostile)[-1]==3 and len(hostile)>N,'outer-join order hostile')
    Erows=subst(C3,[tt(2),tt(2),tt(3)])
    print(json.dumps({'status':'PASS','scope':'Exact old grammar through4095; single ternary nodes with old children order<=5, not full expanded grammar','cutoff':X,'old_bound':N,'seed_presentations':seed_checks,'seed_profiles':len(seeds),'old_scalar_count':len(reachable),'old_scalar_sha256':digest,'ternary_candidate_triples':candidate_checks,'new_count':len(found),'new_values':sorted(found),'new_witnesses':witnesses,'transitive_kernel_controls':exact_controls,'transitive_literal_controls':literal_controls,'escape_formula_literal_controls':formula_checks,'general_template_count':template_checks,'general_template_literal_controls':general_checks,'escape_137_rows':Erows,'escape_family_m1_through10':[4*3**m+5*2**m-11 for m in range(1,11)],'nineteen_still_excluded':'Every genuinely new ternary node has H>=45; old grammar excludes19.'},indent=2,sort_keys=True))
if __name__=='__main__':main()
