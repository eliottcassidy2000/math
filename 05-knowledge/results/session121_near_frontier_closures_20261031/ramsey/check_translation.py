#!/usr/bin/env python3
"""Finite exact controls only; not a finite-table proof of stabilization."""
from fractions import Fraction
from functools import lru_cache
from itertools import product, combinations
import json

def need(ok,msg):
    if not ok: raise RuntimeError(msg)

@lru_cache(None)
def w(a):
    if a<2:return Fraction(a)
    h=a//2
    return Fraction(5,2)*w(h) if a%2==0 else Fraction(3,2)*w(h)+w(h+1)

@lru_cache(None)
def patterns(k,step,maximal):
    out=[]
    for t in product(range(1,k),repeat=5):
        if any(t[i]+t[(i+step)%5]>k for i in range(5)):continue
        if maximal and any(t[i]+t[(i+step)%5]<k and t[i]+t[(i-step)%5]<k for i in range(5)):continue
        out.append(t)
    return out

def table(amax,bmax,maximal):
    f={};plans={};count=0
    for a in range(1,amax+1):
        for b in range(1,bmax+1):
            best=1;plan=('v',)
            for x in range(1,a):
                z=f[x,b]+f[a-x,b]
                if z>best:best=z;plan=('j',x)
            for y in range(1,b):
                z=f[a,y]+f[a,b-y]
                if z>best:best=z;plan=('u',y)
            for aa in patterns(a,1,maximal):
                for bb in patterns(b,2,maximal):
                    z=sum(f[x,y] for x,y in zip(aa,bb));count+=1
                    if z>best:best=z;plan=('c',aa,bb)
            f[a,b]=best;plans[a,b]=plan
    return f,plans,count

def row(a,b):
    return (b,5*b//2,4*b-int(b%4!=0),25*(b//4)+(0,4,10,16)[b%4])[a-1]

def literal_checks(f,plans):
    def substitute(children,edges):
        starts=[];s=0
        for c in children:starts.append(s);s+=len(c)
        rows=[]
        for i,c in enumerate(children):
            outside=sum(((1<<len(d))-1)<<starts[j] for j,d in enumerate(children) if edges(i,j))
            rows.extend((r<<starts[i])|outside for r in c)
        return tuple(rows)
    @lru_cache(None)
    def build(a,b):
        p=plans[a,b]
        if p[0]=='v':return (0,)
        if p[0]=='j':return substitute([build(p[1],b),build(a-p[1],b)],lambda i,j:i!=j)
        if p[0]=='u':return substitute([build(a,p[1]),build(a,b-p[1])],lambda i,j:False)
        return substitute([build(x,y) for x,y in zip(p[1],p[2])],lambda i,j:(i-j)%5 in (1,4))
    def clique(rows):
        @lru_cache(None)
        def rec(mask):
            if not mask:return 0
            v=(mask&-mask).bit_length()-1
            return max(rec(mask&~(1<<v)),1+rec(mask&rows[v]))
        return rec((1<<len(rows))-1)
    checked=[]
    for a in range(1,5):
        for b in range(1,7):
            rs=build(a,b);n=len(rs);full=(1<<n)-1
            need(n==f[a,b],'literal order')
            need(all(not(r>>i&1) for i,r in enumerate(rs)),'loops')
            need(all(bool(rs[i]>>j&1)==bool(rs[j]>>i&1) for i in range(n) for j in range(n)),'symmetry')
            om=clique(rs);al=clique(tuple(full^r^(1<<i) for i,r in enumerate(rs)))
            need(om<=a and al<=b,'literal clique and independence')
            checked.append([a,b,n,om,al])
    return checked

def closure(gens,limit,base=(0,0)):
    seen={base};todo=[base]
    while todo:
        x=todo.pop()
        for v in gens:
            y=tuple(a+b for a,b in zip(x,v))
            if max(y)<=limit and y not in seen:seen.add(y);todo.append(y)
    return seen

def star_check():
    # First form each linear set literally inside the box, then close its
    # actual points under addition. Compare the activated-component formula.
    families=[ [((1,1),[(0,1)])],
        [((1,2),[(1,0),(0,2)]),((2,1),[(0,1)])],
        [((2,3),[]),((3,2),[(2,0)]),((1,4),[(0,1)])],
        [((0,2),[(1,0)]),((3,0),[(0,1)])] ]
    receipts=[]
    for parts in families:
        lim=18
        points=set().union(*(closure(v,lim,u) for u,v in parts))
        direct=closure(points,lim)
        compiled={(0,0)}
        for m in range(1,1<<len(parts)):
            ids=[i for i in range(len(parts)) if m>>i&1]
            base=tuple(sum(parts[i][0][j] for i in ids) for j in range(2))
            if max(base)>lim:continue
            gens=[parts[i][0] for i in ids]+[v for i in ids for v in parts[i][1]]
            compiled.update(closure(gens,lim,base))
        need(direct==compiled,'activated semilinear star mismatch')
        receipts.append(len(direct))
    need((0,1) not in closure(closure([(0,1)],18,(1,1)),18),'activation guard missing')
    return receipts

def main():
    f,plans,cnt=table(6,16,True)
    small,_,cnt2=table(4,6,False)
    need(all(f[k]==v for k,v in small.items()),'full vs maximal-pattern DP')
    for a in range(1,5):
        for b in range(1,17):need(f[a,b]==row(a,b),'inherited exact row')
    inequalities=0
    for (a,b),n in f.items():
        need(n<=min(w(a)*b,w(b)*a),'capacity upper bound')
        L=1<<(a-1).bit_length();C=w(a)*L
        need(C.denominator==1,'sharp integer')
        if b%L==0:need(n==w(a)*b,'sharp progression')
        if b+L<=16:
            need(f[a,b+L]>=n+C,'union translation inequality');inequalities+=1
        for c in range(1,17-b):need(f[a,b+c]>=n+f[a,c],'union superadditivity')
    periods=[]
    for a in range(1,5):
        L=1<<(a-1).bit_length()
        residues=[w(a)*b-row(a,b) for b in range(1,L+1)]
        d=next(d for d in range(1,L+1) if L%d==0 and all(residues[i]==residues[(i+d)%L] for i in range(L)))
        need(d%w(a).denominator==0,'period denominator guard');periods.append(d)
    hostiles=0
    for M in (3,11,101):
        def hostile(b):return 2*b-int(b%2==1 and b<M)
        for x in range(1,120):
            for y in range(1,120):need(hostile(x+y)>=hostile(x)+hostile(y),'plateau hostile superadditivity');hostiles+=1
        need(hostile(2)==4 and hostile(M)==2*M,'delayed improvement')
    out={'scope':'bounded controls, not a proof of an unseen eventual tail',
         'capacity_rectangle':[6,16],'maximal_pattern_candidates':cnt,
         'independent_full_pattern_rectangle':[4,6],'full_pattern_candidates':cnt2,
         'capacities':[[f[a,b] for b in range(1,17)] for a in range(1,7)],
         'translation_inequalities':inequalities,'literal_graphs':literal_checks(f,plans),
         'known_exact_periods_rows1to4':periods,'star_box18_counts':star_check(),
         'delayed_plateau_superadditivity_checks':hostiles}
    print(json.dumps(out,sort_keys=True,indent=2))
if __name__=='__main__':main()
