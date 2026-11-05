#!/usr/bin/env python3
"""Standalone exact controls for the single-crossing sector theorem."""
from collections import Counter
from functools import lru_cache
from hashlib import sha256
from itertools import combinations, permutations
from pathlib import Path
import json
HERE=Path(__file__).resolve().parent

def require(q,msg):
    if not q:raise RuntimeError(msg)

def induced(rows,S):
    vs=[v for v in range(len(rows)) if S>>v&1]
    return tuple(sum(1<<j for j,y in enumerate(vs) if rows[x]>>y&1) for x in vs),vs

def converse(rows):
    F=(1<<len(rows))-1
    return tuple(F^x^(1<<i) for i,x in enumerate(rows))

@lru_cache(maxsize=20000)
def profile(rows):
    n=len(rows);N=1<<n;F=N-1
    start=[[0]*n for _ in range(N)];end=[[0]*n for _ in range(N)];h=[0]*N;h[0]=1
    inc=converse(rows)
    for S in range(1,N):
        for a in range(n):
            if not S>>a&1:continue
            Q=S^(1<<a)
            if Q==0:start[S][a]=end[S][a]=1
            else:
                start[S][a]=sum(start[Q][b] for b in range(n) if Q>>b&1 and rows[a]>>b&1)
                end[S][a]=sum(end[Q][b] for b in range(n) if Q>>b&1 and inc[a]>>b&1)
        h[S]=sum(start[S])
        require(h[S]==sum(end[S]),'forward/backward path counts disagree')
    R=[sum(start[S][a]*h[F^S] for S in range(N)) for a in range(n)]
    L=[sum(end[S][a]*h[F^S] for S in range(N)) for a in range(n)]
    return h,start,end,R,L

def components(rows,S):
    if S==0:return []
    reach=[(rows[i]&S)|(1<<i) for i in range(len(rows))]
    for k in range(len(rows)):
        if S>>k&1:
            for i in range(len(rows)):
                if reach[i]>>k&1:reach[i]|=reach[k]
    groups=[];remaining=S
    while remaining:
        i=(remaining&-remaining).bit_length()-1
        C=sum(1<<j for j in range(len(rows)) if remaining>>j&1 and reach[i]>>j&1 and reach[j]>>i&1)
        groups.append(C);remaining^=C
    snapshot=tuple(groups)
    groups.sort(key=lambda C:-sum(bool(rows[(C&-C).bit_length()-1]&D) for D in snapshot if D!=C))
    for i,C in enumerate(groups):
        for D in groups[i+1:]:
            require(all(rows[a]&D==D for a in range(len(rows)) if C>>a&1),'bad SCC order')
    return groups

def g(n):
    if n<4:raise ValueError(n)
    if n<=8:return (2,4,6,10,16)[n-4]
    k,r=divmod(n-9,3)
    return (26,46,78)[r]*5**k

def literal_h(rows):
    return sum(all(rows[a]>>b&1 for a,b in zip(p,p[1:])) for p in permutations(range(len(rows))))

def analyse(rows,v,check_remainder=True):
    n=len(rows);F=(1<<n)-1;B=rows[v];A=F^B^(1<<v)
    h,_,_,_,_=profile(rows);delta=h[F]-h[F^(1<<v)]
    ar,av=induced(rows,A);br,bv=induced(rows,B)
    ha,sa,ea,ra,la=profile(ar);hb,sb,eb,rb,lb=profile(br)
    arcs=[(b,a) for a in av for b in bv if rows[b]>>a&1]
    terms=[ra[av.index(a)]*lb[bv.index(b)]+ea[-1][av.index(a)]*sb[-1][bv.index(b)] for b,a in arcs]
    require(all(t>=2 and t%2==0 for t in terms),'nonpositive or odd single-crossing summand')
    lower=sum(terms);require(delta>=lower and (delta-lower)%2==0,'single-crossing inequality failed')
    star=not arcs or len({a for b,a in arcs})==1 or len({b for b,a in arcs})==1
    if star:require(delta==lower,'star equality failed')
    ca=components(rows,A);cb=components(rows,B);k=0;extra5=0
    strong=len(components(rows,F))==1
    if strong and n>=3:
        require(A and B,'strong root has empty side')
        k=sum(bool(ca[0]>>a&1 and cb[-1]>>b&1) for b,a in arcs)
        if n>=4 and k:
            require(delta>=k*g(n)+2*(len(arcs)-k),'extreme sector floor failed')
            for (b,a),term in zip(arcs,terms):
                if ca[0]>>a&1 and cb[-1]>>b&1:require(term>=g(n),'extreme term floor failed')
        if k==0:
            require(len(ca)>1 and len(cb)>1,'residual has a strong side')
            for b,a in arcs:
                if not ca[0]>>a&1:continue
                for d,c in arcs:
                    if not cb[-1]>>d&1:continue
                    p=(v,b,a,d,c)
                    require(len(set(p))==5 and all(rows[p[i]]>>p[(i+1)%5]&1 for i in range(5)),'residual 5-cycle failed')
                    extra5+=2*h[F^sum(1<<x for x in p)]
            if check_remainder:require(delta>=lower+extra5,'forced five-cycle charge failed')
    return dict(delta=delta,lower=lower,triangles=len(arcs),extreme_arcs=k,star=star,strong=strong,forced_five_cycle_charge=extra5)

def literal_single_crossing(rows,v):
    n=len(rows);F=(1<<n)-1;B=rows[v];A=F^B^(1<<v)
    h,*_=profile(rows);total=0;cycles=0
    others=[i for i in range(n) if i!=v]
    for length in range(2,n,2):
        for p in permutations(others,length):
            cyc=(v,)+p
            if not all(rows[cyc[i]]>>cyc[(i+1)%len(cyc)]&1 for i in range(len(cyc))):continue
            crossings=sum(bool(B>>a&1 and A>>b&1) for a,b in zip(p,p[1:]))
            if crossings!=1:continue
            mask=sum(1<<x for x in p)
            total+=2*h[A&~mask]*h[B&~mask];cycles+=1
    return total,cycles

def labelled(n):
    pp=list(combinations(range(n),2))
    for z in range(1<<len(pp)):
        rows=[0]*n
        for k,(a,b) in enumerate(pp):
            if z>>k&1:rows[a]|=1<<b
            else:rows[b]|=1<<a
        yield tuple(rows)

def check_parity(rows):
    h,_,end,R,L=profile(rows);F=(1<<len(rows))-1;count=0
    for a in range(len(rows)):
        odd=even=0
        def visit(path,S):
            nonlocal odd,even,count
            count+=1
            if len(path)%2:odd+=h[F^S]
            else:even+=h[F^S]
            nxt=rows[path[-1]]&F
            for b in range(len(rows)):
                if nxt>>b&1 and not S>>b&1:visit(path+(b,),S|(1<<b))
        visit((a,),1<<a)
        require(odd+even==R[a] and odd-even==end[-1][a],'literal parity identity failed')
    return count

def residual_family(r,s):
    n=r+s+3;A=list(range(1,r+2));B=list(range(r+2,n));rows=[0]*n
    for a in A:
        rows[a]|=1
        for b in B:rows[a]|=1<<b
    rows[0]=sum(1<<b for b in B)
    for V in (A,B):
        for i,a in enumerate(V):
            for b in V[i+1:]:rows[a]|=1<<b
    for b,a in ((B[0],A[0]),(B[-1],A[-1])):
        rows[b]|=1<<a;rows[a]^=1<<b
    return tuple(rows)

def main():
    counts=Counter();bankstats=Counter();residual_min={};sector_min={}
    for n in range(1,6):
        for rows in labelled(n):
            counts['labelled_tournaments']+=1
            counts['literal_rooted_paths']+=check_parity(rows)
            require(profile(rows)[0][-1]==literal_h(rows),'literal Hamilton count failed')
            for v in range(n):
                z=analyse(rows,v);literal,cycles=literal_single_crossing(rows,v)
                require(z['lower']==literal,'literal single-crossing object count failed')
                counts['literal_single_crossing_root_cycles']+=cycles
                counts['labelled_root_checks']+=1
    path=HERE/'inherited_tournaments_through8.txt';lines=path.read_text().splitlines()
    require(int(lines[0])==len(lines)-1==7412,'unexpected bank size')
    for line in lines[1:]:
        n,*rows=map(int,line.split());rows=tuple(rows)
        require(len(rows)==n,'bad bank row')
        if n<4 or len(components(rows,(1<<n)-1))!=1:continue
        for v in range(n):
            z=analyse(rows,v);bankstats[f'order_{n}_strong_roots']+=1
            sector='extreme' if z['extreme_arcs'] else 'residual'
            bankstats[f'order_{n}_{sector}']+=1
            dst=sector_min if sector=='extreme' else residual_min
            if str(n) not in dst or z['delta']<dst[str(n)]['delta']:dst[str(n)]={**z,'rows':rows,'root':v}
    hostile=(48,5,9,3,14,30);hz=analyse(hostile,0)
    require(hz['delta']==hz['lower']==30 and hz['extreme_arcs']==3,'H33 hostile changed')
    require(literal_h(hostile)==33,'literal H33 failed')
    hdel=induced(hostile,((1<<6)-1)^1)[0];require(literal_h(hdel)==3,'literal Hminus33 failed')
    edits=[]
    for b in (4,5):
        for a in (1,2,3):
            rr=list(hostile);rr[a]^=1<<b;rr[b]^=1<<a;rr=tuple(rr)
            zz=analyse(rr,0);require(zz['delta']==(32 if b==4 else 34),'smoothing hostile failed')
            edits.append({'arc':[b,a],**zz})
    residuals=[]
    for r,s in ((1,1),(2,2),(3,3)):
        rr=residual_family(r,s);z=analyse(rr,0)
        require(z['lower']==2**r+2**s and z['extreme_arcs']==0,'residual family incorrect')
        if r<=2:require(literal_h(rr)-literal_h(induced(rr,((1<<len(rr))-1)^1)[0])==z['delta'],'literal residual failed')
        residuals.append({'r':r,'s':s,'rows':rr,**z})
    return dict(status='PASS',labelled_controls=dict(counts),bank_sha256=sha256(path.read_bytes()).hexdigest(),bank_counts=dict(bankstats),extreme_sector_minima=sector_min,residual_sector_minima=residual_min,H33={**hz,'single_arc_edits':edits},residual_examples=residuals)

if __name__=='__main__':print(json.dumps(main(),indent=2,sort_keys=True))
