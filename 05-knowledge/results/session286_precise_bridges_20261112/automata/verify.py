#!/usr/bin/env python3
"""Exact physical full-central fiber, binary Ising and semigroup checks.
Standard library only. No project imports; explicit gates survive python -O.
"""
from collections import defaultdict, deque
from functools import lru_cache
from itertools import permutations, product
from math import comb
from pathlib import Path
import json

BASE=Path(__file__).resolve().parent
D=(18,8,3,5,14)
KEYS=('H','h','Q','U','V','d','e')

def require(ok,msg):
    if not ok: raise RuntimeError(msg)
def attach(c):
    return (c<<1,)+tuple((row<<1)|(1-((c>>i)&1)) for i,row in enumerate(D))
def reverse(T):
    return tuple(sum(1<<u for u in range(len(T)) if T[u]>>v&1) for v in range(len(T)))
def hp(T):
    n=len(T); incoming=reverse(T)
    @lru_cache(None)
    def end(S,v):
        R=S^(1<<v)
        if not R:return 1
        bits=R&incoming[v];ans=0
        while bits:
            b=bits&-bits;bits-=b;ans+=end(R,b.bit_length()-1)
        return ans
    H=[1]+[sum(end(S,v) for v in range(n) if S>>v&1) for S in range(1,1<<n)]
    return H,end

def packet(T):
    H,end=hp(T);RH,rend=hp(reverse(T));F=(1<<len(T))-1;DD=F^1;vals=[]
    for HH,ee in ((RH,rend),(H,end)):
        mass=signed=0
        for S in range(1,F+1,2):
            x=ee(S,0)*HH[F^S];mass+=x;signed+=x*(1 if S.bit_count()%2 else -1)
        vals.extend((mass,signed))
    U,d,V,e=vals
    return dict(H=H[F],h=H[DD],Q=sum(H[S]*H[DD^S] for S in range(F+1) if not S&1),U=U,V=V,d=d,e=e)
def perm_packet(T):
    n=len(T);F=(1<<n)-1;paths={0:[()]}
    for S in range(1,F+1):
        paths[S]=[p for p in permutations([i for i in range(n) if S>>i&1]) if all(T[a]>>b&1 for a,b in zip(p,p[1:]))]
    U=V=d=e=0
    for S in range(1,F+1,2):
        sign=1 if S.bit_count()%2 else -1;outside=len(paths[F^S])
        a=sum(p[0]==0 for p in paths[S])*outside;b=sum(p[-1]==0 for p in paths[S])*outside
        U+=a;V+=b;d+=sign*a;e+=sign*b
    return dict(H=len(paths[F]),h=len(paths[F^1]),Q=sum(len(paths[S])*len(paths[(F^1)^S]) for S in range(F+1) if not S&1),U=U,V=V,d=d,e=e)
def multiply(P,R):
    h=P['h'];r=R['h'];hh=h*r
    return dict(H=hh+h*(R['H']-r)+(P['H']-h)*r+(P['U']-h)*(R['V']-r)+(P['d']-h)*(R['e']-r),h=hh,Q=P['Q']*R['Q'],U=hh+h*(R['U']-r)+(P['U']-h)*R['Q'],V=hh+P['Q']*(R['V']-r)+(P['V']-h)*r,d=h*R['d'],e=P['e']*r)
def word_packet(w,Ps):
    require(bool(w),'nonempty word');out=Ps[w[0]]
    for c in w[1:]:out=multiply(out,Ps[c])
    return out

def fuse(A,B):
    n=len(A);m=len(B);T=[0]*(n+m-1)
    for C,ids in ((A,list(range(n))),(B,[0]+list(range(n,n+m-1)))):
        for i,u in enumerate(ids):
            for j,v in enumerate(ids):
                if C[i]>>j&1:T[u]|=1<<v
    for i in range(1,n):
        for j in range(n,n+m-1):T[i]|=1<<j
    return tuple(T)
def literal_word(w):
    T=attach(w[0])
    for c in w[1:]:T=fuse(T,attach(c))
    return T

def accepted(w):
    state=0
    for c in w:
        if state==0:
            if c==9:pass
            elif c==3:state=1
            else:return False
        else:
            if c==18:pass
            elif c==24:state=0
            else:return False
    return state==0

def spin_word(bits):
    seq=(0,)+tuple(bits)+(0,);table={(0,0):9,(0,1):3,(1,1):18,(1,0):24}
    return tuple(table[a,b] for a,b in zip(seq,seq[1:]))
def energy(bits):return sum(bits)-9*sum(a*b for a,b in zip(bits,bits[1:]))
def semigroup(t):return t>=0 and t>=17*(t%8)
def spectrum(k):return {0}|{m-8*r for m in range(1,k//2+1) for r in range(k-2*m+1)}
def inverse(k,j):
    if j==0:return (9,)*k
    t=j-(17-8*k)
    if not (0<=t<=8*k-17+k//2 and semigroup(t)):return None
    u=t%8;a=(t-17*u)//8
    u+=8*(a//17);a%=17;r=k-2-2*u-a
    require(r>=0,'resource bound')
    return (3,)+(18,)*r+(24,)+(3,24)*u+(9,)*a

def baseline_graph(ai,Ps):
    A=Ps[ai];start=('start',);seen={start};todo=deque([start]);edges={}
    while todo:
        st=todo.popleft();out=[]
        for c,R in enumerate(Ps):
            if st==start:
                x=y=0;prev=None
                if R['e']!=A['e']:continue
            else:x,y,prev=st
            num=x+R['V']-A['V'];ny=8*y-R['U']+A['U']
            if num%8 or abs(num//8)>=13 or abs(ny)>=13:continue
            ns=(num//8,ny,c);wt=13*(R['H']-A['H'])
            if prev is not None:wt+=-(A['U']-13)*x-(R['V']-13)*y+(Ps[prev]['d']-13)*(R['e']-13)-(A['d']-13)*(A['e']-13)
            out.append((ns,c,wt))
            if ns not in seen:seen.add(ns);todo.append(ns)
        edges[st]=out
    final={st for st in seen if st!=start and st[:2]==(0,0) and Ps[st[2]]['d']==A['d']}
    rev=defaultdict(set)
    for st,ee in edges.items():
        for ns,c,wt in ee:rev[ns].add(st)
    good=set(final);todo=deque(final)
    while todo:
        for st in rev[todo.popleft()]:
            if st not in good:good.add(st);todo.append(st)
    slim={(st[:2] if st!=start else(0,0),c,ns[:2],wt) for st in good for ns,c,wt in edges[st] if ns in good}
    return len(seen),len(good),slim

def main():
    Ts=[attach(c) for c in range(32)];Ps=[packet(T) for T in Ts]
    for T,P in zip(Ts,Ps):require(P==perm_packet(T),'independent six-vertex packet')
    require(all((P['h'],P['Q'])==(13,104) for P in Ps),'core pair')
    A=Ps[9];start=('start',);seen={start};todo=deque([start]);edges={}
    while todo:
        st=todo.popleft();out=[]
        for c,R in enumerate(Ps):
            if st==start:
                x=y=0;prev=None
                if R['e']!=5:continue
            else:x,y,prev=st
            num=x+R['V']-71;ny=8*y-R['U']+46
            if num%8:continue
            nx=num//8
            if abs(nx)>=13 or abs(ny)>=13:continue
            ns=(nx,ny,c);wt=13*(R['H']-37)
            if prev is not None:wt+=-33*x-(R['V']-13)*y+(Ps[prev]['d']-13)*(R['e']-13)-40
            out.append((ns,c,wt))
            if ns not in seen:seen.add(ns);todo.append(ns)
        edges[st]=out
    final={st for st in seen if st!=start and st[:2]==(0,0) and Ps[st[2]]['d']==8}
    rev=defaultdict(set)
    for s,ee in edges.items():
        for t,c,wt in ee:rev[t].add(s)
    good=set(final);todo=deque(final)
    while todo:
        for s in rev[todo.popleft()]:
            if s not in good:good.add(s);todo.append(s)
    require(len(seen)==7 and len(good)==5,'exact carry graph')
    slim={(s[0:2] if s!=start else(0,0),c,t[:2],wt) for s,ee in edges.items() if s in good for t,c,wt in ee if t in good}
    require(slim=={((0,0),9,(0,0),0),((0,0),3,(-1,1),-52),((-1,1),18,(-1,1),-16),((-1,1),24,(0,0),54)},'trim transitions')
    baseline_results={}
    exceptions={9:(3,18,24,-1,1,-52,-16,54),13:(7,22,28,-1,1,-52,16,86),18:(24,9,3,1,-1,78,16,-44),22:(28,13,7,1,-1,78,-16,-76)}
    for ai in range(32):
        reachable,trimmed,trans=baseline_graph(ai,Ps)
        if ai in exceptions:
            bi,ei,ci,x,y,wb,we,wc=exceptions[ai]
            require(trans=={((0,0),ai,(0,0),0),((0,0),bi,(x,y),wb),((x,y),ei,(x,y),we),((x,y),ci,(0,0),wc)},'exceptional baseline exact grammar')
        else:require(trans=={((0,0),ai,(0,0),0)},'rigid baseline')
        baseline_results[str(ai)]=dict(reachable=reachable,trimmed=trimmed,transitions=sorted(trans))
    word_controls=0
    for k in (1,2,3):
        base=word_packet((9,)*k,Ps)
        for w in product(range(32),repeat=k):
            P=word_packet(w,Ps);same=all(P[a]==base[a] for a in KEYS if a!='H')
            require(same==accepted(w),'complete all32 word control');word_controls+=1
    spin_controls=0
    for k in range(1,15):
        base=word_packet((9,)*k,Ps);counts=defaultdict(int)
        for bits in product((0,1),repeat=k-1):
            w=spin_word(bits);j=energy(bits);P=word_packet(w,Ps)
            require(all(P[a]==base[a] for a in KEYS if a!='H'),'spin packet')
            require(P['H']==base['H'] if k==1 else P['H']-base['H']==2*13**(k-2)*j,'spin energy')
            counts[j]+=1;spin_controls+=1
        require(set(counts)==spectrum(k),'spin spectrum')
        for j,n in counts.items():
            mcount=int(j==0)+sum(comb(k-r-m,m)*comb(r+m-1,m-1) for m in range(1,k//2+1) for r in range(k-2*m+1) if m-8*r==j)
            require(n==mcount,'multiplicity formula')
        if k>=3:
            size=2**(k-1);sj=sum(j*n for j,n in counts.items());sj2=sum(j*j*n for j,n in counts.items())
            require(4*sj==(16-7*k)*size,'mean')
            require(16*(sj2*size-sj*sj)==(337*k-832)*size*size,'variance')
    gaps=[t for t in range(112) if not semigroup(t)]
    require(len(gaps)==56 and gaps[-1]==111,'Apery gaps')
    inverse_controls=0
    for k in range(2,151):
        vals=spectrum(k);low=17-8*k;upper=k//2
        predicted={0}|{low+t for t in range(upper-low+1) if semigroup(t)}
        require(vals==predicted,'truncated semigroup iff')
        if k>=17:
            require({j-low for j in range(low,upper+1) if j not in vals}==set(gaps),'persistent56holes')
            require(len(vals)==8*k+k//2-72,'linear spectrum size')
        for j in range(low-1,upper+2):
            w=inverse(k,j);require((w is not None)==(j in vals),'exact inverse accept')
            if w is not None:
                require(len(w)==k and accepted(w),'inverse physical word')
                # Compute the action directly from excursions, not large matrix products.
                require(sum(c==3 for c in w)-8*sum(c==18 for c in w)==j,'inverse action');inverse_controls+=1
    flow_controls=0
    for k in range(1,13):
        flow_values=set()
        for aa in range(k+1):
            for bb in range(k+1-aa):
                for cc in range(k+1-aa-bb):
                    ee=k-aa-bb-cc
                    # Four edge types00,01,10,11, pinned start/end0.
                    if bb!=cc:continue
                    if bb==0 and ee>0:continue  # disconnected1-loop
                    flow_values.add(bb-8*ee);flow_controls+=1
        require(flow_values==spectrum(k),'connected integer flow iff')
    require(-8 not in spectrum(1),'disconnected circulation hostile')
    literal=[]
    for w in ((9,9),(3,24),(3,18,24),(9,9,9)):
        T=literal_word(w);H=hp(T)[0][-1];P=word_packet(w,Ps)
        require(H==P['H'],'literal full HP count')
        if len(w)==2:require(packet(T)==P,'literal elevenvertex packet')
        literal.append(dict(word=list(w),vertices=len(T),H=H))
    for j in range(-300,301):
        if j>0:k=2*j
        elif j<0:k=2*j+17*((1-j+7)//8)
        else:k=1
        require(j in spectrum(k) and (k==1 or j not in spectrum(k-1)),'minimum length')
    for k in range(2,101):
        pos={0}|{17*m+8*r for m in range(1,k//2+1) for r in range(k-2*m+1)}
        require(pos=={0}|{17+t for t in range(17*(k//2)+8*(k%2)-16) if semigroup(t)},'positive sector semigroup')
    (BASE/'baseline_graphs.json').write_text(json.dumps(baseline_results,sort_keys=True,indent=2)+'\n')
    table={str(i):P for i,P in enumerate(Ps)}
    (BASE/'attachment_table.json').write_text(json.dumps(table,sort_keys=True,indent=2)+'\n')
    result=dict(status='PASS',scope='Complete fixed-core central fiber; no absolute odd target coverage',core=list(D),reachable_states=len(seen),trimmed_states=len(good),all_baseline_classification={'rigid':28,'two_state':[9,13,18,22]},all32_word_controls=word_controls,spin_word_controls=spin_controls,integer_flow_controls=flow_controls,semigroup_max_k=150,inverse_controls=inverse_controls,Apery_gaps=gaps,literal_counts=literal,normalized_support_at_k17=sorted(spectrum(17)),gf='(1-z*t^-8)/((1-z)*(1-z*t^-8)-z^2*t)')
    print(json.dumps(result,sort_keys=True,indent=2))
if __name__=='__main__':main()
