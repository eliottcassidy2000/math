"""Standalone replay of eea15f5d968 ordered-run witnesses; no repo imports."""
from itertools import permutations
from functools import lru_cache
import json

def need(ok,msg):
    if not ok:raise RuntimeError(msg)
def arc(t,u,v):return (t[u]>>v)&1
def tournament(t):
    return all(not arc(t,i,i) and all(arc(t,i,j)+arc(t,j,i)==1 for j in range(len(t)) if i!=j) for i in range(len(t)))
def converse(t):
    return tuple(((1<<len(t))-1)^(1<<v)^r for v,r in enumerate(t))
def strong(t):
    for root in range(len(t)):
        reach=1<<root
        while True:
            more=reach
            for v in range(len(t)):
                if reach>>v&1:more|=t[v]
            if more==reach:break
            reach=more
        if reach!=(1<<len(t))-1:return False
    return True
def is_auto(t,g):
    return all(arc(t,u,v)==arc(t,g[u],g[v]) for u in range(len(t)) for v in range(len(t)))
def c3_substitute(children):
    starts=[0,len(children[0]),len(children[0])+len(children[1])]
    out=[]
    for i,t in enumerate(children):
        nxt=(i+1)%3;outbits=((1<<len(children[nxt]))-1)<<starts[nxt]
        out.extend((r<<starts[i])|outbits for r in t)
    return tuple(out)
def partner(t,p):
    out=list(converse(t))
    need(sorted(p)==list(range(len(t))),'path visits each vertex')
    for u,v in zip(p,p[1:]):
        need(arc(t,u,v),'legal source Hamilton path')
        out[u]^=1<<v;out[v]^=1<<u
    need(tournament(out),'partner is tournament')
    return tuple(out)
def forward(t):
    n=len(t);a=[[0]*n for _ in range(1<<n)]
    for v in range(n):a[1<<v][v]=1
    for mask in range(1,1<<n):
        for v in range(n):
            c=a[mask][v]
            if not c:continue
            unused=t[v]&~mask
            while unused:
                b=unused&-unused;unused-=b
                a[mask|b][b.bit_length()-1]+=c
    return sum(a[-1])
def backward(t):
    n=len(t)
    @lru_cache(None)
    def tails(v,remain):
        if not remain:return 1
        return sum(tails(w,remain^(1<<w)) for w in range(n) if remain>>w&1 and arc(t,v,w))
    return sum(tails(v,((1<<n)-1)^(1<<v)) for v in range(n))
def parity(p):return sum(p[i]>p[j] for i in range(len(p)) for j in range(i+1,len(p)))%2

def main():
    A=(10,12,9,16,7);D=(2,4,9,3);T4=(14,12,8,0)
    need(tournament(A) and strong(A) and strong(D),'strong children')
    G=[g for g in permutations(range(5)) if is_auto(A,g)]
    need(len(G)==3 and all(g[3]==3 and g[4]==4 for g in G),'full automorphism group C3 fixes outer ports')
    need(all(is_auto(converse(A),g) for g in G),'same automorphisms on forest complement')
    U=(3,0,1,2,4);V=(3,0,2,1,4)
    witnesses=[]
    for u,v in zip(U,V):
        maps=[g for g in G if g[3]==3 and g[4]==4 and g[u]==v]
        need(bool(maps),'each indexed cursor response matches on both sheets')
        witnesses.append(maps[0])
    need(not any(tuple(g[u] for u in U)==V for g in G),'no simultaneous transporter')
    need(arc(A,U[1],U[2])!=arc(A,V[1],V[2]),'one joint arc detects parity')
    source_rows={};entries=[]
    for name,child,expected in [('transitive',T4,(2483993,2491877)),('strong',D,(2899675,2906091))]:
        S=c3_substitute((A,child,child));source_rows[name]=S
        need(strong(S),'assembled source strong')
        for mid in permutations(range(3)):
            W=(3,)+mid+(4,);P=[]
            for i,v in enumerate(W):
                P.append(v)
                if i<4:P.extend((5+i,9+i))
            B=partner(S,P);h=forward(B)
            need(h==backward(B),'independent forward/backward subset recurrence')
            need(h==expected[parity(mid)],'exact supplied parity count')
            need(strong(B),'assembled partner strong')
            orbit=min(tuple(g[v] for v in W) for g in G)
            entries.append(dict(children=name,middle=mid,parity=parity(mid),joint_orbit=orbit,path=P,partner_rows=B,H=h))
    need(len({tuple(e['joint_orbit']) for e in entries})==2,'exactly two order orbits')
    print(json.dumps(dict(source_pin='eea15f5d968',A=A,D=D,T4=T4,automorphisms=G,U=U,V=V,cursor_maps=witnesses,source_rows=source_rows,entries=entries),indent=2))
if __name__=='__main__':main()
