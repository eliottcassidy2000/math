"""Polynomial-in-order path loads for the full one-apex/one-apex F1 family."""
from collections import Counter
from functools import lru_cache
from itertools import product, permutations
import argparse
import json
from selector import profile, check_rows
from small_arm import g

def core(e):
    v=e-1
    n=v.bit_length()
    R=[((1<<i)-1)|((1<<n) if not v>>i&1 else 0) for i in range(n)]+[v]
    return R

def features(e,kind):
    """kind P (all paths), E (paths ending at root), or C (root-start covers).
    'last' tracks latest assigned vertex in each of two ascending slots.
    Slot strings are reversed when read as actual paths."""
    v=e-1;n=v.bit_length()
    def advance(i,last,b):
        prev=last[b]
        state=(i,last[1]) if b==0 else (last[0],i)
        es=[('a',i,prev)] if prev>=0 else []
        if kind=='E' and b==0:return None
        if prev<0:
            if kind in ('P','E'):
                if b==0:es.append(('end',i))
                else:
                    if v>>i&1:return None
                    es.append(('a',i,n))
            else:es.append(('pe' if b==0 else 'se',i))
        return state,es
    def finish(last):
        a,b=last
        if kind in ('P','E'):
            if a>=0 and not v>>a&1:return None
            return [('a',n,a) if a>=0 else ('end',n),
                    ('start',b if b>=0 else n)]
        if b>=0 and not v>>b&1:return None
        return [('a',n,b) if b>=0 else ('se',n)]
    @lru_cache(None)
    def back(i,last):
        if i==n:return int(finish(last) is not None)
        return sum(back(i+1,r[0]) for b in (0,1)
                   if (r:=advance(i,last,b)) is not None)
    S={(-1,-1):1};answer=Counter()
    for i in range(n):
        nxt=Counter()
        for last,p in S.items():
            for b in (0,1):
                t=advance(i,last,b)
                if t is None:continue
                state,ev=t
                q=back(i+1,state)
                if not q:continue
                nxt[state]+=p
                for e0 in ev:answer[e0]+=p*q
        S=nxt
    for last,p in S.items():
        ev=finish(last)
        if ev is not None:
            for f in ev:answer[f]+=p
    return back(0,(-1,-1)),answer

def analyze(b,e):
    assert b>=1 and e>=1
    L,R=core(b),core(e)
    return _combine(b,e,L,R,{t:features(b,t) for t in ('P','E','C')})

def explicit_features(rows):
    """Small canonical core, last vertex marked (any rooted tournament)."""
    check_rows(rows); n=len(rows);root=n-1
    sums={t:(0,Counter()) for t in ('P','E','C')}
    def add(t,events):
        a,c=sums[t];c.update(events);sums[t]=a+1,c
    for p in permutations(range(n)):
        edges=list(zip(p,p[1:]))
        bad={e for e in edges if not rows[e[0]]>>e[1]&1}
        if not bad:
            es=[('a',*e) for e in edges]+[('start',p[0]),('end',p[-1])]
            add('P',es)
            if p[-1]==root:add('E',es)
        j=p.index(root);cut=p[j-1:j+1]
        if not bad or j and bad=={cut}:
            es=[('a',*e) for e in edges if not j or e!=cut]+[('se',p[-1])]
            if j:es.append(('pe',p[j-1]))
            add('C',es)
    return sums

def analyze_bank(rows,e):
    D=explicit_features(rows)
    return _combine(D['C'][0],e,rows,core(e),D)

def _combine(b,e,L,R,D):
    l,r=len(L),len(R);n=l+r;x=l-1;y=l+r-1
    # Left module is the converse of its canonical core; right beats left except x→y.
    masks=[sum(1<<j for j in range(l) if L[j]>>i&1) for i in range(l)]
    masks += [(m<<l)|((1<<l)-1) for m in R]
    masks[y]^=1<<x;masks[x]|=1<<y
    h,U=D['P'];k,V=features(e,'P')
    c,UE=D['E'];d,VE=features(e,'E')
    bc,UC=D['C'];ec,VC=features(e,'C')
    assert (k,d,bc,ec)==(g(e-1),e%2,b,e)
    z=Counter()
    def weight(s,scale,side):
        for p,num in s.items():
            if p[0]=='a':
                _,u,v=p
                z[(v,u) if side=='L' else (u+l,v+l)]+=num*scale
    weight(U,k,'L');weight(UE,-d,'L')
    weight(V,h,'R');weight(VE,-c,'R')
    weight(UC,e,'L');weight(VC,b,'R')
    for i,j in product(range(r),range(l)):
        if i!=r-1 or j!=l-1:z[(i+l,j)]+=V['end',i]*U['end',j]
        z[(i+l,j)]+=VC['pe',i]*UC['se',j]+VC['se',i]*UC['pe',j]
    z[x,y]=b*e
    H=h*k-c*d+b*e
    arcs=[(u,v,z[u,v]) for u in range(n) for v in range(n) if masks[u]>>v&1]
    assert all(not load or masks[u]>>v&1 for (u,v),load in z.items())
    assert all(load>=0 for u,v,load in arcs)
    assert sum(load for u,v,load in arcs)==(n-1)*H
    return dict(order=n,H=H,b=b,e=e,rows=masks,unused=sum(not w for u,v,w in arcs),
                zero_arcs=[[u,v] for u,v,w in arcs if not w],arcs_uv_load=arcs)

if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('pair',nargs='*',help='b:e')
    p.add_argument('--check',action='store_true')
    a=p.parse_args()
    if a.check:
        for b,e in product(range(1,17),repeat=2):
            P=analyze(b,e);Q=profile(P['rows'])
            assert P['H']==Q['H']
            assert P['arcs_uv_load']==[arc[:3] for arc in Q['arcs_uv_load_reverse_load']]
        print("256 F1 tables agree with generic full load calculation")
        for e in range(1,17):
            P=analyze_bank((2,4,9,3,15),e);Q=profile(P['rows'])
            assert P['arcs_uv_load']==[arc[:3] for arc in Q['arcs_uv_load_reverse_load']]
        print("16 F2 tables agree with generic full load calculation")
    for item in a.pair:
        l,r=item.split(':')
        print(json.dumps(analyze_bank(tuple(map(int,l.split(','))),int(r))
                         if ',' in l else analyze(int(l),int(r))))
