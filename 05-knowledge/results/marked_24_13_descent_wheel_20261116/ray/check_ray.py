from collections import defaultdict, deque
from itertools import product

ORIGIN=(0,0)
DEBTS=(0,-13,24,11)
PATHS={
  ((0,0),'133220'):(-8,1),
  ((0,0),'311103311002'):(0,0),
  ((-8,1),'032220'):(0,0),
  ((-8,1),'2101012022331'):(-11,11),
  ((-11,11),'022230022331'):(-11,11),
  ((-11,11),'200113'):(-3,10),
  ((-3,10),'301113'):(-11,11),
  ((-3,10),'1232321311002'):(0,0),
}
P='133220032220'; Q='311103311002'
U='1332202101012022331'; V='022230022331'
W='200113301113'; Z='2001131232321311002'

def transition(pos,z):
    f,c=pos
    d=DEBTS[z]
    F,C=2*f+d,c+11-d
    if -23<=F<=12 and not(C%2):
        return F,C//2

def graph():
    front=deque([ORIGIN]); ed={}
    while front:
        p=front.popleft()
        if p in ed: continue
        m={str(z):q for z in range(4) if (q:=transition(p,z)) is not None}
        ed[p]=m
        front.extend(m.values())
    core={ORIGIN}
    while (b:=core|{p for p,m in ed.items() if core.intersection(m.values())})!=core:
        core=b
    return len(ed),{p:{z:q for z,q in ed[p].items() if q in core} for p in core}

def compressed(g):
    ret={}
    for p,m in g.items():
        if len(m)<=1: continue
        for key,q in m.items():
            while len(g[q])==1:
                z,q=next(iter(g[q].items())); key+=z
            ret[p,key]=q
    return ret

def address(w):
    a=b=c=d=0
    for j,z in enumerate(w):
        i,k=divmod(int(z),2)
        a,b=(a<<1)+i,(b<<1)+k
        c,d=c+((1-i)<<j),d+((1-k)<<j)
    return a,b,c,d

def contracted_gap(w):
    tape=[divmod(int(z),2) for z in w]
    a,b=zip(*tape)
    def cross(j,i):
        return (1-a[i])*b[j]-(1-b[i])*a[j]
    return (b[0]+b[-1]-a[0]-a[-1]
         +sum((1<<(i-j-1))*cross(j,i) for i in range(len(w)) for j in range(i))
         +sum(cross(j+1,j) for j in range(len(w)-1)))

def compressed_score(w,state):
    prev=None; score=0
    for z in w:
        a,b=divmod(int(z),2); f,_=state
        state=transition(state,int(z))
        score-=a*state[1]+f*(1-a)
        if prev is not None:
            x,y=prev
            score+=13*((1-x)*b-(1-y)*a)
        prev=a,b
    return state,score

if __name__=='__main__':
    states,g=graph()
    assert (states,len(g))==(66,58)
    assert compressed(g)==PATHS
    assert compressed_score(U,ORIGIN)==((-11,11),-2)
    assert compressed_score(V,(-11,11))==((-11,11),-26)
    assert compressed_score(W,(-11,11))==((-11,11),26)
    assert compressed_score(Z,(-11,11))==(ORIGIN,2)
    assert [contracted_gap(w) for w in (P,Q)]==[2,-2]
    for size in range(6):
        for loops in product((V,W),repeat=size):
            w=U+''.join(loops)+Z
            a,b,c,d=address(w)
            assert 24*a==13*b and 24*c==13*d
            assert contracted_gap(w)==2*(loops.count(W)-loops.count(V))
    print('66 reachable states, 58 on return walks, 4 branch states.')
    print('First returns: P,Q; or U{V,W}*Z. Both loops weighted ±2.')

    pell='13221321333110313000202102220231121001012103212323321220131113213133302032200021021102'
    a,b,c,d=address(pell)
    assert len(pell)==86 and (a,b)==(c,d)
    assert 324*a==187*b and contracted_gap(pell)==0
    assert all(b*b-3*a*a==69 for b,a in [(9,2),(24,13),(87,50),(324,187),(1209,698)])
    assert 2*23**12*25**7+23**19<25**19 and 25**8<2*23**8
