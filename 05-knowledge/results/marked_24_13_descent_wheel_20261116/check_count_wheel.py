from itertools import product
from collections import deque

P='133220032220'; Q='311103311002'
U='1332202101012022331'; V='022230022331'; W='200113301113'; Z='2001131232321311002'
ATOM=[(P,1),(Q,-1),(U+Z,0),(U+V+Z,-1),(U+W+Z,1)]

# Two copies of a three-run automaton, one for each apex order.
# State (f,l,m): run0 ends illegally into first apex (f);
# run2 is live (l); run1 empty or ends legally/illegally into last apex
# (m=0,1,2). All runs are increasing in spine order.
ST=list(product(range(2),range(2),range(3)))
def step(old,a,b):
    v=[0]*12
    for f,l,m in ST:
        val=old[6*f+3*l+m]
        v[6*a+3*l+m]+=val
        if a or m:v[6*f+3*l+1+b]+=val
        if b or l:v[6*f+3+m]+=val
    return tuple(x%120 for x in v)

def scan(old,w,swap):
    for z in w:
        a,b=divmod(int(z),2)
        old=step(old,b,a) if swap else step(old,a,b)
    return old

ids=[tuple(int(i==j) for i in range(12)) for j in range(12)]
def prepare(w):
    return [tuple(tuple(ij for ij in enumerate(col) if ij[1]) for col in map(lambda x:scan(x,w,sw),ids)) for sw in (False,True)]
def apply(state,M):
    dest=[]
    for cur,cols in zip(state,M):
        out=[0]*12
        for x,col in zip(cur,cols):
            if x:
                for j,v in col:out[j]+=x*v
        dest.append(tuple(x%120 for x in out))
    return tuple(dest)
def read(state):
    x,y=state
    return (x[0]+x[3]+x[1]+x[4]+y[1]+y[4])%120,(y[0]+y[3]+y[1]+y[4]+x[1]+x[4])%120

def bracket_witnesses(state,limit=6):
    seen={}
    mats=[prepare(w) for w,d in ATOM]
    todo=deque([(state,0,'')])
    while todo and len(seen)<60:
        s,d,seq=todo.popleft()
        if d==1:
            above,below=read(s)
            assert (above-below)%120==2
            assert below%2==1
            seen.setdefault(below,seq)
        if len(seq)<limit:
            for i,M in enumerate(mats):
                dn=d+ATOM[i][1]
                if abs(1-dn)>limit-len(seq)-1:continue
                todo.append((apply(s,M),dn,seq+str(i)))
    return seen
def main():
    base=(ids[0],ids[0])
    assert len(bracket_witnesses(base))==60
    pump=prepare(P+Q)
    start=base
    prev={}
    for j in range(122):
        if start in prev:
            t0=prev[start];period=j-t0
            assert (t0,period)==(1,120)
            break
        prev[start]=j
        start=apply(start,pump)
    else:
        raise AssertionError('expected repeat')
    wheel=bracket_witnesses(start)
    assert set(wheel)==set(range(1,120,2))
    # Check all returned words from their explicit spine descriptions.
    for r,seq in sorted(wheel.items()):
        w=''.join(ATOM[int(i)][0] for i in seq)
        assert sum(ATOM[int(i)][1] for i in seq)==1
        direct=tuple(scan(base[j],P+Q+w,bool(j)) for j in (0,1))
        assert read(direct)==((r+2)%120,r)
        print(r,seq)
    print('All 60 residues verified; neutral-prefix state repeats from t=1 with period 120.')
if __name__=='__main__':main()
