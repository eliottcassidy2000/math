"""Build and independently check the signed residue board."""
from collections import Counter, defaultdict
from pathlib import Path
from slot_interface_check import BRANCHES, read_board
import sys


def compressed(board):
    _, raw=read_board(board)
    mark=raw[0][7]
    D=[{y:t for y,t in ds.items() if all(a!=b for a,b in zip(t,mark))} for ds in raw[1:]]
    cap=[[0]*28 for _ in range(3)]
    for g in range(3):
        for s in {t[g] for r in D for t in r.values()}:
            cap[g][s%28]|=1<<(s//28)
    options=[]
    for h,row in zip(BRANCHES[1:],D):
        bucket={}
        special=h in (27,49)
        for y,c in row.items():
            r=y if special else y%28
            b=0 if special else y//28
            rep=tuple((s+28*b)%56 for s in c)
            if r not in bucket:bucket[r]=[rep,0]
            assert bucket[r][0]==rep
            bucket[r][1]|=1<<b
        options.append([(r,*c,m) for r,(c,m) in bucket.items()])
    assert sum(map(len,D))==1741 and sum(map(len,options))==1064
    assert Counter(o[4].bit_count() for r in options for o in r)=={1:387,2:677}
    assert Counter(m.bit_count() for g in cap for m in g)=={0:1,1:10,2:73}
    return D,cap,options


def check_profile(board,profile):
    D,cap,O=compressed(board)
    pick=list(map(int,Path(profile).read_text().split()))
    assert len(pick)==52
    T=[o[j] for o,j in zip(O,pick)]
    slots=defaultdict(list)
    G=[[] for _ in range(53)]
    def add(i,j,b):
        G[i].append((j,b));G[j].append((i,b))
    f=0
    for i,(_,a,b,c,mask) in enumerate(T):
        if mask!=3:add(i,52,mask//2);f+=1
        for k,s in enumerate((a,b,c)):slots[k,s%28].append((i,s//28))
    for g in range(3):
        for r,m in enumerate(cap[g]):
            P=slots[g,r]
            assert len(P)==m.bit_count()
            if len(P)==1:
                i,b=P[0];add(i,52,b^(m//2))
            if len(P)==2:
                (i,a),(j,b)=P;add(i,j,1^a^b)
    seen={}
    n=0;bad=[]
    for i in [52,*range(52)]:
        if i in seen:continue
        n+=1
        Q=[i];seen[i]=0;broken=False
        for a in Q:
            for b,x in G[a]:
                v=seen[a]^x
                if b not in seen:
                    seen[b]=v;Q.append(b)
                elif seen[b]!=v:broken=True
        if broken:bad.append(tuple(Q))
    print("components",n,"pinned edges",f,"cycle rank",30+f+n,
          "lift count",0 if bad else 1<<(n-1))
    if not bad:
        y=[r if BRANCHES[i+1] in (27,49) else r+28*seen[i] for i,(r,*_) in enumerate(T)]
        P=[D[i][j] for i,j in enumerate(y)]
        assert all(len({x[k] for x in P})==52 for k in range(3))
        print("Physical word:",7,*y)


if __name__=='__main__':
    if len(sys.argv)==3:check_profile(*sys.argv[1:])
    else:
        _,cap,options=compressed(sys.argv[1])
        print(*sum(cap,[]))
        for O in options:
            print(len(O))
            for v in O:print(*v)
