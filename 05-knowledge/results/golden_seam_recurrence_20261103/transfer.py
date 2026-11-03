"""Hamilton paths of literal golden complements via a one-vertex seam."""
from collections import defaultdict
from functools import lru_cache
from itertools import combinations
from pathlib import Path
from importlib.util import module_from_spec, spec_from_file_location
import json

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
START, TAIL = -2,-1
STATES = (((START,TAIL),), ((START,0),), ((0,TAIL),), ((0,0),),
          ((START,TAIL),(0,0)), ((START,TAIL),(0,TAIL)))


def need(ok, msg):
    if not ok:
        raise RuntimeError(msg)


def chart(n,bits):
    rows=[0]*n
    j=0
    for a,b in combinations(range(n),2):
        forward=b==a+1 or bits>>j&1
        j+=b!=a+1
        if forward: rows[a]|=1<<b
        else: rows[b]|=1<<a
    return tuple(rows)


ATOMS={'a':(0,), 'b':chart(3,0), **{f'q{i}':chart(4,i) for i in (0,1,2,4,5)}}
COST={k:len(v) for k,v in ATOMS.items()}
POSITIVE={
    'a':((1,0,1),(2,0,1),(0,1,0)),
    'b':((1,2,2),(4,3,4),(4,2,3)),
    'q0':((1,4,4),(8,9,10),(8,8,9)),
    'q1':((3,3,3),(12,6,9),(12,5,8)),
    'q2':((5,3,4),(8,7,8),(6,6,7)),
    'q4':((3,6,6),(6,8,9),(6,5,6)),
    'q5':((5,4,5),(10,5,8),(8,4,5)),
}
# Twice the projection from the six physical states to the three coordinates.
PROJECTION=((2,0,0),(2,-1,1),(2,1,1),(2,0,2),(0,1,1),(0,1,1))


def complement(rows):
    n=len(rows)
    out=[(((1<<n)-1)^row)&~(1<<u) for u,row in enumerate(rows)]
    for u in range(n-1):
        need(rows[u]>>(u+1)&1,'legal physical path')
        out[u]|=1<<(u+1)
        out[u+1]&=~(1<<u)
    return tuple(out)


def reversed_rows(rows):
    n=len(rows)
    return tuple(sum(1<<b for b in range(n) if rows[n-1-a]>>(n-1-b)&1) for a in range(n))


def joined(word):
    rows=()
    for letter in word:
        atom=ATOMS[letter]
        k,d=len(rows),len(atom)
        rows=tuple(row|(((1<<d)-1)<<k) for row in rows)+tuple(row<<k for row in atom)
    return rows


def advance(state, rows, old=True):
    """Introduce a whole reversed target block; keep its last vertex active.

    Only the old boundary vertex0 has exceptional future incidence.
    Old anonymous tails point to every new vertex. A missing incoming edge
    is declared the unique global source only when its last chance expires.
    """
    d=len(rows)
    defects={(a+1,b+1) for a,b in combinations(range(d),2) if rows[b]>>a&1}
    if old:
        defects.add((0,1))
    last=[-1]*(d+1)
    for a,b in defects: last[a]=max(last[a],b)
    last[d]=d+1
    bank={state:1}
    for v in range(1,d+1):
        following=defaultdict(int)
        for components,weight in bank.items():
            incoming=[None]+[i for i,(_,t) in enumerate(components) if t==TAIL or (t,v) not in defects]
            outgoing=[None]+[i for i,(h,_) in enumerate(components) if h>=0 and (h,v) in defects]
            for left in incoming:
                for right in outgoing:
                    if left is not None and left==right: continue
                    head=components[left][0] if left is not None else v
                    tail=components[right][1] if right is not None else v
                    new=[item for i,item in enumerate(components) if i!=left and i!=right]+[(head,tail)]
                    expired=[]
                    starts=0
                    for h,t in new:
                        if h>=0 and last[h]<=v: h=START
                        if t>=0 and last[t]<=v: t=TAIL
                        starts+=h==START
                        expired.append((h,t))
                    if starts>1: continue
                    following[tuple(sorted(expired))]+=weight
        bank=following
    result=defaultdict(int)
    for state,weight in bank.items():
        renamed=tuple((0 if h==d else h,0 if t==d else t) for h,t in state)
        need(renamed in STATES,('closed seam state',renamed))
        result[renamed]+=weight
    return tuple(result[s] for s in STATES)


@lru_cache(None)
def matrices():
    mats={}
    initial={}
    for name,atom in ATOMS.items():
        rows=reversed_rows(complement(atom))
        mats[name]=tuple(advance(state,rows) for state in STATES)
        initial[name]=advance((),rows,False)
    return mats,initial


def multiply(v,M):
    return tuple(sum(v[i]*M[i][j] for i in range(len(v))) for j in range(len(v)))


def vector(word):
    if not word: return None
    mats,initial=matrices()
    v=initial[word[-1]]
    for letter in reversed(word[:-1]):
        v=multiply(v,mats[letter])
    return v


def count(word):
    v=vector(word)
    return 1 if v is None else sum(v[:4])


def fast_vector(word):
    v=(1,0,0)
    for letter in reversed(word):
        m=POSITIVE[letter]
        v=tuple(sum(v[i]*m[i][j] for i in range(3)) for j in range(3))
    return v


def verify_projection():
    mats,initial=matrices()
    for name,M in mats.items():
        P=POSITIVE[name]
        for i in range(6):
            for j in range(3):
                left=sum(M[i][k]*PROJECTION[k][j] for k in range(6))
                right=sum(PROJECTION[i][k]*P[k][j] for k in range(3))
                need(left==right,('all-word intertwining',name,i,j))
        for j in range(3):
            need(sum(initial[name][i]*PROJECTION[i][j] for i in range(6))==2*P[0][j],
                 ('all-word initialization',name,j))


def subset_count(rows):
    n=len(rows)
    if n==0: return 1
    dp=[[0]*n for _ in range(1<<n)]
    for u in range(n): dp[1<<u][u]=1
    for mask in range(1,1<<n):
        for u,w in enumerate(dp[mask]):
            if not w: continue
            available=rows[u]&~mask
            while available:
                bit=available&-available
                available-=bit
                dp[mask|bit][bit.bit_length()-1]+=w
    return sum(dp[-1])


def words(n):
    if n==0:
        yield ()
    for a,c in COST.items():
        if c<=n:
            for tail in words(n-c): yield (a,)+tail


def main():
    verify_projection()
    mats,initial=matrices()
    checks=0
    for n in range(11):
        for word in words(n):
            if not word: continue
            actual=subset_count(complement(joined(word)))
            need(count(word)==actual,('Hamilton recurrence',word,count(word),actual,vector(word)))
            need(fast_vector(word)[0]==actual,('positive recurrence',word))
            checks+=1
    print(json.dumps({'checks':checks,'states':STATES,'matrices':mats,'initial':initial},indent=2))


if __name__=='__main__': main()
