#!/usr/bin/env python3
from collections import Counter
from functools import lru_cache

CORE=(974,1020,1016,976,929,201,528,832,608,32)
MARKS=(1,9)
CYCLE=(0,1,2,3,4,7,8,6,9,5)


def require(test,message):
    if not test:
        raise RuntimeError(message)


def neighbors(mask):
    while mask:
        b=mask&-mask
        yield b.bit_length()-1
        mask-=b


def compressed(rows, marks):
    n=len(rows); a,b=marks
    unmarked=((1<<n)-1)^(1<<a)^(1<<b)

    @lru_cache(maxsize=None)
    def visit(seen, x, i, j):
        result=Counter({(i,j):1}) if seen==unmarked and i and j else Counter()
        for y in neighbors(rows[x]):
            if y==a and i<n-1:
                result.update(visit(seen,y,i+1,j))
            elif y==b and j<n-1:
                result.update(visit(seen,y,i,j+1))
            elif ((unmarked^seen)>>y)&1:
                result.update(visit(seen|(1<<y),y,i,j))
        return result
    ans=Counter()
    for x in range(n):
        ans.update(visit(unmarked&(1<<x),x,int(x==a),int(x==b)))
    return dict(ans)


def inflate(rows, marks, m, k):
    sizes=[1]*len(rows)
    sizes[marks[0]]=m
    sizes[marks[1]]=k
    ends=[0]
    for s in sizes:
        ends.append(ends[-1]+s)
    blocks=[((1<<ends[i+1])-(1<<ends[i])) for i in range(len(rows))]
    g=[]
    for i,row in enumerate(rows):
        outside=sum(blocks[j] for j in neighbors(row))
        for x in range(ends[i],ends[i+1]):
            g.append(outside+((1<<ends[i+1])-(1<<(x+1))))
    return g


def strong(rows):
    All=(1<<len(rows))-1
    for x in range(len(rows)):
        done=0; Q=1<<x
        while Q:
            b=Q&-Q; Q-=b
            done|=b
            Q|=rows[b.bit_length()-1]&~done
        if done!=All:
            return False
    return True


def H(rows):
    n=len(rows)
    h=[[0]*n for _ in range(1<<n)]
    for v in range(n):
        h[1<<v][v]=1
    for S, entries in enumerate(h):
        for i,value in enumerate(entries):
            if value:
                for j in neighbors(rows[i]&~S):
                    h[S|(1<<j)][j]+=value
    return sum(h[-1])


if __name__ == "__main__":
    require(all(CORE[a]>>b&1 for a,b in zip(CYCLE,CYCLE[1:]+CYCLE[:1])),"cycle")
    sig=compressed(CORE,MARKS)
    require(sig=={(1,1):899,(1,2):551,(2,1):418,(2,2):256},"compressed profile")
    print("compressed signature:",sorted(sig.items()))
    for p in (6,7,8,31):
        m,k=p-5,p-4
        g=inflate(CORE,MARKS,m,k)
        E=(1<<(p-1))*((1<<p)-1)
        x,y=1<<m,1<<k
        v=899+551*(y-2)+418*(x-2)+256*(x-2)*(y-2)
        require(strong(g),"inflation is not strong")
        require(v==E-15,"perfect-offset formula")
        if p<=8:
            require(H(g)==v,"literal subset-DP check")
        print(p,len(g),v)
