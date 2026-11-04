"""Exact ordered path-bin counts for an actual small-feedback tournament.

No subset-size exponential table. The supplied feedback deletion must be
transitive; endpoint constraints promote their named vertices to feedback.
"""
from functools import lru_cache
from itertools import permutations


def need(ok, message):
    if not ok:
        raise RuntimeError(message)


def compositions(n, k):
    if k == 1:
        yield (n,)
    else:
        for first in range(n+1):
            for rest in compositions(n-first,k-1):
                yield (first,)+rest


@lru_cache(None)
def count(rows, feedback, bins=1, starts=(), ends=()):
    """Count ordered possibly empty bins, with (bin,vertex) constraints."""
    n=len(rows)
    need(bins>=1,'positive actual bin count')
    need(all(0<=v<n for v in feedback),'feedback indices')
    need(all(0<=row<1<<n and not(row>>v&1) for v,row in enumerate(rows)),
         'tournament row bounds and diagonal')
    need(all(((rows[v]>>u)&1)+((rows[u]>>v)&1)==1
             for v in range(n) for u in range(v)), 'tournament pairs')
    start=dict(starts);end=dict(ends)
    need(len(start)==len(starts) and len(end)==len(ends),'one constraint per bin and direction')
    need(all(0<=b<bins and 0<=v<n for b,v in starts+ends),'endpoint constraint')
    ports=tuple(sorted(set(feedback)|set(start.values())|set(end.values())))
    rest=[v for v in range(n) if v not in ports]
    chain=tuple(sorted(rest,key=lambda v:-sum(rows[v]>>u&1 for u in rest)))
    need(all(rows[u]>>v&1 for i,u in enumerate(chain) for v in chain[i+1:]),
         'feedback deletion must be transitive')
    total=0
    for order in permutations(ports):
        for sizes in compositions(len(ports),bins):
            words=[];pos=0;okay=True
            for b,size in enumerate(sizes):
                w=order[pos:pos+size];pos+=size
                if b in start and (not w or w[0]!=start[b]):okay=False;break
                if b in end and (not w or w[-1]!=end[b]):okay=False;break
                words.append(w)
            if not okay:continue
            gaps=[]
            for b,w in enumerate(words):
                for j in range(len(w)+1):
                    left=w[j-1] if j else None
                    right=w[j] if j<len(w) else None
                    forbidden=(j==0 and b in start) or (j==len(w) and b in end)
                    gaps.append((left,right,forbidden))
            powers=[3**j for j in range(len(gaps))]
            states={0:1}
            for v in chain:
                nxt={}
                for state,weight in states.items():
                    for j,(left,right,forbidden) in enumerate(gaps):
                        if forbidden:continue
                        old=state//powers[j]%3
                        if old==0 and left is not None and not(rows[left]>>v&1):continue
                        new=1 if right is None or rows[v]>>right&1 else 2
                        target=state+(new-old)*powers[j]
                        nxt[target]=nxt.get(target,0)+weight
                states=nxt
            for state,weight in states.items():
                okay=True
                for j,(left,right,_) in enumerate(gaps):
                    digit=state//powers[j]%3
                    if digit==2 or (digit==0 and left is not None and right is not None
                                    and not(rows[left]>>right&1)):
                        okay=False;break
                if okay:total+=weight
    return total


def induced(rows, vertices):
    return tuple(sum(1<<j for j,v in enumerate(vertices) if rows[u]>>v&1)
                 for u in vertices)


def deletion_h(rows,feedback,root):
    vertices=[v for v in range(len(rows)) if v!=root]
    if not vertices:return 1
    return count(induced(rows,vertices),tuple(vertices.index(v) for v in feedback if v!=root))


def pair_fields(rows,feedback,r,s):
    """Match inherited endpoint.cpp order H,P2,[E1,S1,q,E2,S2]*2,W2."""
    rows=tuple(rows);feedback=tuple(feedback)
    out=[count(rows,feedback),count(rows,feedback,2)]
    for p in (r,s):
        out += [count(rows,feedback,ends=((0,p),)),
                count(rows,feedback,starts=((0,p),)),
                deletion_h(rows,feedback,p),
                count(rows,feedback,2,ends=((0,p),)),
                count(rows,feedback,2,starts=((0,p),))]
    for p in (r,s):
        for q in (r,s):
            out.append(0 if p==q else count(rows,feedback,2,starts=((1,q),),ends=((0,p),)))
    return out


def nine_fields(rows,feedback,w):
    """h,s,q,F,e,P2,S2,E3,Q2 for one-polarity R/N component moves."""
    rows=tuple(rows);feedback=tuple(feedback)
    vs=[v for v in range(len(rows)) if v!=w]
    q2=count(induced(rows,vs),tuple(vs.index(v) for v in feedback if v!=w),2) if vs else 1
    return [count(rows,feedback),count(rows,feedback,starts=((0,w),)),
            deletion_h(rows,feedback,w),count(rows,feedback,2,ends=((0,w),)),
            count(rows,feedback,ends=((0,w),)),count(rows,feedback,2),
            count(rows,feedback,2,starts=((0,w),)),count(rows,feedback,3,ends=((0,w),)),q2]
