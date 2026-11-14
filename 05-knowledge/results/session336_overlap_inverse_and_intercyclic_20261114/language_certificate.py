from itertools import product,combinations
from collections import deque

alphabet=tuple(range(4))
pairs=tuple(combinations(range(4),2))
bad={w for w in product(alphabet,repeat=4) for i,j in pairs
     if (lambda k,l: w[i]&2 and not w[j]&2 and w[k]&1 and not w[l]&1)(*(sorted(set(range(4))-{i,j})))}
prefix={w[:i] for w in bad for i in range(4)}
monitor=lambda seen,d: (seen | {w+(d,) for w in seen if w+(d,) in prefix})
# state flags: partial bad-subsequences, last (initial -1), X-then-nonX,
# Y-then-nonY, reach Y->X, seen bits, no deletion & <=1 deletion states
start=(frozenset({()}),-1,False,False,False,0,0,frozenset())
def plain(q,d):
    return d if q is not None and d!=1 and q<=d else None
def tr(s,d):
    if s is None:return None
    seen,last,ix,iy,ret,bits,legal,gone=s
    if last==-1 and d==0:return None
    if any(w+(d,) in bad for w in seen):return None
    return (frozenset(monitor(seen,d)), d, bool(ix or (bits&2 and not d&2)),
           bool(iy or (bits&1 and not d&1)), ret or d==1 or (d==0 and bool(bits&4)),
           bits|d|(4 if d==3 else 0),plain(legal,d),
           frozenset(y for x in gone if (y:=plain(x,d)) is not None)|({legal} if legal is not None else set()))
def acc(s):
    if s is None:return False
    seen,last,ix,iy,ret,bits,legal,gone=s
    return last in (0,1,2) and ix and iy and ret and legal is None and not gone

P=["1+20","321+","1+21+","1+3+0","2+10","30+1+","312+","310","320"]
def parse(s):return tuple((int(d),i+1<len(s) and s[i+1]=='+') for i,d in enumerate(s) if d!='+')
P=list(map(parse,P))
q0=frozenset((i,0) for i in range(len(P)))
def nr(s,d):
    t=set()
    for i,p in s:
        a=P[i]
        if p<len(a) and a[p][0]==d:t.add((i,p+1))
        if p and a[p-1]==(d,True):t.add((i,p))
    return frozenset(t)
def na(s):return any(j==len(P[i]) for i,j in s)
todo=deque([(start,q0,"")]);seen={(start,q0)}
while todo:
    q,r,w=todo.popleft()
    assert acc(q)==na(r),w
    for d in alphabet:
        a,b=tr(q,d),nr(r,d)
        if (a,b) not in seen:
            seen.add((a,b));todo.append((a,b,w+str(d)))
print('all-length word language equality over',len(seen),'reachable paired states;',
      len(bad),'forbidden length-four patterns')
