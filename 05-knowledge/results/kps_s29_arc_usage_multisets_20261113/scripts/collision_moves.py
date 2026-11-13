# classify arc-usage multiset collisions of merged classes (orders 7-9) by the moves of THM-5154 (A):
# (i) singleton redistribution among nonempty runs, (ii) converse flips of cores flanked by nonempty runs, (iii) global converse
# inputs (scratch dir): cls_n.txt (classgen: code H c3 strong SC aut ccode | ...), arcms_n.txt (scripts/arcms.c output)
import collections, itertools, functools
def decode(code,n):
    out=[0]*n; k=0
    for i in range(n):
        for j in range(i+1,n):
            if code>>k&1: out[i]|=1<<j
            else: out[j]|=1<<i
            k+=1
    return out
def sccs_ordered(out):
    n=len(out); reach=[0]*n
    for v in range(n):
        seen=fr=1<<v
        while fr:
            nf=0
            for w in range(n):
                if fr>>w&1: nf|=out[w]
            nf&=~seen; seen|=nf; fr=nf
        reach[v]=seen
    cs=[]; done=0
    for v in range(n):
        if done>>v&1: continue
        c=sum(1<<w for w in range(n) if reach[v]>>w&1 and reach[w]>>v&1); cs.append(c); done|=c
    def beats(c,d):
        v=(c&-c).bit_length()-1; w=(d&-d).bit_length()-1; return out[v]>>w&1
    cs.sort(key=functools.cmp_to_key(lambda c,d: -1 if beats(c,d) else 1)); return cs
def sub(out,mask):
    idx=[v for v in range(len(out)) if mask>>v&1]; pos={v:i for i,v in enumerate(idx)}
    return tuple(sum(1<<pos[w] for w in idx if out[v]>>w&1) for v in idx)
def iso(A,B):
    n=len(A)
    if sorted(bin(x).count('1') for x in A)!=sorted(bin(x).count('1') for x in B): return False
    return any(all(((A[i]>>j)&1)==((B[p[i]]>>p[j])&1) for i in range(n) for j in range(n) if i!=j) for p in itertools.permutations(range(n)))
def conv(A):
    n=len(A); return tuple(sum(1<<v for v in range(n) if A[v]>>w&1) for w in range(n))
def signature(code,n):
    out=decode(code,n); cores=[]; runs=[0]
    for c in sccs_ordered(out):
        if bin(c).count('1')==1: runs[-1]+=1
        else: cores.append(sub(out,c)); runs.append(0)
    return cores,runs
def equivalent(A,B):
    (ca,ra),(cb,rb)=A,B
    if len(ca)!=len(cb) or [r>0 for r in ra]!=[r>0 for r in rb] or sum(ra)!=sum(rb): return False
    return all(iso(x,y) or (ra[i]>0 and ra[i+1]>0 and iso(conv(x),y)) for i,(x,y) in enumerate(zip(ca,cb)))
for n in (7,8,9):
    cc={}
    for line in open('cls_%d.txt'%n):
        p=line.split('|')[0].split(); cc[int(p[0])]=int(p[6])
    ms={}
    for line in open('arcms_%d.txt'%n):
        left,right=line.split('|'); a=left.split(); ms[int(a[0])]=(int(a[2]),tuple(map(int,right.split())))
    reps=sorted({min(c,cc[c]) for c in cc})
    M=collections.defaultdict(list)
    for r in reps: M[ms[r]].append(r)
    groups=[g for g in M.values() if len(g)>=2]
    ok=sum(1 for g in groups if all(equivalent(signature(g[0],n),signature(r,n)) or equivalent(signature(g[0],n),signature(cc[r],n)) for r in g[1:]))
    print('order',n,'collision groups',len(groups),'explained by the moves',ok)
