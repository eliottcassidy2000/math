# enumerate G (ordinal + cyclic C3 joins) up to isomorphism via canonical modular-decomposition trees; emit rows
import sys
NMAX=int(sys.argv[1])
# a tree: ('v',) | ('L', children tuple) linear (ordinal chain, length>=2, no child is 'L') | ('C', (a,b,c)) rotation-canonical
def key(t): return repr(t)
def size(t):
    if t[0]=='v': return 1
    return sum(size(c) for c in t[1])
trees={1:[('v',)]}
def lin(a,b):
    ca=a[1] if a[0]=='L' else (a,); cb=b[1] if b[0]=='L' else (b,)
    return ('L',ca+cb)
def cyc(a,b,c):
    rots=[(a,b,c),(b,c,a),(c,a,b)]
    return ('C',min(rots,key=lambda r:key(r)))
for n in range(2,NMAX+1):
    S={}
    for a in range(1,n):
        for x in trees[a]:
            for y in trees[n-a]:
                t=lin(x,y); S[key(t)]=t
    for a in range(1,n-1):
        for b in range(1,n-a):
            c=n-a-b
            for x in trees[a]:
                for y in trees[b]:
                    for z in trees[c]:
                        t=cyc(x,y,z); S[key(t)]=t
    trees[n]=list(S.values())
    print('order',n,'G classes',len(trees[n]),file=sys.stderr)
def build(t):
    """return (n, rows) with vertices 0..n-1"""
    if t[0]=='v': return 1,[0]
    parts=[build(c) for c in t[1]]
    n=sum(p[0] for p in parts); rows=[0]*n; offs=[]; o=0
    for p in parts: offs.append(o); o+=p[0]
    for k,(m,r) in enumerate(parts):
        for i in range(m):
            rows[offs[k]+i]|=r[i]<<offs[k]
    if t[0]=='L':
        for k in range(len(parts)):
            for l in range(k+1,len(parts)):
                for i in range(parts[k][0]):
                    for j in range(parts[l][0]): rows[offs[k]+i]|=1<<(offs[l]+j)
    else:
        for k in range(3):
            l=(k+1)%3
            for i in range(parts[k][0]):
                for j in range(parts[l][0]): rows[offs[k]+i]|=1<<(offs[l]+j)
    return n,rows
for n in range(1,NMAX+1):
    for t in trees[n]:
        m,rows=build(t); print(m,*rows)
