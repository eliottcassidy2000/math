import itertools, json, math, hashlib
from pathlib import Path

def ck(x):
    if not x: raise RuntimeError('independent control failed')
def partitions(r,k):
    out=[]
    def rec(t):
        if len(t)==r: out.append(t); return
        for v in range(min(k,1+max(t,default=-1)+1)):rec(t+(v,))
    rec(())
    return out

def stirling(n):
    row=[1]
    for _ in range(n):row=[0]+[j*(row[j] if j<len(row) else 0)+row[j-1] for j in range(1,len(row)+1)]
    return row

def chrom(n,edges,k):
    # Independent inclusion-exclusion over monochromatic edge events.
    val=0
    for mask in range(1<<len(edges)):
        parent=list(range(n))
        def root(x):
            while parent[x]!=x:x=parent[x]
            return x
        for i,(a,b) in enumerate(edges):
            if mask>>i&1:parent[root(a)]=root(b)
        val+=(-1 if mask.bit_count()%2 else 1)*k**len({root(x) for x in range(n)})
    return val

def sig(n,edges,k,p):
    return sum(all(c[a]!=c[b] for a,b in edges) for tail in itertools.product(range(k),repeat=n-len(p)) for c in [p+tail])
counts=[];products=0;indicators=0;distinguish=0
for r in range(6):
 for k in range(1,7):
    ps=partitions(r,k);ck(len(ps)==sum(stirling(r)[:k+1]));counts.append([r,k,len(ps)])
    pairs=list(itertools.combinations(range(r),2))
    for p in ps:
      within=[e for e in pairs if p[e[0]]==p[e[1]]]
      cross=[e for e in pairs if p[e[0]]!=p[e[1]]]
      for q in ps:
        # Expand each partition indicator as signed boundary-test monomials.
        expanded=0
        if all(q[a]!=q[b] for a,b in cross):
          for mask in range(1<<len(within)):
            if all(not(mask>>i&1) or q[a]!=q[b] for i,(a,b) in enumerate(within)):
              expanded+=(-1)**mask.bit_count()
        ck(expanded==int(p==q));indicators+=1
        if p!=q:ck(any((p[a]==p[b])!=(q[a]==q[b]) for a,b in pairs));distinguish+=1
# All fragments with one private vertex, all r<=3, all k<=4.
for r in range(4):
 pairs=list(itertools.combinations(range(r+1),2))
 fam=[tuple(e for i,e in enumerate(pairs) if mask>>i&1) for mask in range(1<<len(pairs))]
 for k in range(1,5):
    ps=partitions(r,k);ss=[[sig(r+1,e,k,p) for p in ps] for e in fam]
    cache={}
    for i,g in enumerate(fam):
      for j,h in enumerate(fam):
        joined=tuple(sorted(set(g)|{tuple(r+1 if x==r else x for x in e) for e in h}))
        if joined not in cache:cache[joined]=chrom(r+2,joined,k)
        expect=sum(math.factorial(k)//math.factorial(k-len(set(p)))*ss[i][a]*ss[j][a] for a,p in enumerate(ps))
        ck(expect==cache[joined]);products+=1
for k in range(1,9):
 ck(chrom(3,[(0,1)],k)==k*k*(k-1));ck(chrom(3,[(0,1),(0,2)],k)==k*(k-1)**2)
# Explicit loops give zero, repeated edges do not alter counts.
for k in range(1,5):
 ck(chrom(2,[(0,0)],k)==0)
 ck(chrom(2,[(0,1),(0,1)],k)==chrom(2,[(0,1)],k))
files=['PORT_SIGNATURES.md','verify_ports.py']
print(json.dumps({'status':'PASS','partition_counts':counts,'expanded_indicator_checks':indicators,'single_edge_separations':distinguish,'gluing_pairs':products,'source_hashes':{f:hashlib.sha256((Path('session25-ports')/f).read_bytes()).hexdigest() for f in files}},indent=2))
