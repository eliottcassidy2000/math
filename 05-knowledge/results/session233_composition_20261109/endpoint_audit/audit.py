"""Independent verifier: pair sets, permutation/cut covers, recursive endpoint counts.
No imports of producer code and no assert statements. --fail tests live guards.
"""
import itertools as it
import functools
import json
import sys

def require(truth, label):
    if not truth:
        raise RuntimeError('AUDIT FAILURE: '+label)

if '--fail' in sys.argv:
    require(False, 'deliberate failure sentinel')

def base(n):
    return {(i,j) for i in range(n) for j in range(i+1,n) if (i,j)!=(0,n-1)} | {(n-1,0)}

def extension(n,N):
    arcs=base(n)
    for v in range(n):
        arcs.add((v,n) if 1<=v<n-1 and (N//2**(v-1))%2 else (n,v))
    return arcs

def literal(n,arcs):
    ends=[0]*n; response=[0]*n
    for p in it.permutations(range(n)):
        bad={i+1 for i in range(n-1) if (p[i],p[i+1]) not in arcs}
        if not bad:
            ends[p[-1]]+=1
        # A cut c gives first (rooted) path p[:c], remaining path p[c:].
        # First path nonempty; second may be empty. Only cut adjacency is ignored.
        for c in range(1,n+1):
            if bad <= {c}:
                response[p[0]]+=1
    return ends,response

def recursive(n,arcs):
    @functools.cache
    def e(S,v):
        rest=S-{v}
        return 1 if not rest else sum(e(rest,u) for u in rest if (u,v) in arcs)
    V=frozenset(range(n))
    return [e(V,v) for v in V]

def strong(n,arcs):
    # Independent Floyd-Warshall transitive closure.
    r=[[i==j or (i,j) in arcs for j in range(n)] for i in range(n)]
    for k in range(n):
        for i in range(n):
            for j in range(n):
                r[i][j] = r[i][j] or (r[i][k] and r[k][j])
    return all(all(row) for row in r)

def aut(n,arcs,root=None):
    count=0
    for p in it.permutations(range(n)):
        if root is not None and p[root]!=root:continue
        if all((p[i],p[j]) in arcs for i,j in arcs):count+=1
    return count

counts={'base_literal':0,'selector_literal':0,'selector_recursive':0,'arbitrary_extension_response':0}
for n in range(4,8):
    a=base(n); e,B=literal(n,a)
    want=[1]+[2**(i-1) for i in range(1,n-1)]+[1]
    require(e==want,'base endpoint vector')
    require(aut(n,a)==1,'base rigidity')
    require(strong(n,a),'base strongness')
    counts['base_literal']+=1
    if n<=6:
        rows=[]
        for N in range(1,2**(n-2)):
            u=extension(n,N); eu,Bu=literal(n+1,u)
            require(eu==recursive(n+1,u),'literal vs recursive endpoints')
            require(eu[-1]==N,'selector endpoint')
            require(aut(n+1,u,n)==1,'actual root stabilizer')
            require(strong(n+1,u),'selector strongness')
            expected=sum(e)+sum(B[v] for v in range(n) if (n,v) in u)
            require(Bu[n]==expected,'literal response sidecar')
            rows.append([N,sum(eu),Bu[n]])
            counts['selector_literal']+=1
        if n==6:
            require(B==[31,69,44,28,18,17],'reported base response')
            require([r[1] for r in rows]==[27,41,33,59,51,65,39,85,71,91,59,109,77,97,47],'reported H profile')
            require([r[2] for r in rows]==[155,180,111,196,127,152,83,206,137,162,93,178,109,134,65],'reported B profile')

for n in range(4,10):
    for N in range(1,2**(n-2)):
        u=extension(n,N)
        require(recursive(n+1,u)[-1]==N,'larger selector endpoint')
        require(strong(n+1,u),'larger selector strongness')
        counts['selector_recursive']+=1

witnesses=[]
for N,H in [(1,9),(7,23),(21,187),(37,461),(101,1573),(257,9899)]:
    n=max(4,N.bit_length()+2);u=extension(n,N); e=recursive(n+1,u)
    require(e[n]==N,'named endpoint')
    require(sum(e)==H,'named total')
    require(n+1==max(5,N.bit_length()+3),'integer order bound')
    witnesses.append({'N':N,'order':n+1,'H':sum(e)})

# The sidecar identity is asserted for ANY extension: exhaust all labelled
# tournaments A of order 1..4 and all outgoing sets of x, strong or otherwise.
for n in range(1,5):
    pairs=list(it.combinations(range(n),2))
    for mask in range(2**len(pairs)):
        a={(i,j) if mask>>k&1 else (j,i) for k,(i,j) in enumerate(pairs)}
        e,B=literal(n,a)
        for outmask in range(2**n):
            u=a|{(n,v) if outmask>>v&1 else (v,n) for v in range(n)}
            eu,Bu=literal(n+1,u)
            require(Bu[n]==sum(e)+sum(B[v] for v in range(n) if outmask>>v&1),'arbitrary extension identity')
            counts['arbitrary_extension_response']+=1
print(json.dumps({'status':'PASS','counts':counts,'witnesses':witnesses},sort_keys=True,indent=2))
