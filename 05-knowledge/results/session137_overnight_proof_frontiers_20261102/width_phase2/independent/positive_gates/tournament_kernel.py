"""Independent Session137 audit kernel; stdlib, exact integers only.

No producer module is imported. Tournament matrix A[i][j]=1 iff i beats j.
Code bit k corresponds to (i,j), i<j, in lexicographic order; set means i->j.
"""
from functools import lru_cache
from itertools import permutations, combinations

def require(condition, message):
    if not condition:
        raise RuntimeError(message)

def matrix(n, code):
    a = [[0]*n for _ in range(n)]
    for k,(i,j) in enumerate(combinations(range(n),2)):
        a[i][j] = (code>>k)&1
        a[j][i] = 1-a[i][j]
    return tuple(map(tuple,a))

def encoding(a):
    return sum(a[i][j]<<k for k,(i,j) in enumerate(combinations(range(len(a)),2)))

def induced(a, vs):
    return tuple(tuple(a[i][j] for j in vs) for i in vs)

def hp(a):
    """Backward endpoint recurrence, rather than producer's forward DP."""
    n = len(a)
    if not n:
        return 1
    incoming = [sum(a[i][j]<<i for i in range(n)) for j in range(n)]
    @lru_cache(None)
    def ending(s,j):
        t=s^(1<<j)
        if not t:
            return 1
        return sum(ending(t,i) for i in range(n) if (t&incoming[j])>>i&1)
    return sum(ending((1<<n)-1,j) for j in range(n))

def literal_hp(a):
    return sum(all(a[p[i]][p[i+1]] for i in range(len(a)-1))
               for p in permutations(range(len(a))))

def strong(a):
    n=len(a)
    if not n:
        return False
    for start in range(n):
        reached={start}; frontier=[start]
        while frontier:
            u=frontier.pop()
            for v in range(n):
                if a[u][v] and v not in reached:
                    reached.add(v); frontier.append(v)
        if len(reached)!=n:
            return False
    return True

def append(a,cut):
    """new root beats i iff i in cut."""
    n=len(a)
    return tuple(a[i]+(1-((cut>>i)&1),) for i in range(n))+(tuple((cut>>i)&1 for i in range(n))+(0,),)

def odd_cycles(a):
    """Independent subset/permutation enumerator, not directed DFS.

    A directed cyclic ordering is rotation-quotiented by fixing smallest vertex.
    Reversal is not quotiented: it is not a directed cycle in a tournament.
    """
    out=[];n=len(a)
    for size in range(3,n+1,2):
        for vs in combinations(range(n),size):
            for tail in permutations(vs[1:]):
                p=(vs[0],)+tail
                if all(a[p[i]][p[(i+1)%size]] for i in range(size)):
                    out.append((p,sum(1<<v for v in p)))
    return out

def cycle_packing_polynomial(cycles):
    """Enumerate vertex-disjoint cycle subfamilies directly."""
    counts=[0]*(1+max((mask.bit_count() for _,mask in cycles),default=0))
    def visit(start,used,size):
        counts[size]+=1
        for i in range(start,len(cycles)):
            mask=cycles[i][1]
            if not (mask&used):
                visit(i+1,used|mask,size+1)
    visit(0,0,0)
    while counts and counts[-1]==0:
        counts.pop()
    return counts

def exposed(a):
    """Direct old-word cut enumeration, endpoints and ordered exposed slots."""
    n=len(a);s=[0]*n;e=[0]*n;q=[[0]*n for _ in range(n)]
    for p in permutations(range(n)):
        valid=[bool(a[p[k]][p[k+1]]) for k in range(n-1)]
        if all(valid):
            s[p[0]]+=1;e[p[-1]]+=1
        for k in range(n-1):
            if all(valid[j] for j in range(n-1) if j!=k):
                q[p[k]][p[k+1]]+=1
    return s,e,q

def canon(a):
    """Score-sorted class code. Exact permutations within equal-score blocks."""
    from itertools import product
    groups={}
    for i,row in enumerate(a):groups.setdefault(sum(row),[]).append(i)
    return min(encoding(induced(a,sum((tuple(p) for p in blocks),())))
               for blocks in product(*(permutations(groups[d]) for d in sorted(groups))))

def class_bank(max_n):
    reps={0};yield 0,reps
    for n in range(1,max_n+1):
        reps={canon(append(matrix(n-1,c),cut))
              for c in reps for cut in range(1<<(n-1))}
        yield n,reps
