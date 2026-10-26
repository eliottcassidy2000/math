"""Exact tournament primitives; row i contains its outgoing neighbours."""
from collections import defaultdict
from functools import lru_cache
from itertools import permutations, product, chain

def require(ok, message):
    if not ok:
        raise RuntimeError(message)

def edge(t, i, j): return (t[i] >> j) & 1
def converse(t): return tuple(((1 << len(t))-1) ^ r ^ (1 << i) for i, r in enumerate(t))
def relabel(t, p): return tuple(sum(edge(t, p[i], p[j]) << j for j in range(len(t))) for i in range(len(t)))
def sigma(t): return relabel(converse(t), tuple(reversed(range(len(t)))))
def complement(t, p=None):
    p = tuple(range(len(t))) if p is None else p
    rows = list(converse(t))
    for a, b in zip(p, p[1:]):
        require(edge(t,a,b), 'complement needs a directed path')
        rows[a] ^= 1 << b; rows[b] ^= 1 << a
    return tuple(rows)

def presentation(n, mask):
    t = [0]*n; k = 0
    for i in range(n):
        for j in range(i+1,n):
            a = j == i+1 or (mask >> k)&1
            if j > i+1: k += 1
            if a: t[i] |= 1 << j
            else: t[j] |= 1 << i
    return tuple(t)

def free_mask(t):
    require(all(edge(t,i,i+1) for i in range(len(t)-1)), 'not a presentation')
    return sum(edge(t,i,j) << k for k,(i,j) in enumerate(( (i,j) for i in range(len(t)) for j in range(i+2,len(t)) )))

@lru_cache(None)
def canonical(t, roots=()):
    fixed = list(dict.fromkeys(roots)); pattern = tuple(fixed.index(v) for v in roots)
    cells = defaultdict(list)
    for v,r in enumerate(t):
        if v not in fixed: cells[r.bit_count()].append(v)
    choices = [tuple(permutations(cells[d])) for d in sorted(cells)]
    model = min(relabel(t, tuple(fixed)+tuple(chain.from_iterable(c))) for c in product(*choices))
    return model, pattern

def is_strong(t):
    # Input here is in Hamilton-path order; SCC boundaries are consecutive.
    n=len(t)
    return n==1 or all(not all(edge(t,i,j) for i in range(k) for j in range(k,n)) for k in range(1,n))

def is_prime(t):
    n=len(t)
    for mask in range(1,1<<n):
        if not 1 < mask.bit_count() < n: continue
        inside=[v for v in range(n) if mask>>v&1]
        if all(len({edge(t,v,u) for u in inside})==1 for v in range(n) if not mask>>v&1): return False
    return n>=3

def paths(t):
    n=len(t); full=(1<<n)-1
    def visit(p,used):
        if used==full:
            yield p; return
        allowed=full^used if not p else t[p[-1]] & (full^used)
        while allowed:
            bit=allowed & -allowed; allowed^=bit; v=bit.bit_length()-1
            yield from visit(p+(v,),used|bit)
    yield from visit((),0)

def hamilton(t):
    @lru_cache(None)
    def walk(mask,last):
        if mask==(1<<len(t))-1: return 1
        return sum(walk(mask|(1<<j),j) for j in range(len(t)) if not mask>>j&1 and edge(t,last,j))
    return sum(walk(1<<i,i) for i in range(len(t)))

def automorphisms(t):
    cells=defaultdict(list)
    for v,r in enumerate(t): cells[r.bit_count()].append(v)
    keys=sorted(cells); order=tuple(chain.from_iterable(cells[k] for k in keys))
    out=[]
    for cp in product(*(permutations(cells[k]) for k in keys)):
        p=[0]*len(t)
        for a,b in zip(order,chain.from_iterable(cp)):p[a]=b
        if relabel(t,p)==t:out.append(tuple(p))
    return out

def substitute(q, blocks):
    starts=[0]
    for b in blocks:starts.append(starts[-1]+len(b))
    t=[];owner=[]
    for i,b in enumerate(blocks):
        ext=sum(((1<<len(blocks[j]))-1)<<starts[j] for j in range(len(blocks)) if edge(q,i,j))
        t.extend((r<<starts[i])|ext for r in b);owner.extend([i]*len(b))
    return tuple(t),tuple(owner)
