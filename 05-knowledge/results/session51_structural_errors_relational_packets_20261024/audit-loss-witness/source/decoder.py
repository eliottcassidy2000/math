"""Exact domain-checked prime cyclic stable-cell automorphism decoder."""
from math import prod, isqrt

def stable_cells(rows):
    n = len(rows)
    if any(r >> n or r & (1 << i) for i, r in enumerate(rows)):
        raise ValueError('invalid loop or bit')
    if any(((rows[i] >> j) & 1) + ((rows[j] >> i) & 1) != 1
           for i in range(n) for j in range(i)):
        raise ValueError('not a tournament')
    colors = [0] * n
    while True:
        cells = [[i for i in range(n) if colors[i] == c]
                 for c in sorted(set(colors))]
        keys = [(colors[i], tuple(sum((rows[i] >> j) & 1 for j in C)
                                 for C in cells)) for i in range(n)]
        palette = {k: j for j, k in enumerate(sorted(set(keys)))}
        new = [palette[k] for k in keys]
        if new == colors:
            return cells
        colors = new

def prime(p):
    return p >= 2 and all(p % d for d in range(2, isqrt(p)+1))

def cyclic_order(rows, cell):
    p = len(cell)
    if not prime(p) or p == 2:
        raise ValueError('cell not odd prime')
    m = (p-1)//2
    successor = {}
    for x in cell:
        choices = [y for y in cell if (rows[x] >> y) & 1
                   and sum(((rows[x] >> z) & 1) * ((rows[y] >> z) & 1)
                           for z in cell) == m-1]
        if len(choices) != 1:
            raise ValueError('cell not consecutive cyclic')
        successor[x] = choices[0]
    order = [min(cell)]
    for _ in range(p-1):
        order.append(successor[order[-1]])
    if len(set(order)) != p or successor[order[-1]] != order[0]:
        raise ValueError('successor not full cycle')
    if any(((rows[order[i]] >> order[j]) & 1) != (0 < (j-i)%p <= m)
           for i in range(p) for j in range(p) if i != j):
        raise ValueError('cyclic template failed')
    return order

def stabilizer(B):
    p = len(B)
    return [(a,b) for a in range(p) for b in range(p)
            if all(B[x][y] == B[(x+a)%p][(y+b)%p]
                   for x in range(p) for y in range(p))]

def decode(rows):
    cells = stable_cells(rows)
    orders = [C if len(C)==1 else cyclic_order(rows,C) for C in cells]
    graph = {i: [] for i,C in enumerate(cells) if len(C)>1}
    pinned = set()
    edges = []
    for i,X in enumerate(orders):
        for j in range(i):
            Y = orders[j]
            B = [[(rows[x] >> y) & 1 for y in Y] for x in X]
            if len(set(v for r in B for v in r)) == 1:
                continue
            if len(X) != len(Y) or len(X)==1:
                raise RuntimeError('coprime equitability failed')
            p = len(X)
            L = stabilizer(B)
            if len(L)==1:
                pinned.update((i,j)); edges.append((i,j,'rigid'))
            elif len(L)==p:
                slope = [b for a,b in L if a==1]
                if len(slope)!=1 or slope[0]==0:
                    raise RuntimeError('invalid line')
                slope=slope[0]
                graph[i].append((j,slope))
                graph[j].append((i,pow(slope,-1,p)))
                edges.append((i,j,slope))
            else:
                raise RuntimeError('invalid nonuniform stabilizer')
    seen=set(); free=[]; frozen=[]; generators=[]
    for root in graph:
        if root in seen: continue
        p=len(cells[root]); values={root:1}; todo=[root]; balanced=True
        while todo:
            i=todo.pop()
            for j,s in graph[i]:
                v=values[i]*s%p
                if j in values:
                    balanced &= values[j]==v
                else:
                    values[j]=v; todo.append(j)
        seen.update(values)
        if balanced and not pinned.intersection(values):
            free.append((p,sorted(values)))
            g=list(range(len(rows)))
            for i,k in values.items():
                C=orders[i]
                for x in range(p): g[C[x]]=C[(x+k)%p]
            generators.append(g)
        else:
            frozen.extend(values)
    orbit_sizes = [1]*len(rows)
    for p, indices in free:
        for i in indices:
            for v in cells[i]: orbit_sizes[v]=p
    error=sum(len(cells[i])*(len(cells[i])-1) for i in frozen)
    return {'cells':cells,'orders':orders,'edges':edges,
            'free_components':free,'frozen_cells':sorted(frozen),
            'group_order':prod(p for p,_ in free),
            'orbit_sizes':orbit_sizes,'error':error,'generators':generators}
