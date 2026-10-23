#!/usr/bin/env python3
"""Lossless finite-domain pair-compatibility to bounded integer equations.
States must be distinct hashable objects. Cell domains and edges retain order.
No relaxation, state quotient, implicit mark removal, or inferred native ports.
"""
def compile_pairs(domains, relations):
    domains=tuple(tuple(D) for D in domains)
    if any(len(set(D))!=len(D) for D in domains):raise ValueError('duplicate state')
    rels={}
    for (i,j),R in relations.items():
        if not(0<=i<j<len(domains)):raise ValueError('edge must have i<j')
        R=tuple(R)
        if len(set(R))!=len(R):raise ValueError('duplicate allowed pair')
        if any(a not in domains[i] or b not in domains[j] for a,b in R):raise ValueError('unknown state')
        rels[i,j]=R
    labels=[('z',i,a) for i,D in enumerate(domains) for a in D]
    labels += [('w',i,j,a,b) for (i,j),R in rels.items() for a,b in R]
    idx={v:k for k,v in enumerate(labels)}
    A=[];rhs=[]
    def row(items,b):
        r=[0]*len(labels)
        for label,c in items:r[idx[label]]+=c
        A.append(tuple(r));rhs.append(b)
    for i,D in enumerate(domains):row([(('z',i,a),1) for a in D],1)
    for (i,j),R in rels.items():
        for a in domains[i]:row([(('z',i,a),-1)]+[(('w',i,j,aa,b),1) for aa,b in R if aa==a],0)
        for b in domains[j]:row([(('z',j,b),-1)]+[(('w',i,j,a,bb),1) for a,bb in R if bb==b],0)
    return {'A':tuple(A),'b':tuple(rhs),'l':(0,)*len(labels),'u':(1,)*len(labels),
            'variables':tuple(labels),'domains':domains,'relations':rels}

def encode(model,assignment):
    if len(assignment)!=len(model['domains']):raise ValueError('assignment size')
    return tuple(int(assignment[v[1]]==v[2]) if v[0]=='z' else
                 int(assignment[v[1]]==v[3] and assignment[v[2]]==v[4]) for v in model['variables'])

def feasible(model,x):
    return len(x)==len(model['variables']) and all(type(a) is int and 0<=a<=1 for a in x) and all(sum(c*a for c,a in zip(r,x))==b for r,b in zip(model['A'],model['b']))

def decode(model,x):
    if not feasible(model,x):raise ValueError('not a feasible binary model point')
    chosen={v[1]:v[2] for v,a in zip(model['variables'],x) if v[0]=='z' and a}
    return tuple(chosen[i] for i in range(len(model['domains'])))

def block_frame_relations(A,B,cells,frames):
    """Return exact allowed frame-index pairs after internal references.
    A/B are full 0/1 matrices on identical actual labels. Each frame is a
    tuple giving images of cell vertices in cell order. Frames must already
    preserve the required actual marks; this routine does not infer marks.
    """
    n=len(A)
    if len(B)!=n or any(len(r)!=n for r in A) or any(len(r)!=n for r in B):raise ValueError('matrix dimensions')
    if sorted(v for C in cells for v in C)!=list(range(n)):raise ValueError('not a partition')
    if len(frames)!=len(cells):raise ValueError('frame domain count')
    for C,F in zip(cells,frames):
        if len({tuple(g) for g in F})!=len(F):raise ValueError('duplicate whole-cell frame')
        for g in F:
            if sorted(g)!=sorted(C):raise ValueError('frame not a cell permutation')
            if any(A[u][v]!=B[g[a]][g[b]] for a,u in enumerate(C) for b,v in enumerate(C)):raise ValueError('frame fails internal comparison')
    R={}
    for i,C in enumerate(cells):
        for j in range(i+1,len(cells)):
            D=cells[j]
            R[i,j]=tuple((a,b) for a,g in enumerate(frames[i]) for b,h in enumerate(frames[j])
                        if all(A[u][v]==B[g[k]][h[l]] for k,u in enumerate(C) for l,v in enumerate(D)))
    return tuple(tuple(range(len(F))) for F in frames),R
