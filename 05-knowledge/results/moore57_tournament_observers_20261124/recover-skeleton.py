#!/usr/bin/env python3
"""Recover stars and row groups from a Moore grid's known-row Gram residual.

The input to recover is the off-diagonal bit adjacency of a residual H H^T,
with the known diagonal r deleted. It may reject an inconsistent skeleton;
acceptance is not a Moore completion certificate.
"""
def entries(mask):
    while mask:
        low=mask & -mask
        yield low.bit_length()-1
        mask^=low

def pack(ix):
    return sum(1 << i for i in ix)

def recover(R,q,r):
    s=q-r
    assert r>=2 and s>r*r-r+1 and q-2>r*r
    N=s*q
    assert len(R)==N
    used=[0]*N
    stars=[]
    incident=[[] for _ in R]
    for e in range(N):
        while unseen:=R[e]&~used[e]:
            f=next(entries(unseen))
            W=(R[e]&R[f])|(1<<e)|(1<<f)
            star=[v for v in entries(W) if (R[v]&W).bit_count()>=s-1]
            if len(star)!=s or e not in star or f not in star:
                raise ValueError("no large star about residual pair")
            mask=pack(star)
            if any((mask^(1<<v))&~R[v] or used[v]&mask for v in star):
                raise ValueError("star is not a residual clique with fresh edges")
            if len({v//q for v in star})!=s:
                raise ValueError("star misses an old row")
            k=len(stars)
            stars.append(frozenset(star))
            for v in star:
                used[v]|=mask^(1<<v)
                incident[v].append(k)
    if len(stars)!=r*q or any(len(g)!=r for g in incident):
        raise ValueError("wrong recovered incidence degrees")
    Q=r*q
    meet=[0]*Q
    for ids in incident:
        mask=pack(ids)
        for a in ids:
            others=mask^(1<<a)
            if meet[a]&others:
                raise ValueError("recovered hypergraph is not linear")
            meet[a]|=others
    allmask=(1<<Q)-1
    opposite=[allmask^meet[i]^(1<<i) for i in range(Q)]
    seen=set(); groups=[]
    for x in range(Q):
        if x in seen:continue
        row={x}
        row.update(y for y in range(Q) if y!=x and (opposite[x]&opposite[y]).bit_count()>=q-2)
        if len(row)!=q or not row.isdisjoint(seen):
            raise ValueError("same-row test is inconsistent")
        seen.update(row);groups.append(row)
    if len(groups)!=r:
        raise ValueError("wrong number of new rows")
    masks=list(map(pack,groups))
    for g,mask in zip(groups,masks):
        for x in g:
            cand=opposite[x]
            if (cand&mask).bit_count()!=q-1 or any((cand&v).bit_count()!=r for v in masks if v!=mask):
                raise ValueError("incorrect cross degree")
    return stars,groups,incident

def mul8(a,b):
    res=0
    while b:
        if b&1:res^=a
        b>>=1;a<<=1
        if a&8:a^=11 # x^3+x+1
    return res

def eval_product(x,j,t):
    return 8*((x//8+j*(t//8))%7)+((x%8)^mul8(j,t%8))

def test(r):
    q=56;s=q-r
    true_stars=[set() for _ in range(r*q)]
    for t in range(s):
        for x in range(q):
            for j in range(r):
                a=j*q+eval_product(x,j,t)
                true_stars[a].add(t*q+x)
    R=[0]*(s*q)
    for star in true_stars:
        mask=pack(star)
        for v in star:R[v]|=mask^(1<<v)
    stars,rows,ids=recover(R,q,r)
    reference={frozenset(x):i for i,x in enumerate(true_stars)}
    assert set(stars)==set(reference)
    assert {frozenset(reference[stars[i]]//q for i in row) for row in rows}=={frozenset([u]) for u in range(r)}
    if r==6:
        # A seventh slope gives an old-column coordinate unused by the r stars.
        col=[eval_product(x,r,t) for t in range(s) for x in range(q)]
        assert all(len({col[v] for v in star})==s for star in stars)
    print(f"width {q}, {r} missing rows: recovered {len(stars)} stars in {len(rows)} row groups")

if __name__=="__main__":
    test(6);test(7)
