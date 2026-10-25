"""Independent affine proof audit, using bit-mask Gaussian elimination."""
import itertools
import json
import math


def need(ok, msg):
    if not ok:
        raise RuntimeError(msg)


def add(basis, row):
    while row:
        i = row.bit_length()-1
        if i not in basis:
            basis[i] = row
            return
        row ^= basis[i]


def reduced(basis, row):
    while row:
        i = row.bit_length()-1
        if i not in basis:
            return row
        row ^= basis[i]
    return 0


def setup(n, h, incoming=None):
    pairs = list(itertools.combinations(range(n), 2))
    # Constant has LOW bit, so elimination of variables leaves 1 on inconsistency.
    ix = {p: 1 << (i+1) for i, p in enumerate(pairs)}
    def bit(u,v):
        return ix[tuple(sorted((u,v)))]
    def force(basis,u,v,val=1):
        add(basis, bit(u,v) | (val ^ (u>v)))
    def blue(basis,q,flip):
        r = dict(zip(q, reversed(q)))
        for u,v in pairs:
            a,b=r[u],r[v]
            e,f=bit(u,v),bit(a,b)
            rhs=1 ^ (a>b) ^ bool(flip&e) ^ bool(flip&f)
            add(basis,e^f^rhs)
    p=tuple(range(n))
    basis={}
    for u,v in zip(p,p[1:]):
        force(basis,u,v)
    ports={h,n-3,n-2,n-1} if incoming is None else set(incoming)
    for v in range(1,n):
        force(basis,0,v,int(v not in ports))
    blue(basis,p,0)
    return pairs,ix,bit,force,blue,basis


def check_ledger(n,h,words,guard_arcs,incoming=None):
    pairs,ix,bit,force,blue,base=setup(n,h,incoming)
    for arc in guard_arcs:
        force(base,*arc)
    if 0 in base:
        return False
    full=sum(ix.values())
    flip=0
    for q in words:
        need(sorted(q)==list(range(n)), 'word does not visit every vertex')
        edges=0
        for u,v in zip(q,q[1:]):
            e=bit(u,v)
            req=e ^ (1 ^ (u>v) ^ bool(flip&e))
            need(reduced(base,req)==0, ('join not universally legal',n,h,q,u,v))
            edges |= e
        trial=dict(base)
        blue(trial,q,flip)
        need(0 in trial, ('blue word still feasible',n,h,q))
        flip ^= full ^ edges
    target=full ^ sum(bit(i,i+1) for i in range(n-1))
    need(flip==target, ('literal target differs',n,h))
    return True


def all_reversals(vertices):
    if not vertices:
        yield {}
        return
    a=vertices[0]
    if len(vertices)%2:
        for rest in all_reversals(vertices[1:]):
            yield {a:a,**rest}
    for j in range(1,len(vertices)):
        b=vertices[j]
        for rest in all_reversals(vertices[1:j]+vertices[j+1:]):
            # Remaining odd block must have a fixed point; already handled recursively.
            yield {a:b,b:a,**rest}


def main():
    cyclic=0
    for n in range(11,46):
        for h in range(2,n-3):
            w=tuple(range(1,h)); tail=tuple(range(h,n-3))
            q1=(n-1,0)+w+tail+(n-3,n-2)
            q2=(n-1,0)+tail+w+(n-3,n-2)
            q3=(0,)+tail+w+(n-3,n-2,n-1)
            cyclic += check_ledger(n,h,[q1,q2,q3],[(n-2,1),(3,n-2),(n-h,2)])
    overlap=0
    for n in range(11,41):
        for ell in range(1,(n-3)//2+1,2):
            d=tuple(range(3,ell+3)); u=tuple(range(ell+3,n-1))
            words=[(0,2)+d+u+(n-1,1), d+(0,2)+u+(n-1,1),d+(0,1,2)+u+(n-1,)]
            for h in range(2,n-3):
                if h in {2,ell+2} or (ell>=3 and h in {4,ell-1}):
                    continue
                overlap += check_ledger(n,h,words,[(ell+3,2)])
    involutions=0
    for n in range(2,10):
        pairs,ix,bit,force,blue,_=setup(n,2)
        r=dict(enumerate(reversed(range(n))))
        # Enumerate actual matching involutions, not arbitrary permutations.
        for s in all_reversals(tuple(range(n))):
            if sum(i==s[i] for i in s)>1:
                continue
            b={}
            for rev in (r,s):
                for u,v in pairs:
                    a,c=rev[u],rev[v]
                    add(b,bit(u,v)^bit(a,c)^(1^(a>c)))
            g={i:r[s[i]] for i in range(n)}
            seen=set(); lengths=[]
            for i in range(n):
                if i in seen:
                    continue
                j=i; k=0
                while j not in seen:
                    seen.add(j);k+=1;j=g[j]
                lengths.append(k)
            odd=all(k%2 for k in lengths)
            need((0 not in b)==odd, ('reflection criterion',n,s,lengths))
            involutions+=1
    tail=0
    for n in range(11,61):
        h=n-4
        d=(3,); u=tuple(range(4,n-1))
        old=[(0,2)+d+u+(n-1,1),d+(0,2)+u+(n-1,1),d+(0,1,2)+u+(n-1,)]
        tail+=check_ledger(n,h,old,[(4,2)])
        if n>=12:
            u=tuple(range(4,n-4))
            new=[(0,1,2)+u+(n-4,n-3,n-2,n-1,3),
                 (1,2)+u+(0,n-4,n-3,n-2,n-1,3),
                 (1,2,3)+u+(0,n-4,n-3,n-2,n-1)]
            tail+=check_ledger(n,h,new,[(2,4)])
    # Literal, independent order-eleven exceptional detour.
    source=(126,220,472,976,1952,1102,1936,1825,1571,1063,15)
    n=11
    arcs={(i,j) for i in range(n) for j in range(n) if source[i]>>j&1}
    need(all(((i,j) in arcs)^((j,i) in arcs) for i,j in itertools.combinations(range(n),2)), 'exception not tournament')
    original=set(arcs)
    words=[(0,1,2,3,4,7,5,6,8,9,10),(0,1,2,3,4,7,5,9,10,6,8),(0,1,2,3,4,5,9,10,6,7,8)]
    for q in words:
        edges=set(zip(q,q[1:])); need(edges<=arcs,'exception illegal join')
        r=dict(zip(q,reversed(q)))
        need(any((r[v],r[u]) not in arcs for u,v in arcs),'exception blue move')
        arcs=edges|{(v,u) for u,v in arcs-edges}
    path=set(zip(range(10),range(1,11)))
    need(arcs==path|{(v,u) for u,v in original-path}, 'exception wrong target')
    general_tail=0
    for k in range(4,13):
        for n in range(2*k+4,2*k+15):
            e=n-k; incoming=range(e,n)
            d=(3,); v=tuple(range(4,n-1))
            old=[(0,2)+d+v+(n-1,1),d+(0,2)+v+(n-1,1),d+(0,1,2)+v+(n-1,)]
            general_tail+=check_ledger(n,e,old,[(4,2)],incoming)
            u=tuple(range(4,e)); end=tuple(range(e,n))
            new=[(0,1,2)+u+end+(3,), (1,2)+u+(0,)+end+(3,), (1,2,3)+u+(0,)+end]
            general_tail+=check_ledger(n,e,new,[(2,4)],incoming)
    print(json.dumps({'status':'PASS','cyclic_prefix_nonempty_affine_cells':cyclic,
        'cyclic_range':[11,45],'odd_prefix_nonempty_affine_cells':overlap,
        'odd_prefix_range':[11,40],'two_reflection_checks':involutions,
        'two_reflection_range':[2,9],
        'terminal_port_affine_cells':tail,'terminal_port_range':[11,60],
        'order_eleven_exception_literal_checks':3,
        'general_tail_affine_cells':general_tail,'general_tail_k_range':[4,12],
        'general_tail_n_range':'2k+4 through 2k+14',
        'scope':'all completions in each stated affine cell; independent bit-mask elimination'},indent=2))


if __name__=='__main__':
    main()
