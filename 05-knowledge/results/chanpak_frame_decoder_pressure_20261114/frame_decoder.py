"""Exact whole-count/locked-frame inverse for a fixed two-hub block skeleton.

A segment is (digit, width), where width None denotes a positive variable.
The anchors, if used, are four increasing segment indices of fixed width 1.
The target M, if used, is the contribution of the *one* marked quadruple
across both hub orders and all its seam matchings, not the whole H.
"""
from argparse import ArgumentParser
from itertools import permutations, product
from counter36 import H


def expansion(length, s):
    return (1, (1 << length) - 2, 3**length - 3*(1 << length) + 3)[s-1]


def word(segments, k):
    x = iter(k)
    return ''.join(str(d)*(s if s is not None else next(x)) for d,s in segments)


def coefficient_tensor(segments):
    d = sum(s is None for _,s in segments)
    keys = list(product((1,2,3), repeat=d))
    c = {k:H(word(segments,k)) for k in keys}
    for j in range(d):
        for k in keys:
            if k[j] != 1:
                continue
            k2, k3 = k[:j]+(2,)+k[j+1:], k[:j]+(3,)+k[j+1:]
            a,b,e = c[k], c[k2], c[k3]
            assert (b-a)%2==0 and (e-3*b+2*a)%6==0
            c[k2],c[k3] = (b-a)//2, (e-3*b+2*a)//6
    assert min(c.values()) >= 0
    return c


def evaluate(c, k):
    es = [[expansion(a,s) for s in (1,2,3)] for a in k]
    total = 0
    for t,v in c.items():
        if v:
            for j,s in enumerate(t):
                v *= es[j][s-1]
            total += v
    return total


def outgoing(d, x):
    return (d>>(1-x))&1


def site_grades(segments, anchors):
    i,j,k,l = anchors
    assert 0<=i<j<k<l<len(segments)
    di,dj,dk,dl = [segments[a][0] for a in anchors]
    assert all(segments[a][1]==1 for a in anchors)
    A = sum(outgoing(di,x)*outgoing(dj,1-x) for x in (0,1))
    B = sum((1-outgoing(dk,x))*(1-outgoing(dl,1-x)) for x in (0,1))
    if not (A and B):
        return None
    grade = [1 if i<t<j or k<t<l else 2 if j<t<k else 0
             for t in range(len(segments))]
    fixed = tuple(sum(s for (_,s),g in zip(segments,grade)
                      if s is not None and g==m) for m in (1,2))
    variable = [g for (_,s),g in zip(segments,grade) if s is None]
    return 2*A*B, fixed, variable


def valuation(n,p):
    c=0
    while n%p==0:
        n//=p
        c+=1
    return c,n


def decode(N, segments, anchors=None, M=None):
    """Return all active length tuples, choosing width one at neutral segments."""
    if not N>0 or (M is not None and not 0<M<=N):
        return []
    if M is not None and anchors is None:
        raise ValueError('a marked target needs the four anchor segments')
    d=sum(s is None for _,s in segments)
    c=coefficient_tensor(segments)
    live=[j for j in range(d) if c[(1,)*j+(2,)+(1,)*(d-j-1)]>0]
    # The theorem shows c_2>0 at minimum other widths for each active variable.
    for j in range(d):
        if j not in live:
            assert all(k[j]==1 or not a for k,a in c.items())
    wanted=None
    fixed=(0,0)
    vg=[0]*d
    if anchors is not None:
        g=site_grades(segments,anchors)
        if g is None:
            return []
        mult,fixed,vg=g
        if M is not None:
            if M%mult:
                return []
            a,t=valuation(M//mult,2)
            b,t=valuation(t,3)
            if t!=1 or a<fixed[0] or b<fixed[1]:
                return []
            wanted=(a,b)
    if not live:
        solutions=[]
        k=(1,)*d
        g=tuple(fixed[m-1]+sum(k[j] for j in range(d) if vg[j]==m)
                for m in (1,2))
        return [k] if evaluate(c,k)==N and (wanted is None or g==wanted) else solutions
    L=(N+2).bit_length()-1
    # Use one linear grade equation per nonempty constrained group to eliminate
    # a variable, and use monotonicity to invert one remaining variable.
    solve={}
    if wanted is not None:
        for g in (1,2):
            batch=[j for j in live if vg[j]==g]
            if batch:
                solve[g]=batch[-1]
    free=[j for j in live if j not in solve.values()]
    root=free.pop() if free else None
    results=[]
    for prefix in product(range(1,L+1),repeat=len(free)):
        k=[1]*d
        for j,p in zip(free,prefix):
            k[j]=p
        if wanted is not None:
            good=True
            for g in (1,2):
                v=solve.get(g)
                used=fixed[g-1]+sum(k[j] for j in range(d) if vg[j]==g and j!=v)
                if v is not None:
                    k[v]=wanted[g-1]-used
                    good &= 1<=k[v]<=L
                else:
                    good &= used==wanted[g-1]
            if not good:
                continue
        if root is not None:
            # All grade-assigned variables were solved above; a free root with
            # prescribed grades cannot occur in either constrained group.
            assert wanted is None or not vg[root]
            lo,hi=1,L
            while lo<hi:
                mid=(lo+hi)//2
                k[root]=mid
                if evaluate(c,k)<N:
                    lo=mid+1
                else:
                    hi=mid
            k[root]=lo
        if evaluate(c,k)==N:
            results.append(tuple(k))
    return results


def adjacency(w):
    m=len(w)
    out=[sum(1<<j for j in range(i+1,m)) for i in range(m)]+[1<<(m+1),0]
    for i,d in enumerate(map(int,w)):
        for h in (0,1):
            if outgoing(d,h):
                out[m+h]|=1<<i
            else:
                out[i]|=1<<(m+h)
    return out


def adjacency_invariant(out):
    n=len(out)
    blocks=[[v for v,a in enumerate(out) if a.bit_count()==d] for d in range(n)]
    def code(order):
        for i,a in enumerate(order):
            for b in order[i+1:]:
                yield bool(out[a]&(1<<b))
    return n, min(tuple(code(sum(p,()))) for p in product(*[tuple(permutations(B)) for B in blocks if B]))


def choose_class(N, candidates):
    """Input each item as (segments, anchor indices or None, M or None)."""
    words={word(S,k) for S,A,M in candidates for k in decode(N,S,A,M)}
    return min((adjacency_invariant(adjacency(w)),w) for w in words)[1] if words else None


if __name__=='__main__':
    p=ArgumentParser(description=__doc__)
    p.add_argument('target',type=int)
    p.add_argument('pattern',help='space-separated segments: 1 0* 3 0 2* 2 1')
    p.add_argument('--anchors',type=int,nargs=4)
    p.add_argument('--marked',type=int)
    z=p.parse_args()
    S=[(int(t[0]),None if t.endswith('*') else int(t[1:] or 1))
       for t in z.pattern.split()]
    for k in decode(z.target,S,z.anchors,z.marked):
        print(k,word(S,k))
