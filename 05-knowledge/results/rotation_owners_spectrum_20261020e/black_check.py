"""Exact labelled replays for all-order rotation-symmetric endcap bypasses.

Research continuation: October 20, 2026. No third-party dependencies.
Explicit checks survive python -O; no output files are written by this code.
"""
from itertools import combinations, permutations


def need(value, message):
    if not value:
        raise RuntimeError(message)


def arc(t, u, v):
    return (t[u] >> v) & 1


def circulant(r, steps):
    need(r % 2 == 1, "odd circulant order")
    need(all((d in steps) != (r-d in steps) for d in range(1,r)), "skew step set")
    return tuple(sum(1 << v for v in range(r) if (v-u) % r in steps) for u in range(r))


def step_sets(r, require_one=True):
    start = 2 if require_one else 1
    for bits in range(1 << (r//2-start+1)):
        yield ({1} if require_one else set()) | {
            d if bits >> (d-start) & 1 else r-d for d in range(start,r//2+1)}


def transpose(t):
    n=len(t)
    return tuple(sum(arc(t,n-1-v,n-1-u) << v for v in range(n)) for u in range(n))


def join(parts):
    n=sum(map(len,parts)); offset=0; rows=[]
    for t in parts:
        later=sum(1 << v for v in range(offset+len(t),n))
        rows.extend((row << offset) | later for row in t)
        offset+=len(t)
    return tuple(rows)


def fixed_path(n, code):
    rows=[1 << (u+1) if u+1<n else 0 for u in range(n)]
    free=[(u,v) for u,v in combinations(range(n),2) if v>u+1]
    for bit,(u,v) in enumerate(free):
        a,b=(u,v) if code >> bit & 1 else (v,u)
        rows[a] |= 1 << b
    return tuple(rows)


def is_path(t,p):
    return sorted(p)==list(range(len(t))) and all(arc(t,u,v) for u,v in zip(p,p[1:]))


def blue(t,p):
    n=len(t)
    return all(arc(t,p[i],p[j])==arc(t,p[n-1-j],p[n-1-i])
               for i in range(n) for j in range(i+1,n))


def reverse_path(t,p):
    need(is_path(t,p), "invalid directed Hamiltonian path")
    rows=list(t)
    for u,v in zip(p,p[1:]):
        rows[u] ^= 1 << v
        rows[v] ^= 1 << u
    return tuple(rows)


def edge_set(p):
    return {tuple(sorted(e)) for e in zip(p,p[1:])}


def automorphism(t,q):
    return all(arc(t,u,v)==arc(t,q[u],q[v]) for u,v in combinations(range(len(t)),2))


def cut(steps,r):
    return next(((t,k) for t in sorted(steps) if t>=2
                 for k in range(t,r-1) if k not in steps),None)


def words(prefix, right, steps):
    """Return three actual path words and the branch/cut used."""
    r=len(right); chosen=cut(steps,r); A=tuple(prefix)
    if chosen:
        t,k=chosen
        ell=t-1; u=k-ell+1
        a=right[0]; C=tuple(right[1:u]); K=tuple(right[u:k+1]); D=tuple(right[k+1:])
        need(C and K and D, "extraction blocks must be nonempty")
        ps=(A+K+(a,)+C+D,
            (a,)+K[::-1]+A[::-1]+D[::-1]+C[::-1],
            C+K+D+A+(a,))
        return ps,("extraction",t,k)
    threshold={1}|set(range((r+1)//2,r-1))
    need(steps==threshold, "uncut step set must equal the threshold pattern")
    if r==5:
        x=right
        ps=(A+(x[0],x[1],x[4],x[2],x[3]),
            (x[1],x[2],x[0])+A[::-1]+(x[3],x[4]),
            (x[3],)+A+(x[0],x[2],x[4],x[1]))
        return ps,("inherited-five",)
    need(r>=7, "threshold repair needs nonempty F")
    m=r//2; a,b=right[:2]; z=right[-1]
    F=tuple(right[2:m]); G=tuple(right[m:-1])
    ps=(A+(a,b)+F+(z,)+G,
        G[::-1]+(z,)+F[::-1]+(a,)+A[::-1]+(b,),
        (b,)+A+(a,)+F+G+(z,))
    return ps,("threshold",)


def replay(r,steps,middle):
    b=circulant(r,steps)
    source=join((transpose(b),middle,b)); n=len(source); boundary=n-r
    ps,branch=words(tuple(range(boundary)),tuple(range(boundary,n)),steps)
    target=reverse_path(source,tuple(range(n)))
    need(edge_set(ps[0]) ^ edge_set(ps[1]) ^ edge_set(ps[2]) == edge_set(tuple(range(n))),
         "path-edge parity does not match the original target")
    need(blue(source,tuple(range(n))) == blue(middle,tuple(range(len(middle)))),
         "middle blue equivalence")
    state=source; endpoint_sums=[]
    for p in ps:
        need(is_path(state,p), "lifted word is not a directed path")
        need(not blue(state,p), "detour contains a blue step")
        endpoint_sums.append(state[p[0]].bit_count()+state[p[-1]].bit_count())
        state=reverse_path(state,p)
    need(state==target, "labelled target mismatch")
    if branch[0]=="extraction":
        need(endpoint_sums[1:]==[r-1,r-1], "regular extraction endpoint charges")
    elif branch[0]=="threshold":
        need(endpoint_sums[1:]==[r,r], "threshold endpoint charges")
    return branch


def main():
    print("Exact replay of PROVED all-order circulant-endcap macros; no global connectivity claim.")
    total=0
    print("All marked circulants (order, step sets containing1, extraction, threshold/base):")
    for r in range(5,24,2):
        branches={"extraction":0,"threshold":0,"inherited-five":0}
        for steps in step_sets(r):
            branch=replay(r,steps,())
            branches[branch[0]]+=1;total+=1
        print(r,sum(branches.values()),branches["extraction"],branches["threshold"]+branches["inherited-five"])
    print("Total marked empty-middle circulant presentations:",total)

    # Independent exhaustive small test of the proper-interval automorphism lemma.
    tests=0
    for r in range(3,12,2):
        for steps in step_sets(r,False):
            b=circulant(r,steps)
            for L in range(2,r):
                for shift in range(1,L):
                    q=tuple(range(L-shift,L))+tuple(range(L-shift))+tuple(range(L,r))
                    need(not automorphism(b,q), "proper-interval rotation automorphism")
                    tests+=1
    print("All small skew step sets / proper-interval rotations checked:",tests)

    # Every fixed-path middle, with no blue filter, in three distinct mechanisms.
    banks=((7,{1,2,3}), (11,{1,3,4,5,9}), (15,{1}|set(range(8,14))))
    middle_count=0; replay_count=0
    for size in range(7):
        free=(size-1)*(size-2)//2 if size>=2 else 0
        for code in range(1 << free):
            middle=fixed_path(size,code);middle_count+=1
            for r,steps in banks:
                replay(r,steps,middle);replay_count+=1
    print("Arbitrary middle presentations through order6:",middle_count,"replays:",replay_count)
    print("Paley11: no unit d has d,2d in S; retained block cut",cut(banks[1][1],11))
    need(not any(d in banks[1][1] and 2*d%11 in banks[1][1] for d in range(1,11)),
         "Paley11 must obstruct the single-step-d criterion")

    # Long threshold words test the variable-length F/G ledger, not a new theorem by extrapolation.
    for r in range(7,102,2):
        replay(r,{1}|set(range((r+1)//2,r-1)),fixed_path(3,0))
    print("Long threshold controls: all48 odd orders7..101 with a three-vertex middle.")

    c3=circulant(3,{1}); low=join((c3,c3)); count=0
    for p in permutations(range(6)):
        if is_path(low,p):
            need(blue(low,p), "order-three endcap hostile lost pure-blue property")
            count+=1
    need(count==9,"C3->C3 path inventory")
    print("Lower-bound hostile C3->C3: all9 Hamiltonian paths are blue.")
    print("PASS exact directions, blackness, XOR ledgers, labels, symmetry guards and hostile controls.")


if __name__=="__main__":
    main()
