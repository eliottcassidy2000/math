#!/usr/bin/env python3
"""Exact, standard-library checks for the signed Schubert/threshold packet.
No repository imports. All checks remain active under python -O.
"""
import itertools, json, math
from functools import lru_cache
from pathlib import Path


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def mul(a, b):
    out = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i+j] += x*y
    return out


def qint(k, sign=1, step=1):
    out = [0] * (step*(k-1)+1)
    for i in range(k):
        out[step*i] = sign**i
    return out


def qfactorial(n):
    out = [1]
    for k in range(1, n+1):
        out = mul(out, qint(k))
    return out


def edges_from_mask(n, mask):
    return frozenset(e for i, e in enumerate(itertools.combinations(range(n), 2))
                     if mask >> i & 1)


def direct(n, edges):
    """Ranks p(v), matching epsilon(p,X) in Session58."""
    out = [0] * (n*(n-1)//2+1)
    for p in itertools.permutations(range(n)):
        length = sum(p[i] > p[j] for i in range(n) for j in range(i+1,n))
        parity = sum(p[i] > p[j] for i,j in edges) % 2
        out[length] += (-1)**parity
    return out


def subset_dp(n, edges):
    """Independent last-vertex recurrence on permutation words p^{-1}."""
    upper_edges = [sum(1 << u for u in range(v+1,n) if (v,u) in edges)
                   for v in range(n)]
    @lru_cache(None)
    def calc(s):
        if not s:
            return (1,)
        size=s.bit_count()
        out=[0]*(size*(size-1)//2+1)
        for v in range(n):
            if s >> v & 1:
                rest=s^(1 << v)
                shift=(rest >> (v+1)).bit_count()
                sign=(-1)**((rest & upper_edges[v]).bit_count())
                for k,c in enumerate(calc(rest)):
                    out[k+shift] += sign*c
        return tuple(out)
    return list(calc((1 << n)-1))


def taylor_at_one(poly):
    return [sum(c*math.comb(i,k) for i,c in enumerate(poly) if i>=k)
            for k in range(len(poly))]


def valuation_at_one(poly):
    return next(i for i,c in enumerate(taylor_at_one(poly)) if c)


def threshold(r,j):
    """Construction I^r U I U ... I U; every new vertex has largest label."""
    n=r
    edges=set()
    for t in range(j):
        if t:
            n+=1  # new isolated vertex
        edges.update((i,n) for i in range(n))
        n+=1     # new universal vertex
    return n,frozenset(edges)


def threshold_formula(r,j):
    out=qfactorial(r)
    for t in range(1,j+1):
        out=mul(out,qint(r+2*t-1,-1))
    for t in range(1,j):
        out=mul(out,qint(r+2*t))
    return out


def effective_remainder(r,j):
    require(r%2==1,'effective formula uses odd r')
    out=qfactorial(r)
    for t in range(1,j+1):
        out=mul(out,qint((r+2*t-1)//2,step=2))
    for t in range(1,j):
        out=mul(out,qint(r+2*t))
    return out


def royle_even(n,edges):
    count=0
    for p in itertools.permutations(range(n)):
        mapped=frozenset(tuple(sorted((p[i],p[j]))) for i,j in edges)
        if mapped==edges:
            count+=1
            if sum(p[i]>p[j] for i,j in edges)%2:
                return False,count
    return True,count


def evaluate(poly,q):
    return sum(c*q**i for i,c in enumerate(poly))


def profile_q(n,edges,base_order):
    # Here a tournament's backward-edge set in the standard order is `edges`.
    pairs=list(itertools.combinations(range(n),2))
    permutations=list(itertools.permutations(range(n)))
    def graph_key(mask):
        graph=edges_from_mask(n,mask)
        codes=[]
        for p in permutations:
            rel=frozenset(tuple(sorted((p[i],p[j]))) for i,j in graph)
            code=sum(1<<k for k,e in enumerate(pairs) if e in rel)
            codes.extend([code,((1<<len(pairs))-1)^code])
        return min(codes)
    ranks={v:i for i,v in enumerate(base_order)}
    out={}
    for order in permutations:
        pos={v:i for i,v in enumerate(order)}
        distance=sum((pos[u]-pos[v])*(ranks[u]-ranks[v])<0 for u,v in pairs)
        mask=sum(1<<k for k,(u,v) in enumerate(pairs)
                 if (((u,v) in edges) != (pos[u]>pos[v])))
        key=graph_key(mask)
        out.setdefault(key,[0]*(len(pairs)+1))[distance]+=1
    return out


def main():
    tested_graphs=0
    for n in range(1,6):
        total=qfactorial(n)
        for mask in range(1<<(n*(n-1)//2)):
            edges=edges_from_mask(n,mask)
            dp=subset_dp(n,edges)
            require(dp==direct(n,edges),f'direct/DP mismatch {n},{mask}')
            require(dp[0]==1,'constant coefficient')
            require(dp[::-1]==[(-1)**len(edges)*x for x in dp],'reciprocity')
            require(all((a+b)%2==0 and (a+b)>=0 and (a-b)>=0
                        for a,b in zip(total,dp)),'positive parity cells')
            if len(edges)%2==0 and sum(dp)==0:
                require(valuation_at_one(dp)%2==0,'even-edges vanishing parity')
            tested_graphs+=1
    family=[]
    for r in range(2,8):
        for j in range(1,5):
            n,edges=threshold(r,j)
            if n>12:
                continue
            poly=subset_dp(n,edges)
            require(poly==threshold_formula(r,j),'threshold factorization')
            require(len(edges)==j*(r+j-1),'threshold edge count')
            if n<=7:
                even,aut=royle_even(n,edges)
                require(even and aut==math.factorial(r),'threshold automorphisms')
            if r%2:
                tay=taylor_at_one(poly)
                require(valuation_at_one(poly)==j,'threshold exact zero order')
                require(tay[j]*2**j==(-1)**j*math.factorial(n),'leading jet')
                rem=effective_remainder(r,j)
                factor=[(-1)**k*math.comb(j,k) for k in range(j+1)]
                require(poly==mul(factor,rem),'effective factorization')
                require(all(c>=0 for c in rem),'positive remainder')
                if j%2==0:
                    torus_coeff=taylor_at_one(poly)
                    require(all(c>=0 for c in torus_coeff),'effective torification')
            else:
                denom=math.prod(r+2*t-1 for t in range(1,j+1))
                require(sum(poly)*denom==math.factorial(n),'even r positive control')
            family.append({'r':r,'j':j,'n':n,'edges':len(edges),
                           'valuation':valuation_at_one(poly),'sum':sum(poly)})
    # Deep zero: j=6 is invisible through the fifth derivative, with n=14.
    n,edges=threshold(3,6)
    deep=subset_dp(n,edges)
    require(deep==threshold_formula(3,6),'deep factorization')
    require(valuation_at_one(deep)==6,'deep jet order')
    family.append({'r':3,'j':6,'n':n,'edges':len(edges),'valuation':6,
                   'sum':sum(deep),'leading_jet':taylor_at_one(deep)[6]})
    # Original labelled F_r presentation, distinct from threshold construction order.
    f_records=[]
    for r in range(1,9):
        n=r+3
        edges=frozenset([(0,1),(0,n-1)]+
                        [(v,z) for z in range(2,r+2) for v in (0,1)])
        poly=subset_dp(n,edges)
        require(poly==threshold_formula(r,2),'original F_r gauge factorization')
        if n<=6:
            require(poly==direct(n,edges),'original F_r direct')
        if r%2:
            require(taylor_at_one(poly)[2]*4==math.factorial(n),'F_r quadratic jet')
        f_records.append({'r':r,'n':n,'valuation':valuation_at_one(poly),
                          'coefficients':poly})
    invisible=[]
    for mask in [255,2015]:
        edges=edges_from_mask(6,mask)
        poly=subset_dp(6,edges)
        require(poly==direct(6,edges),'n6 direct')
        require(royle_even(6,edges)[0],'n6 Royle-even')
        require(valuation_at_one(poly)==2,'n6 exact second order')
        invisible.append({'mask':mask,'coefficients':poly,
                          'second_derivative_at_one':2*taylor_at_one(poly)[2]})
    require(invisible[0]['second_derivative_at_one']==360,'mask255 derivative')
    require(invisible[1]['second_derivative_at_one']==408,'mask2015 derivative')
    # Smallest even-edge Royle-even gauge hostile.
    star_a=direct(3,frozenset([(0,1),(0,2)]))
    star_b=direct(3,frozenset([(0,1),(1,2)]))
    require(star_a==[1,0,0,1] and star_b==[1,-2,-2,1],'gauge hostile')
    require(star_a!=star_b and star_a!=[-c for c in star_b],'not even sign-gauge covariance')
    # Effectivity and vanishing depth themselves depend on the label gauge.
    p=[0,2,1,3,4,5]
    relabelled=frozenset(tuple(sorted((p[i],p[j]))) for i,j in edges_from_mask(6,255))
    require(relabelled==edges_from_mask(6,1599),'effectivity hostile relabelling')
    hostile_poly=subset_dp(6,relabelled)
    require(hostile_poly==direct(6,relabelled),'effectivity hostile direct')
    require(evaluate(hostile_poly,2)==-1953,'negative finite-field-count hostile')
    require(taylor_at_one(hostile_poly)[2]==-60,'negative minimum-rank coefficient')
    n,deep_edges=threshold(3,4)
    p=list(range(n));p[0],p[8]=p[8],p[0]
    relabelled_deep=frozenset(tuple(sorted((p[i],p[j]))) for i,j in deep_edges)
    changed_poly=subset_dp(n,relabelled_deep)
    require(valuation_at_one(changed_poly)==2,'depth drops under relabelling')
    require(taylor_at_one(changed_poly)[2]==5760,'depth-drop coefficient')
    # Reference-order symmetrization, literal graph/complement observer at n=3,4.
    observer_records=[]
    for n in [3,4]:
        ps=list(itertools.permutations(range(n)))
        for mask in range(1<<(n*(n-1)//2)):
            edges=edges_from_mask(n,mask)
            profiles=[profile_q(n,edges,b) for b in ps]
            keys=set().union(*(p.keys() for p in profiles))
            uniform={k:sum(profiles[0].get(k,[])) for k in keys}
            for k in keys:
                summed=[sum(p.get(k,[0]*(n*(n-1)//2+1))[d] for p in profiles)
                        for d in range(n*(n-1)//2+1)]
                require(summed==[uniform[k]*x for x in qfactorial(n)],'reference averaging')
            if mask==0:
                observer_records.append({'n':n,'transitive_base_profiles_differ':
                                         any(p!=profiles[0] for p in profiles[1:])})
    # Negative control: nonzero polynomial does not imply nonzero Euler value.
    require(evaluate(invisible[0]['coefficients'],2)>0 and
            evaluate(invisible[0]['coefficients'],1)==0,'Euler hostile')
    result={'status':'PASS','all_graphs_direct_dp_through_n':5,
            'all_graph_count':tested_graphs,'threshold_cases':family,
            'original_F_r':f_records,'n6_invisible':invisible,
            'gauge_hostile':{'star_center0':star_a,'star_center1':star_b,
                             'F3_relabelled_mask1599':hostile_poly,
                             'mask1599_at_q2':-1953,
                             'H_3_4_after_transposition_0_8':{'valuation':2,'jet2':5760}},
            'reference_order_symmetrization':observer_records,
            'checks_use_assert':False}
    Path(__file__).with_name('certificate.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'status':'PASS','graphs_checked':tested_graphs,
                      'threshold_cases':len(family),'original_F_r_cases':len(f_records),
                      'n6_second_derivatives':[x['second_derivative_at_one'] for x in invisible],
                      'deep_zero_order':6},sort_keys=True))

if __name__=='__main__':
    main()
