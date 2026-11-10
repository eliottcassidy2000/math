#!/usr/bin/env python3
"""Exhaustive labelled universes n=2..5; standard library, no external data."""
from itertools import combinations, permutations
from collections import Counter
from fractions import Fraction
from math import factorial, lcm
import json
from pathlib import Path


def check(n):
    edges=list(combinations(range(n),2)); index={e:i for i,e in enumerate(edges)}
    M=1<<len(edges); N=1<<(n-1); full=(1<<n)-1
    perms=list(permutations(range(n)))
    def norm(x): return x^full if x>>(n-1)&1 else x
    def px(x,p): return norm(sum(((x>>i)&1)<<p[i] for i in range(n)))
    cuts=[sum((((x>>i)^(x>>j))&1)<<k for k,(i,j) in enumerate(edges)) for x in range(N)]
    def relabel(t,p):
        ans=0
        for k,(i,j) in enumerate(edges):
            a,b=p[i],p[j]
            bit=((t>>k)&1)^(a>b)
            ans|=bit<<index[tuple(sorted((a,b)))]
        return ans
    rel=[[relabel(t,p) for t in range(M)] for p in perms]
    canon=[min(row[t] for row in rel) for t in range(M)]
    def h(t):
        count=0
        for p in perms:
            count+=all((((t>>index[tuple(sorted((a,b)))])&1)==(a<b)) for a,b in zip(p,p[1:]))
        return count
    hs=[h(t) for t in range(M)]
    direct=Counter(hs[t] for t in set(canon))
    remaining=set(range(M)); total=Counter(); classes=[]
    for t in range(M):
        if t not in remaining: continue
        orbit={t^cut:x for x,cut in enumerate(cuts)}
        whole={row[u] for row in rel for u in orbit}
        remaining.difference_update(whole)
        f=[hs[t^cut] for cut in cuts]
        chars=[s for s in range(1<<n) if s.bit_count()%2==0]
        sign=lambda s,x: -1 if (s&x).bit_count()%2 else 1
        hat={s:Fraction(sum(f[x]*sign(s,x) for x in range(N)),N) for s in chars}
        assert all(v==0 for s,v in hat.items() if s.bit_count()>2*(n//3))
        P=[(p,orbit[row[t]]) for p,row in zip(perms,rel) if row[t] in orbit]
        numerator=Counter(); fixed_types=Counter()
        for p,b in P:
            seen=set(); cycles=[]
            for i in range(n):
                if i in seen: continue
                cy=[]; j=i
                while j not in seen: seen.add(j);cy.append(j);j=p[j]
                cycles.append(cy)
            order=lcm(*(len(c) for c in cycles))
            fixed=[x for x in range(N) if b^px(x,p)==x]
            assert all(f[b^px(x,p)]==f[x] for x in range(N))
            if order%2==0:
                assert not fixed
                fixed_types['even_empty']+=1
                continue
            assert len(fixed)==1<<(len(cycles)-1)
            x0=fixed[0]
            allowed=[s for s in chars if all(sum((s>>v)&1 for v in c)%2==0 for c in cycles)]
            avg=sum(hat[s]*sign(s,x0) for s in allowed)
            assert avg==Fraction(sum(f[x] for x in fixed),len(fixed))
            fixed_types['odd_nonempty']+=1
            numerator.update(f[x] for x in fixed)
        assert all(v%len(P)==0 for v in numerator.values())
        poly=Counter({k:v//len(P) for k,v in numerator.items()})
        # Independent ordinary-isomorphism orbit representatives inside this packet.
        reps={canon[u] for u in orbit}
        assert poly==Counter(hs[u] for u in reps)
        total.update(poly)
        # The adjacency descends to affine P-orbits and inherits the boundary spectrum.
        qreps={}
        for x in range(N): qreps.setdefault(canon[t^cuts[x]],x)
        walks=[]; g=list(f)
        for k in range(n//3+4):
            walks.append(sum(g[x] for x in qreps.values()))
            g=[sum(g[norm(x^(1<<v))] for v in range(n)) for x in range(N)]
        coeff=[1]
        for lam in [n-4*j for j in range(n//3+1)]:
            out=[0]*(len(coeff)+1)
            for j,c in enumerate(coeff): out[j]-=lam*c;out[j+1]+=c
            coeff=out
        assert all(sum(c*walks[k+j] for j,c in enumerate(coeff))==0 for k in range(len(walks)-len(coeff)+1))
        assert sum(f)==factorial(n)
        classes.append({'representative':t,'projection_order':len(P),'fixed_types':dict(fixed_types),'H_enumerator':dict(sorted(poly.items())),'unlabelled_walk_H_sums':walks})
    assert total==direct
    assert sum(total.values())=={2:1,3:2,4:4,5:12}[n]
    return {'n':n,'labelled_universe':M,'switching_isomorphism_classes':len(classes),'A000568':sum(total.values()),'H_enumerator':dict(sorted(total.items())),'classes':classes}

if __name__=='__main__':
    result={'status':'FINITE-EXACT','checks':[check(n) for n in range(2,6)]}
    path=Path(__file__).with_name('output.json');path.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
