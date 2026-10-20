#!/usr/bin/env python3
"""Exact constructive checks of the two quotient-black leaf families."""
import json
from functools import lru_cache
from pathlib import Path
from census import schema, tournament, hp_fibre, reflection


@lru_cache(None)
def words(n):
    if n < 0:
        return ()
    if n == 0:
        return ((),)
    return tuple((a,)+w for a in (1,3) for w in words(n-a))


def stacked_mask(blocks):
    """Blocks are (size, local explorer tiling mask), stacked forward."""
    n = sum(s for s,t in blocks)
    tiles,_ = schema(n)
    ix = {(n-x,n-y):k for k,(x,y) in enumerate(tiles)}
    ans = start = 0
    for s,t in blocks:
        for k,(x,y) in enumerate(schema(s)[0]):
            if (t>>k)&1:
                ans |= 1 << ix[start+s-x,start+s-y]
        start += s
    return ans


def blockword(w):
    return [(s,1 if s==3 else 0) for s in w]


def c3(t,n):
    adj=tournament(t,n,schema(n)[0])
    import math
    return math.comb(n,3)-sum(math.comb(a.bit_count(),2) for a in adj)


def main():
    rows=[]
    for n in range(4,13):
        full=(1<<len(schema(n)[0]))-1
        a=len(words(n))
        pal=sum(w==w[::-1] for w in words(n))
        ns=[];mixed=[]
        for w in words(n):
            if w>=w[::-1]:
                continue
            t=stacked_mask(blockword(w))
            f,h,aut=hp_fibre(tournament(t,n,schema(n)[0]),n,schema(n)[0])
            assert f=={t} and h==aut==3**w.count(3)
            assert c3(t,n)==w.count(3)
            assert c3(full^t,n)>=n-2
            ns.append({"word":w,"tiling":t,"c3":c3(t,n),
                       "partner_c3":c3(full^t,n),"H":h,"aut":aut})
        if n%2==0:
            r=(n-4)//2
            for w in words(r):
                expected={stacked_mask(blockword(w)+[(4,s)]+blockword(w[::-1]))
                          for s in (1,3,5,6,7)}
                t=min(expected)
                f,h,aut=hp_fibre(tournament(t,n,schema(n)[0]),n,schema(n)[0])
                assert f==expected and len(f)==5 and h==5*aut
                blue={s for s in f if reflection(s,schema(n)[1])==s}
                assert len(blue)==3
                partner_cs={c3(full^s,n) for s in f-blue}
                assert len(partner_cs)==1
                if n>=6:
                    assert min(partner_cs)>=n
                mixed.append({"word":w,"tilings":sorted(f),"blue":sorted(blue),
                              "c3":c3(t,n),"black_partner_c3":next(iter(partner_cs)),
                              "H":h,"aut":aut})
        assert len(ns)==(a-pal)//2
        if n>=5:
            assert all(3*r["c3"]<=n+2 and r["partner_c3"]>=n-2 for r in ns)
            assert all(3*r["c3"]<=n+2 and r["black_partner_c3"]>=n-2 for r in mixed)
        rows.append({"n":n,"a_n":a,"p_n":pal,"NS_leaves":len(ns),
                     "mixed_leaves":len(mixed),
                     "distinct_bridge_lower_bound":1 if n==4 else len(ns)+len(mixed),
                     "NS_witnesses":ns,"mixed_witnesses":mixed})
        print(n,a,pal,len(ns),len(mixed),rows[-1]["distinct_bridge_lower_bound"],flush=True)
    Path(__file__).with_name("family_checks.json").write_text(json.dumps(rows,indent=2)+"\n")


if __name__=="__main__":
    main()
