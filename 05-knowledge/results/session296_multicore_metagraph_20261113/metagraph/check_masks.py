from itertools import combinations, product
from math import comb
import json

def require(condition):
    if not condition:
        raise RuntimeError("mask law failed")

cases = {"cards": 0, "joint_shortcut_pairs": 0, "central_or_off_central_deletions": 0}
for n in range(3,7):
    ix = list(range(n))
    E = list((i,j) for i,j in combinations(ix,2) if j-i>1)
    rev = [E.index((n-1-j,n-1-i)) for i,j in E]
    orbits = {min(i,j) for i,j in enumerate(rev)}
    data = list(product(range(2),repeat=len(E)))
    blue = [v for v in data if all(v[i] == v[j] for i,j in enumerate(rev))]
    masks = []
    for k in range(1,n+1):
        for S in combinations(ix,k):
            P = [j for j,(a,b) in enumerate(E) if a in S and b in S]
            D = [E.index((a,b)) for a,b in zip(S,S[1:]) if b-a>1]
            if D: masks.append(set(D))
            for bits in product(range(2),repeat=len(P)-len(D)):
                p = dict(zip([i for i in P if i not in D],bits))
                p.update({i:1 for i in D})
                def matches(v): return all(v[i] == t for i,t in p.items())
                hits = [v for v in data if matches(v)]
                B = sum(matches(v) for v in blue)
                require(len(hits) == 2**(len(E)-len(P)))
                p_ref = {min(i,rev[i]) for i in p}
                clash = any(p[i] != p[rev[i]] for i in p if rev[i] in p)
                require(B == (0 if clash else 2**(len(orbits)-len(p_ref))))
                if D: require(not any(matches(tuple(1-x for x in v)) for v in hits))
                if k==n-1:
                    x = next(iter(set(ix)-set(S)))
                    y = n-1-x
                    exponent = (n-3)//2 if x==y else int(abs(x-y)>1)
                    require(B == (0 if clash else 2**exponent))
                    cases["central_or_off_central_deletions"]+=1
                cases["cards"]+=1
    if n >= 4:
        cases.setdefault("adjacent_internal_line_pairs",0)
        idx=E.index((0,n-1))
        lines=[v for v in data if v[idx]==0]
        def good_card(v,S):
            bits=tuple(v[E.index((a,b))] for a,b in combinations(S,2)
                if S.index(b)-S.index(a)>1)
            D=[E.index((a,b)) for a,b in zip(S,S[1:]) if b-a>1]
            if D and v[D[0]]==0: return tuple(1-t for t in bits)
            return bits
        for j in range(1,n-2):
            J,K=[i for i in ix if i!=j],[i for i in ix if i!=j+1]
            pair=[(good_card(v,J),good_card(v,K)) for v in lines]
            require(len(set(pair))==len(lines))
            if n%2==0 and j==n//2-1:
                E2=list((i,l) for i,l in combinations(range(n-1),2) if l-i>1)
                r2=[E2.index((n-2-l,n-2-i)) for i,l in E2]
                require(all((a==tuple(b[t] for t in r2)) == (v in blue)
                    for (a,b),v in zip(pair,lines)))
            cases["adjacent_internal_line_pairs"]+=len(lines)
        if n>=5:
            def forget(a): return min(a,tuple(1-t for t in a))
            T=[0,1,n-2,n-1]
            deck=[tuple(forget(good_card(v,S)) for S in (ix[:-1],ix[1:],T)) for v in lines]
            require(len(set(deck))==len(lines))
            cases["three_phase_only_cards"]=cases.get("three_phase_only_cards",0)+len(lines)
    masks = list(map(set, sorted({tuple(sorted(x)) for x in masks})))
    for A,B in combinations(masks,2):
        def cc(mirror):
            C = [{i} for i in range(len(E))]
            for J in [A,B]+([{i,j} for i,j in enumerate(rev)] if mirror else []):
                merged = [x for x in C if x & J]
                C = [x for x in C if not x & J]+[set.union(*merged)]
            return len(C)
        def monochrome(v): return len({v[j] for j in A}) == len({v[j] for j in B}) == 1
        require(sum(map(monochrome,data))==2**cc(False))
        require(sum(map(monochrome,blue))==2**cc(True))
        cases["joint_shortcut_pairs"]+=1
print(json.dumps(cases,sort_keys=True))
