#!/usr/bin/env python3
from itertools import combinations,permutations
from functools import lru_cache
import json
from independent import rows,modules,flip,polynomial_bank,need


def main():
    overlaps=uniform_nested=0;max_h_gap=0
    for n in range(1,6):
        B,bank=polynomial_bank(n);full=(1<<n)-1
        def value(code,a,c):
            mask=a|(c<<n)
            return sum(v for s,v in zip(B,bank[code]) if s&mask==s)
        for code in range(len(bank)):
            adj=rows(n,code);mods=modules(adj)
            @lru_cache(None)
            def h(S):
                names=[u for u in range(n) if S>>u&1]
                return sum(all(adj[u]>>v&1 for u,v in zip(w,w[1:])) for w in permutations(names))
            for A,C in combinations(mods,2):
                I=A&C
                if I and I not in (A,C):
                    U=A|C;need(U in mods,'union module')
                    target=code^flip(n,A)^flip(n,C)
                    quotient=full^U|(U&-U)
                    lower=2*h(quotient)*h(U)
                    gap=bank[target][0]-bank[code][0]
                    need(gap>=lower>0,'strict crossing profile consequence')
                    overlaps+=1;max_h_gap=max(max_h_gap,gap)
                elif I in (A,C):
                    inner,outer=(A,C) if I==A else (C,A)
                    rest=outer^inner
                    orientations={bool(adj[u]>>v&1) for u in range(n) if inner>>u&1 for v in range(n) if rest>>v&1}
                    if len(orientations)!=1:continue
                    forward=orientations.pop()
                    representative=(outer&-outer).bit_length()-1
                    ext=sum(1<<u for u in range(n) if not outer>>u&1 and adj[u]>>representative&1)
                    a=ext|(inner if forward else rest)
                    target=code^flip(n,inner)^flip(n,outer)
                    gap=value(target,a,full)-value(code,a,full)
                    quotient=full^outer|(outer&-outer)
                    lower=2*h(quotient)*h(inner)*h(rest)
                    need(gap>=lower>0,'uniform nested one-port separator')
                    uniform_nested+=1
    need((overlaps,uniform_nested)==(2166,3580),'pair universe')
    return dict(status='PASS',proper_overlapping_original_module_pairs=overlaps,
                uniformly_split_nested_pairs=uniform_nested,
                largest_scalar_overlap_gap=max_h_gap,
                method='Scalar and two-port values from independent literal-word banks; quotient and induced-block Hamilton counts from literal permutations.')


if __name__=='__main__':print(json.dumps(main(),indent=2,sort_keys=True))
