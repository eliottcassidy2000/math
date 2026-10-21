#!/usr/bin/env python3
from collections import Counter
from itertools import product
from pathlib import Path
import json
from check_balanced_group_ports import tournament, path, choice, states, infinity, need, pair, edges, reverse_edges

HERE=Path(__file__).resolve().parent

def run():
    reports=[];total=0
    for k,d in product(range(2,5),range(1,4)):
        m=k+d;D=1<<d;q=infinity(k)+1
        deleted=set(range(k,k+d))|set(range(k+d+1,k+2*d+1))
        shrink=lambda v:v if v<k else v-d if v==k+d else v-2*d
        def delete(E):return frozenset((shrink(a),shrink(b)) for a,b in E if a not in deleted and b not in deleted)
        need(delete(tournament(m))==tournament(k),'actual induced family embedding')
        counts=Counter()
        for s in product(states(m),repeat=2):
            parent=tuple(infinity(k) if x==infinity(m) else x%infinity(k) for x in s)
            p=path(m,s);old=path(k,parent)
            literal=tuple(shrink(v) for v in p if v not in deleted)
            need(literal==old,'all-depth path deletion')
            A=set(range(k));B=set(range(k+1,2*k+1));I,K=(choice(k,x) for x in parent)
            J={2*k-v for v in A-K};Ar=A-I;Br=B-J
            pre=(max(I) if I else None,min(J) if J else k)
            post=(max(Ar) if Ar else k,min(Br) if Br else None)
            E,H=(D-1 if x==infinity(m) else x//infinity(k) for x in s)
            C=set()
            if (E!=0 or H!=D-1) and None not in pre:C.add(pair(*pre))
            if (E!=D-1 or H!=0) and None not in post:C.add(pair(*post))
            lhs=delete(reverse_edges(tournament(m),edges(p)))
            rhs=reverse_edges(tournament(k),edges(old)^C)
            need(lhs==rhs,'all-depth two-seam target identity')
            counts[len(C)]+=1
        expected=Counter({0:D*D+2*(q-2),1:2*(q-2)*D*D+2*D+2*(q-2)**2,
                          2:(q-2)**2*(D*D-2)+2*(q-2)*(D-1)+1})
        need(counts==expected,'all-depth defect census')
        need(sum(counts.values())==(D*(q-1)+1)**2,'all-depth total')
        reports.append(dict(parent_k=k,depth=d,source_order=2*m+1,
                            source_paths=sum(counts.values()),defect_counts=dict(sorted(counts.items()))))
        total+=sum(counts.values())
    result=dict(status='FINITE-EXACT local checker for all-order two-seam theorem',source_path_checks=total,reports=reports)
    (HERE/'all_depth_ports.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))

if __name__=='__main__':run()
