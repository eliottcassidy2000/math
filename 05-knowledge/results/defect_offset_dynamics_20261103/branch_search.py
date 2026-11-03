"""Source-preserving search, deduplicated only for the next enumeration."""
from ear_cuts import *
from itertools import permutations


def canon_order(T):
    groups={}
    for i,row in enumerate(T):groups.setdefault(row.bit_count(),[]).append(i)
    best=None;best_order=None
    for choices in product(*(permutations(groups[s]) for s in sorted(groups))):
        order=sum((tuple(c) for c in choices),())
        code=induced(T,order)
        if best is None or code<best:best=code;best_order=order
    return best,best_order


def canon(T):return canon_order(T)[0]


def extend(bank,targets):
    result={};attempts=0
    for T,history in bank.items():
        for m,q,h in all_two(T,targets):
            U=birth(birth(T,m),q)
            if not strong(U):continue
            C,order=canon_order(U);attempts+=1
            if C not in result:result[C]=history+[dict(source=T,masks=(m,q),h=h,target=C,canonical_order=order)]
    print('to',targets,'literal',attempts,'classes',len(result),flush=True)
    return result


def main():
    T,_,_=source(((1,1),));bank={canon(T):[]}
    output=[]
    for targets in ({17},{53},{161,163}):
        bank=extend(bank,targets)
        output.append([dict(T=t,history=path) for t,path in bank.items()])
        (HERE/'branches.json').write_text(json.dumps(output,indent=2)+'\n')
        if not bank:break
    for T,path in list(bank.items())[:5]:print('example',path,flush=True)


if __name__=='__main__':main()
