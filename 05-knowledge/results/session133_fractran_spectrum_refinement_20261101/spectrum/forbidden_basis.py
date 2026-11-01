#!/usr/bin/env python3
"""Complete minimal induced obstructions to H<=B for B<=21.
Enumeration through order seven, plus exhaustively bounded SCC words.
Normal forms sort degree cells, then minimize over all within-cell permutations.
"""
import argparse,json,math
from collections import Counter,defaultdict
from itertools import permutations,product
from build_atlas import rows_of,perm_code,path_dp,components,need

def degree_perms(a):
    cells=defaultdict(list)
    for v,row in enumerate(a):cells[row.bit_count()].append(v)
    for pp in product(*(permutations(cells[d]) for d in sorted(cells))):
        yield tuple(v for p in pp for v in p)

def canonical(a):return min(perm_code(a,p) for p in degree_perms(a))
def delete(a,v):
    old=[u for u in range(len(a)) if u!=v]
    return [sum(((a[u]>>w)&1)<<j for j,w in enumerate(old)) for u in old]
def join(blocks):
    a=[];offset=0;n=sum(map(len,blocks))
    for b in blocks:
        following=((1<<n)-1)^((1<<(offset+len(b)))-1)
        a.extend((row<<offset)|following for row in b)
        offset+=len(b)
    return a

def atlas7():
    old={0};classes={}; levels=[]
    for n in range(1,8):
        current=set()
        for code in old:
            base=rows_of(n-1,code)
            for mask in range(1<<(n-1)):
                a=base+[0]
                for u in range(n-1):
                    if mask>>u&1:a[u]|=1<<(n-1)
                    else:a[-1]|=1<<u
                current.add(canonical(a))
        records=[];mass=0
        for code in sorted(current):
            a=rows_of(n,code); h=path_dp(a)
            # Isomorphisms preserve labelled degrees; enumerate all such maps.
            cells=defaultdict(list)
            for v,row in enumerate(a):cells[row.bit_count()].append(v)
            aut=0
            for pp in product(*(permutations(c) for c in cells.values())):
                p=list(range(n))
                for source,target in zip(cells.values(),pp):
                    for u,v in zip(source,target):p[u]=v
                aut+=perm_code(a,p)==code
            need(aut>0 and math.factorial(n)%aut==0,'bad automorphism count')
            mass+=math.factorial(n)//aut
            cards=[path_dp(delete(a,v)) for v in range(n)]
            need(all(x<=h for x in cards),'deletion monotonicity failure')
            records.append({'n':n,'degree_code':code,'rows':a,'H':h,
                            'strong':len(components(a))==1 and n>=3,
                            'automorphisms':aut,'deletion_H':cards})
        need(mass==1<<(n*(n-1)//2),'orbit mass fails')
        classes[n]=records;levels.append({'n':n,'classes':len(current),'labelled_mass':mass})
        old=current
    need([v['classes'] for v in levels]==[1,1,2,4,12,56,456],'small class counts fail')
    for n,want_h,want_classes,want_mass in ((6,9,1,80),(7,15,2,3360)):
        reduced=[row for row in classes[n] if all(len(c)>1 for c in components(row['rows']))]
        floor=min(row['H'] for row in reduced)
        extremals=[row for row in reduced if row['H']==floor]
        mass=sum(math.factorial(n)//row['automorphisms'] for row in extremals)
        need((floor,len(extremals),mass)==(want_h,want_classes,want_mass),'exact core-floor control fails')
    return classes,levels

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--threshold',type=int,default=21)
    args=parser.parse_args();B=args.threshold
    need(1<=B<=21,'supported proven bounded range is 1..21')
    classes,levels=atlas7()
    strong=[]
    for records in classes.values():
        for row in records:
            if row['strong'] and row['H']>B and max(row['deletion_H'])<=B:strong.append(row)
    # B<=21 implies every good strong component has order <=6 (Busch).
    atoms=[row for records in classes.values() for row in records if row['strong'] and row['H']<=B]
    r=0;power=1
    while power<=B:r+=1;power*=3
    # A minimal bad word cannot contain a bad atom, and deleting a complete
    # atom yields a proper induced tournament. Hence at most r atoms suffice.
    words=[]
    def walk(seq,prod_h):
        if len(seq)>=2 and prod_h>B:
            hmax=max(math.prod(x['H'] if j!=i else max(x['deletion_H']) for j,x in enumerate(seq)) for i in range(len(seq)))
            if hmax<=B:
                a=join([row['rows'] for row in seq]); cards=[path_dp(delete(a,v)) for v in range(len(a))]
                need(path_dp(a)==prod_h and max(cards)==hmax,'SCC product check fails')
                words.append({'n':len(a),'H':prod_h,'rows':a,'deletion_H':cards,
                    'atoms':[{'n':x['n'],'H':x['H'],'degree_code':x['degree_code']} for x in seq]})
            return
        if len(seq)==r:return
        for row in atoms:walk(seq+[row],prod_h*row['H'])
    walk([],1)
    total=strong+words
    candidates=[]
    for b in range(3):
        if 5**b>B:continue
        a=0
        while 3**(a+1)*5**b<=B:a+=1
        candidates.append(3*a+4*b)
    K=max(candidates)
    for row in total:
        cs=components(row['rows'])
        if all(len(c)==3 for c in cs):
            need(len(cs)==r,'pure triangle obstruction is not the sharp chain')
        else:need(row['n']<=K+1,'exact-core sharpening fails')
    need(all(x['n']<=3*r for x in total),'sharp general obstruction cutoff fails')
    triangle=rows_of(3,2)
    sharp=join([triangle]*r)
    need(path_dp(sharp)==3**r and max(path_dp(delete(sharp,v)) for v in range(3*r))==3**(r-1),'sharpness hostile fails')
    if B==21:
        need(Counter(tuple(x['n'] for x in w['atoms']) for w in words)=={(3,5):6,(5,3):6,(4,4):1,(3,3,3):1},'disconnected classification fails')
    payload={'schema':'session133.minimal-induced-H-obstructions.v1','threshold':B,
             'scope':'All orders, conditional only on cited OCF, Moon strong-deletion theorem, and Busch bound; finite classification exact.',
             'canonical_code':'Minimum upper-triangle bitcode subject to nondecreasing outdegree order.',
             'universe':levels,'strong_obstructions':strong,'nonstrong_obstructions':words,
             'counts':{'strong':len(strong),'nonstrong':len(words),'total':len(total),
               'by_order':dict(sorted(Counter(x['n'] for x in total).items())),
               'strong_by_order':dict(sorted(Counter(x['n'] for x in strong).items())),
               'strong_by_H':dict(sorted(Counter(x['H'] for x in strong).items()))},
             'sharp_general_maximum_order':3*r,'nontriangle_chain_order_bound':K+1,'checks':{'labelled_mass':'all orders 1..7','maximal_order_control':'ordered r triangle chain','all_vertex_deletions':'every listed obstruction'}}
    print(json.dumps(payload,indent=2,sort_keys=True))
if __name__=='__main__':main()
