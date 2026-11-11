"""Exact two-birth search and transitive-module response, standalone via check.py."""
import argparse
import json
import check as a

def all_birth_counts(rows):
    _,end=a.dp(rows)
    n=len(rows)
    full=(1 << n)-1
    reverse=tuple(sum((i != j and not(rows[i] >> j & 1)) << j for j in range(n)) for i in range(n))
    _,start=a.dp(reverse)
    joint=[[0]*n for _ in range(n)]
    for mask in range(1,full):
        for i,e in enumerate(end[mask]):
            if e:
                for j,s in enumerate(start[full^mask]):
                    joint[i][j]+=e*s
    return [sum(start[-1][j] for j in range(n) if out >> j & 1)
            +sum(end[-1][i] for i in range(n) if not(out >> i & 1))
            +sum(joint[i][j] for i in range(n) if not(out >> i & 1)
                 for j in range(n) if out >> j & 1)
            for out in range(full+1)]

def module_runs(rows,module=0):
    """A_r counts paths on other vertices with r separated module placeholders.

    Replacing module by TT_m has H=sum_r A_r*r!*Stirling(m,r), m>=1.
    No adjacent placeholders are allowed, since they are one module run.
    """
    ids=[v for v in range(len(rows)) if v != module]
    rest=a.induced(rows,ids)
    n=len(rest)
    full=(1 << n)-1
    to=[bool(rows[module] >> v & 1) for v in ids]
    back=[bool(rows[v] >> module & 1) for v in ids]
    table={(0,n,1):1}
    for v in range(n):
        table[(1 << v,v,0)]=1
    for mask in range(full+1):
        for r in range(n+1):
            add=sum(table.get((mask,v,r),0) for v in range(n) if back[v])
            if add:
                table[(mask,n,r+1)]=table.get((mask,n,r+1),0)+add
        for v in range(n+1):
            for r in range(n+2):
                count=table.get((mask,v,r),0)
                if not count:
                    continue
                for w in range(n):
                    if not(mask >> w & 1) and (to[w] if v == n else rest[v] >> w & 1):
                        key=(mask | (1 << w),w,r)
                        table[key]=table.get(key,0)+count
    return [sum(table.get((full,v,r),0) for v in range(n+1)) for r in range(n+2)]

def inflate(rows,module,copies):
    labels=[(v,i) for v in range(len(rows)) for i in range(copies if v == module else 1)]
    return tuple(sum((i != j and (u == v and p < q or u != v and bool(rows[u] >> v & 1))) << j
                     for j,(v,q) in enumerate(labels)) for i,(u,p) in enumerate(labels))

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--prefix',type=int,default=1)
    args=parser.parse_args()
    k=args.prefix
    report={'prefix_length':k,'scope':'All ordered two-vertex induced extensions of each specified actual source; no isomorphism filtering.'}
    sources=[]
    for tail in ((2,3,1),(2,2,1)):
        source=a.ready_source(tail,k)
        h=a.dp(source)[0][-1]
        hits={str(3*h+2):[],str(3*h+4):[]}
        for out1 in range(1 << len(source)):
            middle=a.attach(source,out1)
            counts=all_birth_counts(middle)
            for out2,value in enumerate(counts):
                if str(value) not in hits:
                    continue
                target=a.attach(middle,out2)
                a.need(a.dp(target)[0][-1]==value,'quadratic birth enumeration vs independent DP')
                modules=[]
                for m in range(k):
                    runs=module_runs(target,m)
                    a.need(runs[1]==value,'one module copy count')
                    inflated=inflate(target,m,2)
                    a.need(a.dp(inflated)[0][-1]==runs[1]+2*runs[2],'two module copy literal check')
                    good=runs[2]==102*2**(k-1) and not any(runs[3:])
                    modules.append({'prefix_vertex':m,'run_coefficients':runs,'uniform_remaining_prefix_lift':good})
                deletion_counts=[a.dp(a.induced(target,[v for v in range(len(target)) if v != r]))[0][-1]
                                 for r in range(len(target))]
                hits[str(value)].append({'first_outmask':out1,'second_outmask':out2,'rows':target,
                                         'module_responses':modules,'deletion_counts':deletion_counts})
        sources.append({'tail':tail,'rows':source,'h':h,'ordered_two_birth_universe':2**(2*len(source)+1),
                        'hits':hits})
    report['sources']=sources
    print(json.dumps(report,sort_keys=True,indent=2))

if __name__=='__main__':
    main()
