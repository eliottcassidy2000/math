"""Full simultaneous defect and U kernels on the inherited stable core ports."""
from pathlib import Path
from collections import defaultdict
from itertools import product,combinations
import json
import port_transport as pt

HERE=Path(__file__).resolve().parent


def ends(rows,marked):
    n=len(rows);old=[v for v in range(n) if v not in marked]
    order=old+list(marked);m=len(marked);r=len(old);zero=(0,)*m
    # Caller checks the marked subtournament is acyclic. One visit per
    # marked symbol between successive old vertices; two covers need r+2.
    bound=r+2
    dp=[{} for _ in range(1<<r)]
    blank=[0]*n
    for mask in range(1<<r):
        for counts in product(range(bound+1),repeat=m):
            value=[0]*n
            for v in range(r):
                if not mask>>v&1:continue
                rest=mask^(1<<v)
                if rest==0 and counts==zero:value[v]=1
                else:
                    prev=dp[rest].get(counts,blank)
                    value[v]=sum(prev[u] for u in range(n) if rows[order[u]]>>order[v]&1)
            for b in range(m):
                if not counts[b]:continue
                previous=tuple(counts[i]-(i==b) for i in range(m))
                v=r+b
                prev=dp[mask].get(previous,blank)
                value[v]=int(mask==0 and previous==zero)
                value[v]+=sum(prev[u] for u in range(n) if rows[order[u]]>>order[v]&1)
            if any(value):dp[mask][counts]=value
    return order,dp


def kernel(rows,marked,selected):
    order,end=ends(rows,marked)
    n=len(rows);reverse=[sum(1<<b for b in range(n) if rows[b]>>a&1) for a in range(n)]
    _,start=ends(reverse,marked)
    full=len(end)-1
    sel=[i for i,v in enumerate(order) if v in selected]
    d=defaultdict(int);u=defaultdict(int)
    for counts,values in end[full].items():d[counts]+=sum(values)
    for counts,values in start[full].items():u[counts]+=sum(values[i] for i in sel)
    for mask in range(full+1):
        lefts=[(a,sum(vals),sum(vals[i] for i in sel)) for a,vals in end[mask].items()]
        rights=[(b,sum(vals[i] for i in sel)) for b,vals in start[full^mask].items()]
        for a,total,cut in lefts:
            if not total:continue
            for b,s in rights:
                if not s:continue
                counts=tuple(x+y for x,y in zip(a,b))
                d[counts]-=cut*s
                u[counts]+=total*s
    return ({k:v for k,v in d.items() if all(k) and v},
            {k:v for k,v in u.items() if all(k) and v})


def eval_kernel(poly,sizes):
    total=0
    for counts,value in poly.items():
        term=value
        for a,j in zip(sizes,counts):term*=pt.ib.onto(a,j)
        total+=term
    return total


def expand(rows,selected,ports,sizes):
    lengths=dict(zip(ports,sizes))
    labels=[v for v in range(len(rows)) for _ in range(lengths.get(v,1))]
    out=[0]*len(labels)
    for i,j in combinations(range(len(labels)),2):
        forward=labels[i]==labels[j] or bool(rows[labels[i]]>>labels[j]&1)
        out[i if forward else j] |= 1<<(j if forward else i)
    return out,[i for i,v in enumerate(labels) if v in selected],labels


def main():
    source=json.loads((HERE.parent/'endpoint_coverage_compiler_20261111/inflation_bank_7.json').read_text())['records']
    groups=defaultdict(list)
    for item in source:
        if len(pt.cp.components(item['rows']))==1:groups[item['n'],item['code'],item['selected']].append(item)
    records,controls=[],0
    for (n,code,mask),items in groups.items():
        rows=items[0]['rows'];ports=sorted(x['inflated_vertex'] for x in items)
        selected=[v for v in range(n) if mask>>v&1]
        for tri in combinations(ports,3):
            pt.ib.tc.require(not all(sum(bool(rows[v]>>w&1) for w in tri if w!=v)==1 for v in tri),'marked set acyclic')
        d,u=kernel(rows,ports,selected)
        stable=d=={(1,)*len(ports):1}
        sizes=[2]*len(ports)
        expanded,es,_=expand(rows,selected,ports,sizes)
        h,eu,ew,_,_,_=pt.ib.tc.path_tables(expanded)
        actual_d=h-sum(ew[a][b] for a,b in combinations(es,2))
        pt.ib.tc.require(eval_kernel(d,sizes)==actual_d,'full simultaneous defect')
        pt.ib.tc.require(eval_kernel(u,sizes)==sum(eu[v] for v in es),'full simultaneous slope')
        controls+=1
        records.append(dict(n=n,code=code,rows=rows,selected=mask,ports=ports,jointly_stable=stable,
                            defect_terms=[[*k,v] for k,v in sorted(d.items())],
                            slope_terms=[[*k,v] for k,v in sorted(u.items())]))
    report=dict(core_states=len(records),jointly_stable=sum(r['jointly_stable'] for r in records),
                literal_controls=controls,records=records)
    (HERE/'cones.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf8')
    print(json.dumps({k:v for k,v in report.items() if k!='records'},indent=2))
    print(json.dumps({'failures':[r for r in records if not r['jointly_stable']]},indent=2))


if __name__=='__main__':main()
