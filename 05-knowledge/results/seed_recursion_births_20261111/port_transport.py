"""All positive block sizes: when does a second old port become stable?

Convert each mixed defect column from ordered partitions to integer-base
exponentials. A strict dominant-term bound makes the integer zero search
finite and exact, rather than a guessed cutoff.
"""
from pathlib import Path
from collections import defaultdict,Counter
from itertools import combinations
from math import comb
import json,sys

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'endpoint_coverage_compiler_20261111'))
import joint_ports as jp
import core_compression as cp
import inflation_bank as ib


def exponential(column):
    result=defaultdict(int)
    for i,c in column.items():
        for b in range(1,i+1):result[b]+=c*(-1)**(i-b)*comb(i,b)
    return {b:c for b,c in result.items() if c}


def dominance(poly):
    top=max(poly)
    a=1
    while abs(poly[top])*top**a <= sum(abs(c)*b**a for b,c in poly.items() if b!=top):a+=1
    return a


def evaluate(poly,a):return sum(c*b**a for b,c in poly.items())


def main():
    inherited=json.loads((HERE.parent/'endpoint_coverage_compiler_20261111/inflation_bank_7.json').read_text())['records']
    records=[r for r in inherited if len(cp.components(r['rows']))==1]
    groups=defaultdict(list)
    for r in records:groups[r['n'],r['code'],r['selected']].append(r['inflated_vertex'])
    decisions,activated,cyclic,controls=[],[],[],0
    for key,stable in groups.items():
        n,code,mask=key
        rows=ib.cc.decode(n,code)
        selected=[v for v in range(n) if mask>>v&1]
        for tri in combinations(stable,3):
            if all(sum(bool(rows[v]>>w&1) for w in tri if w!=v)==1 for v in tri):
                cyclic.append(dict(n=n,code=code,selected=mask,stable_ports=stable,triangle=tri))
        for v in stable:
            for w in range(n):
                if w==v or w in selected:continue
                tensor=jp.defect_coefficients(rows,[v,w],selected)
                ib.tc.require({i:d for (i,j),d in tensor.items() if j==1}=={1:1},'first port stable')
                columns=defaultdict(dict)
                for (i,j),d in tensor.items():
                    if j>=2:columns[j][i]=d
                polys={j:exponential(col) for j,col in columns.items()}
                polys={j:p for j,p in polys.items() if p}
                if not polys:
                    sizes='all'
                    bound,certificate=None,None
                else:
                    j=min(polys,key=lambda j:dominance(polys[j]))
                    bound=dominance(polys[j])
                    sizes=[a for a in range(1,bound) if all(evaluate(p,a)==0 for p in polys.values())]
                    certificate=dict(column=j,exponential_coefficients=sorted(polys[j].items()))
                record=dict(n=n,code=code,selected=mask,first_port=v,second_port=w,
                            second_initially_stable=w in stable,stable_sizes=sizes,
                            no_roots_from=bound,certificate=certificate,
                            defect_tensor=[[i,j,d] for (i,j),d in sorted(tensor.items())])
                decisions.append(record)
                if w not in stable and sizes:activated.append(record)
                # Literal child at two to check the actual next port, not just delta.
                labels,child=ib.inflate(rows,v,2)
                old_w=labels.index(w)
                sel=[i for i,label in enumerate(labels) if label in selected]
                h,u,W=ib.coefficients(child,old_w)
                d=[h[j]-sum(W[j][a][b] for a,b in combinations(sel,2)) for j in range(1,len(child)+1)]
                actual=d==[1]+[0]*(len(child)-1)
                ib.tc.require(actual==(sizes=='all' or 2 in sizes),'literal transported port')
                controls+=1
        jp.ends.cache_clear()
    report=dict(strong_stable_port_presentations=len(records),marked_core_states=len(groups),
                stable_port_count_histogram=dict(sorted(Counter(len(s) for s in groups.values()).items())),
                old_port_transports=len(decisions),literal_controls=controls,
                activated=activated,cyclic_stable_sets=cyclic,decisions=decisions)
    (HERE/'port_transport.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf8')
    print(json.dumps({k:v for k,v in report.items() if k not in ('decisions','activated')},indent=2))
    print(json.dumps({'activated':activated[:4]},indent=2))


if __name__=='__main__':main()
