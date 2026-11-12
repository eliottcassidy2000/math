"""Full simultaneous stable-port tensors for the complete order8 strong bank."""
from pathlib import Path
from itertools import combinations,product
import json,sys,time,gzip
import numpy as np
import sympy as s
from extend_round import tables
from c3_pair_family import exponentials,evaluate
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'seed_recursion_births_20261111'))
import cones
def need(x,msg):
    if not x:raise RuntimeError(msg)
def run():
    data=json.loads((HERE/'cyclic_port_receipt.json').read_text());records=[];parts={};t=time.time()
    source=[r for r in data['positive_screen_certificates'] if r['n']==8 and r['individually_stable']]
    for count,r in enumerate(source,1):
        rows=r['rows'];ports=r['individually_stable'];O=[v for v in range(8) if r['selected']>>v&1]
        need(not r['cyclic_triples'],'acyclic ports required for finite kernel')
        d,u=cones.kernel(rows,ports,O);de=exponentials(d);ue=exponentials(u)
        stable=d=={(1,)*len(ports):1}
        child,sel,_=cones.expand(rows,O,ports,[2]*len(ports));h,cu,w,_=tables(np.array(child,dtype=np.int64))
        need(cones.eval_kernel(d,[2]*len(ports))==int(h)-sum(int(w[a,b]) for a,b in combinations(sel,2)),'joint literal defect')
        need(cones.eval_kernel(u,[2]*len(ports))==sum(int(cu[v]) for v in sel),'joint literal slope')
        item={**r,'ports':ports,'jointly_stable':stable,'defect_terms':[[*k,v] for k,v in sorted(d.items())],'U_terms':[[*k,v] for k,v in sorted(u.items())],'U_exponential_terms':[[*k,v] for k,v in sorted(ue.items())]}
        records.append(item)
        if stable:
            grid={(1,)*len(ports)}
            for i in range(len(ports)):
                for a in range(1,25):
                    point=[1]*len(ports);point[i]=a;grid.add(tuple(point))
            for i,j in combinations(range(len(ports)),2):
                for a,b in product(range(1,13),repeat=2):
                    point=[1]*len(ports);point[i]=a;point[j]=b;grid.add(tuple(point))
            for a in range(1,13):grid.add((a,)*len(ports))
            for point in grid:
                U=evaluate(ue,point);odd=U;j=0
                while odd%2==0:odd//=2;j+=1
                if odd<=4097:
                    record={'core_record':len(records)-1,'sizes':point,'U':U,'v2_U':j,'order':8+sum(a-1 for a in point)}
                    if odd not in parts or record['order']<parts[odd]['order']:parts[odd]=record
        if count%20==0:print(json.dumps({'processed':count,'total':len(source),'seconds':round(time.time()-t,2)}),flush=True)
    out={'scope':'Every order8 strong marked unit seed having at least one individually stable unselected port, modulo selected-set automorphisms. Full simultaneous finite tensors; every complete port set is acyclic. Literal all-ports-doubled controls. Bounded slope sample: each single length1..24; each pair length1..12 with other lengths1; uniform lengths1..12.','core_states':len(records),'jointly_stable':sum(r['jointly_stable'] for r in records),'records':records,'sampled_odd_parts_le4097':{str(k):v for k,v in sorted(parts.items())}}
    tag='.optimized' if not __debug__ else ''
    raw=(json.dumps(out,separators=(',',':'))+'\n').encode()
    (HERE/('joint_port_extension'+tag+'.json.gz')).write_bytes(gzip.compress(raw,mtime=0))
    forms=[]
    for i,r in enumerate(records):
        poly=0
        for term in r['U_exponential_terms']:
            value=term[-1]
            for j,base in enumerate(term[:-1]):
                for prime,e in s.factorint(base).items():value*=s.Symbol(f'x{j}_{prime}')**e
            poly+=value
        forms.append({'record':i,'U':str(s.factor(poly))})
    (HERE/('joint_port_forms'+tag+'.json')).write_text(json.dumps(forms,indent=2)+'\n')
    print(json.dumps({'states':len(records),'jointly_stable':out['jointly_stable'],'odd_parts':len(parts),'requested':{c:parts.get(c) for c in [1,3,5,7,11,13,19,27,33,1039]},'failures':[(r['rep'],r['selected'],r['ports']) for r in records if not r['jointly_stable']]}),flush=True)
if __name__=='__main__':run()
