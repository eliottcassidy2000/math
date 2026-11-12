"""Complete stable-clone descendants at order9, not all order9 unit seeds.

Completeness follows from the stable-twin contraction lemma: any such
descendant contracts to a strong order8 unit seed with a stable port.
"""
from pathlib import Path
from itertools import combinations
from collections import Counter
import json,gzip,sys,time
import numpy as np
from contraction_dag import key
from cyclic_port_probe import ib,triangles
from extend_round import tables
from joint_port_extension import cones
HERE=Path(__file__).resolve().parent
def need(x,msg):
    if not x:raise RuntimeError(msg)
def run():
    t=time.time();bank=json.loads((HERE/'cyclic_port_receipt.json').read_text())['positive_screen_certificates']
    source=[r for r in bank if r['n']==8 and r['individually_stable']]
    dag=json.loads((HERE/'contraction_dag.json').read_text())['records']
    old={key(r['rows'],r['selected']):r for r in dag if r['n']==8}
    generated={};presentations=0
    for r in source:
        rows=r['rows'];O=[v for v in range(8) if r['selected']>>v&1]
        for v in r['individually_stable']:
            labels,child=ib.inflate(rows,v,2);sel=sum(1<<j for j,q in enumerate(labels) if q in O)
            k=key(child,sel);presentations+=1
            generated.setdefault(k,[]).append({'parent_rows':rows,'parent_selected':r['selected'],'inflated_port':v})
    records=[];hist=Counter();birth_orders=Counter();edges=0
    for idx,(k,parents) in enumerate(sorted(generated.items())):
        n,mask,*rows=k;O=[v for v in range(n) if mask>>v&1]
        H,U,W,_=tables(np.array(rows,dtype=np.int64))
        need(int(H)-sum(int(W[a,b]) for a,b in combinations(O,2))==1,'generated unit defect')
        stable=[];false=[];coefficients={}
        for v in range(n):
            if v in O:continue
            labels,child=ib.inflate(rows,v,2);sel=[j for j,q in enumerate(labels) if q in O]
            h,u,w,_=tables(np.array(child,dtype=np.int64))
            if int(h)-sum(int(w[a,b]) for a,b in combinations(sel,2))!=1:continue
            h,u,w=ib.coefficients(rows,v)
            d=[h[j]-sum(w[j][a][b] for a,b in combinations(O,2)) for j in range(1,len(h))]
            coefficients[v]=d
            if d==[1]+[0]*(n-1):stable.append(v)
            else:false.append(v)
        need(stable,'inherited stable clone port missing')
        cyclic=triangles(rows,stable);joint=None;terms=None
        if not cyclic:
            d,u=cones.kernel(rows,stable,O);joint=d=={(1,)*len(stable):1}
            terms=[[*j,c] for j,c in sorted(d.items())]
            child,sel,_=cones.expand(rows,O,stable,[2]*len(stable));h,cu,w,_=tables(np.array(child,dtype=np.int64))
            need(cones.eval_kernel(d,[2]*len(stable))==int(h)-sum(int(w[a,b]) for a,b in combinations(sel,2)),'joint literal defect control')
        children=set();roots=set();new_port_events=[]
        for a,b in combinations([v for v in range(n) if v not in O],2):
            if not(a in stable or b in stable):continue
            if not all(bool(rows[a]>>v&1)==bool(rows[b]>>v&1) for v in range(n) if v not in [a,b]):continue
            kept=[v for v in range(n) if v!=a]
            child=[sum(1<<j for j,v in enumerate(kept) if rows[u]>>v&1) for u in kept]
            sel=sum(1<<j for j,v in enumerate(kept) if v in O);kk=key(child,sel)
            need(kk in old,'contracted strong stable seed absent from order8 bank')
            rr=old[kk];children.add(rr['record']);roots.update(rr['terminal_birth_records'])
            # Check all surviving stable ports on this actual contracted graph.
            for v in stable:
                if v==a:continue
                vv=kept.index(v);h,u,w=ib.coefficients(child,vv)
                selected=[j for j in range(8) if sel>>j&1]
                dd=[h[j]-sum(w[j][i][l] for i,l in combinations(selected,2)) for j in range(1,len(h))]
                if dd!=[1]+[0]*7:new_port_events.append({'contracted_pair':[a,b],'lost_stable_port':v})
        need(children,'descendant without contraction')
        hist[len(stable)]+=1;edges+=len(children)
        for root in roots:birth_orders[dag[root]['n']]+=1
        records.append({'rows':rows,'selected':mask,'H':int(H),'U':sum(int(U[v]) for v in O),'stable_ports':stable,'cyclic_triples':cyclic,'jointly_stable':joint,'joint_defect_terms':terms,'one_port_defects':coefficients,'false_positive_screens':false,'children_order8_records':sorted(children),'terminal_birth_records':sorted(roots),'stable_ports_lost_on_contraction':new_port_events,'generation_presentations':parents})
        if idx%25==0:print(json.dumps({'at':idx,'total':len(generated),'seconds':round(time.time()-t,2)}),flush=True)
    out={'scope':'Every order9 strong marked unit seed admitting an eligible stable-twin contraction, generated from the complete order8 bank. No enumeration or count of order9 stable-clone births. Full stability tests, acyclic joint tensors where applicable, actual marked-isomorphism quotient and contraction graph.','presentations':presentations,'marked_descendant_types':len(records),'stable_port_histogram':dict(hist),'edges_to_order8_types':edges,'terminal_birth_order_histogram':dict(birth_orders),'nonunique_terminal_births':sum(len(r['terminal_birth_records'])>1 for r in records),'cyclic_stable_sets':sum(bool(r['cyclic_triples']) for r in records),'joint_failures':sum(r['jointly_stable'] is False for r in records),'new_stable_port_events':sum(len(r['stable_ports_lost_on_contraction']) for r in records),'records':records}
    tag='.optimized' if not __debug__ else ''
    (HERE/('descendants_nine'+tag+'.json.gz')).write_bytes(gzip.compress((json.dumps(out,separators=(',',':'))+'\n').encode(),mtime=0))
    print(json.dumps({k:v for k,v in out.items() if k!='records'}),flush=True)
if __name__=='__main__':run()
