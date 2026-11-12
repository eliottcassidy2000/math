"""All-length stable-face cylinder union for the 300 births through nine."""
from pathlib import Path
from collections import defaultdict,Counter
from itertools import combinations
import importlib.util,sys,json,gzip
import sympy as s
from census_nine import read,save,need
from joint_nine import quotient,key
from birth_cones import canonical_orders
HERE=Path(__file__).resolve().parent
OLD=HERE.parent/'five_stress_connection_atlas_20261112'
spec=importlib.util.spec_from_file_location('s288_cylinders',HERE.parent/'session288_twin_lengths_20261112/stable/cones.py')
cc=importlib.util.module_from_spec(spec);spec.loader.exec_module(cc)

def run():
    old=json.loads((OLD/'cyclic_port_receipt.json').read_text())['positive_screen_certificates']
    dag=json.loads((OLD/'contraction_dag.json').read_text())['records']
    roots={key(r['rows'],r['selected']) for r in dag if not r['children']}
    births=[{'rows':r['rows'],'selected':r['selected'],'ports':r['individually_stable'],'faces':[r['individually_stable']]} for r in old if key(r['rows'],r['selected']) in roots]
    nine=read('stable_nine.json.gz')['records'];joint=read('joint_nine.json.gz')['records']
    for r,j in zip(nine,joint):
        if r['eligible_twins']:continue
        P=r['stable_ports'];faces=[P] if j['joint_stable'] else [[p] for p in P]
        births.append({'rows':r['rows'],'selected':r['selected'],'ports':P,'faces':faces})
    need(len(births)==300,'complete birth count')
    groups=defaultdict(list);records=[];hist=Counter()
    for idx,r in enumerate(births):
        rows=r['rows'];mask=r['selected'];q,o,comps=quotient(rows,mask);canonical,orders=canonical_orders(q,o)
        need(len(orders)==1,'rigid marked quotient for this bank')
        order=orders[0];hist[(len(rows),len(r['ports']),len(r['faces']))]+=1
        charts=[]
        for face in r['faces']:
            free={i for i,c in enumerate(comps) if any(v in face for v in c)}
            need(all(len(comps[i])==1 for i in free),'birth free block singleton')
            chart=tuple(0 if i in free else len(comps[i]) for i in order)
            groups[canonical].append((idx,chart));charts.append(chart)
        records.append({**r,'birth_id':idx,'canonical_quotient':canonical,'charts':charts})
    cross=[];within=[];terms=Counter();z=s.Symbol('z')
    for q,charts in groups.items():
        for (i,a),(j,b) in combinations(charts,2):
            c=cc.meet(a,b)
            if c is not None:(within if i==j else cross).append({'first_birth':i,'second_birth':j,'intersection':c,'minimum_order':sum(c)+c.count(0)})
        for c,m in cc.union_terms({c for _,c in charts}).items():terms[(sum(c)+c.count(0),c.count(0))]+=m
    gf=sum(m*z**b/(1-z)**f for (b,f),m in terms.items());series=s.series(gf,z,0,41).removeO().expand();counts=[int(series.coeff(z,n)) for n in range(1,41)]
    need(counts[:9]==[0,0,1,1,5,13,38,103,324],'complete low-order counts')
    out={'scope':'Exact all-length union of maximal stable-face cylinders of all 300 strong marked stable-clone births through nine. Only the initially certified stable faces are varied; adaptive newly stable ports are not assumed absent.','births':len(births),'maximal_cylinders':sum(len(r['faces']) for r in births),'quotient_types':len(groups),'birth_order_rank_faces_histogram':[[*k,v] for k,v in sorted(hist.items())],'cross_birth_intersections':cross,'within_birth_intersections':within,'generating_terms':[[b,f,m] for (b,f),m in sorted(terms.items()) if m],'GF':str(s.factor(gf)),'counts_n1_to40':counts,'records':records}
    if not cross:
        need(all(counts[n-1]==4*n*n+n-7 for n in range(10,41)),'quadratic family tail')
        tail=sum(counts[n-1]*z**n for n in range(1,10))+4*(z*s.diff(z*s.diff(z**10/(1-z),z),z))+z*s.diff(z**10/(1-z),z)-7*z**10/(1-z)
        need(s.cancel(gf-tail)==0,'exact quadratic tail identity')
        out['for_every_n_ge10']='4*n*n+n-7'
    save('birth_union.json.gz',out);print(json.dumps({k:v for k,v in out.items() if k not in ['records','counts_n1_to40']}),flush=True)
if __name__=='__main__':run()
