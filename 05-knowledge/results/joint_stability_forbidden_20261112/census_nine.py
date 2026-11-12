"""Complete marked unit/stable-port census on the inherited order9 universe.

Phase1 reconstructs every endpoint table and actual selection orbits.
Phase2 uses the proved selected-endpoint degree bound to certify whole axes.
All arithmetic is integral; largest literal order is13 in phase2.
"""
from pathlib import Path
from collections import Counter
from itertools import combinations
import csv,gzip,json,hashlib,sys,time
import numpy as np
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'five_stress_connection_atlas_20261112'))
from extend_round import tables
from tournament_tools import strong
from cyclic_port_probe import ib,triangles
def need(x,msg):
    if not x:raise RuntimeError(msg)
def save(name,data):
    if not __debug__:name=name.replace('.json.gz','.optimized.json.gz')
    (HERE/name).write_bytes(gzip.compress((json.dumps(data,separators=(',',':'))+'\n').encode(),mtime=0))
def read(name):return json.loads(gzip.decompress((HERE/name).read_bytes()))
def automorphisms(rows):
    n=len(rows);degrees=[r.bit_count() for r in rows]
    cand=[[v for v in range(n) if degrees[v]==degrees[u]] for u in range(n)]
    order=sorted(range(n),key=lambda u:(len(cand[u]),u));mapping=[-1]*n;out=[]
    def visit(depth,used):
        if depth==n:out.append(tuple(mapping));return
        u=order[depth]
        for v in cand[u]:
            if used>>v&1:continue
            if any(mapping[a]>=0 and bool(rows[u]>>a&1)!=bool(rows[v]>>mapping[a]&1) for a in range(n)):continue
            mapping[u]=v;visit(depth+1,used|1<<v);mapping[u]=-1
    visit(0,0);return out
def move(mask,p):return sum(1<<p[i] for i in range(len(p)) if mask>>i&1)
def phase1():
    source=HERE.parent/'blue_quotient_20261020b/n9_classes.tsv.gz';start=time.time();records=[];mass=0;unit_sources=unit_masks=unit_orbits=0;hvalues=set();aut_checks=0;strong_count=0
    with gzip.open(source,'rt') as f:
        for idx,item in enumerate(csv.DictReader(f,delimiter='\t'),1):
            code=int(item['code']);a=int(item['Aut']);rows=[0]*9;b=35
            for i in range(9):
                for j in range(i+1,9):
                    if code>>b&1:rows[i]|=1<<j
                    else:rows[j]|=1<<i
                    b-=1
            H,U,W,masks=tables(np.array(rows,dtype=np.int64));H=int(H);hvalues.add(H)
            need(H%a==0 and H%2==1,'parity/automorphism control');mass+=H//a
            if int(item['code'])==int(item['converse_code']):need(H==int(item['SC_H']),'inherited self-converse count')
            unit_sources+=bool(masks);unit_masks+=len(masks);need(len(masks)%a==0,'free marking action');unit_orbits+=len(masks)//a
            isstrong=strong(rows);strong_count+=isstrong
            if masks and isstrong:
                aut=automorphisms(rows) if a>1 else [tuple(range(9))];need(len(aut)==a,'actual automorphism group');aut_checks+=a>1
                reps=sorted({min(move(int(m),p) for p in aut) for m in masks});need(len(reps)==len(masks)//a,'marked quotient')
                records.append({'code':code,'rows':rows,'H':H,'Aut':a,'marks':[{'selected':m,'U':sum(int(U[v]) for v in range(9) if m>>v&1)} for m in reps]})
            if idx%10000==0:print(json.dumps({'phase':1,'classes':idx,'strong_unit_classes':len(records),'seconds':round(time.time()-start,2)}),flush=True)
    need(idx==191536 and mass==2**28,'complete ordinary universe/mass')
    need((unit_sources,unit_masks,unit_orbits,len(hvalues))==(17423,19011,18849,1520),'independent inherited C++ aggregates')
    out={'scope':'Every inherited ordinary order9 tournament class, all nonempty selected sets, actual marked isomorphism quotient; output retains all strong unit seeds. No converse merging.','source':str(source.relative_to(HERE.parents[2])).replace('\\','/'),'source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'ordinary_classes':idx,'strong_classes':strong_count,'all_unit_sources':unit_sources,'all_unit_masks':unit_masks,'all_unit_orbits':unit_orbits,'strong_unit_sources':len(records),'strong_unit_orbits':sum(len(r['marks']) for r in records),'labelled_path_mass':mass,'actual_nontrivial_aut_checks':aut_checks,'records':records}
    save('unit_nine.json.gz',out);print(json.dumps({k:v for k,v in out.items() if k!='records'}),flush=True)
def defect(h,w,selected):return int(h)-sum(int(w[a,b]) for a,b in combinations(selected,2))
def axis_values(rows,mask,v,max_m):
    vals=[]
    for m in range(1,max_m+1):
        labels,rr=ib.inflate(rows,v,m);sel=[i for i,q in enumerate(labels) if mask>>q&1];h,u,w,_=tables(np.array(rr,dtype=np.int64));vals.append(defect(h,w,sel))
    return vals
def phase2():
    source=read('unit_nine.json.gz');start=time.time();records=[];hist=Counter();false_screens=[];literal_calls=0;births=0
    for idx,item in enumerate(source['records']):
        rows=item['rows'];n=9;marks=item['marks'];stable={r['selected']:[] for r in marks};cert={r['selected']:{} for r in marks};screen={r['selected']:[] for r in marks}
        for v in range(n):
            candidates=[r['selected'] for r in marks if not r['selected']>>v&1]
            if not candidates:continue
            # The endpoint capacity bound is independent of the mark set.
            R=1+min(rows[v].bit_count(),n-1-rows[v].bit_count());values={m:[1] for m in candidates}
            for size in range(2,R+1):
                labels,rr=ib.inflate(rows,v,size);h,u,w,_=tables(np.array(rr,dtype=np.int64));literal_calls+=1;remaining=[]
                for mask in candidates:
                    sel=[i for i,q in enumerate(labels) if mask>>q&1];d=defect(h,w,sel);values[mask].append(d)
                    if size==2 and d==1:screen[mask].append(v)
                    if d==1:remaining.append(mask)
                    elif size>2:false_screens.append({'code':item['code'],'rows':rows,'selected':mask,'port':v,'R':R,'first_failed_length':size,'values':values[mask]})
                candidates=remaining
                if not candidates:break
            for mask in candidates:stable[mask].append(v);cert[mask][v]={'R':R,'values':values[mask]}
        for itemmark in marks:
            mask=itemmark['selected'];P=stable[mask];hist[len(P)]+=1
            if not P:continue
            twins=[list(ab) for ab in combinations([v for v in range(n) if not mask>>v&1],2) if any(v in P for v in ab) and all(bool(rows[ab[0]]>>w&1)==bool(rows[ab[1]]>>w&1) for w in range(n) if w not in ab)]
            births+=not twins
            records.append({'code':item['code'],'rows':rows,'selected':mask,'H':item['H'],'U':itemmark['U'],'stable_ports':P,'cyclic_triples':triangles(rows,P),'eligible_twins':twins,'axis_certificates':cert[mask]})
        if idx%1000==0:print(json.dumps({'phase':2,'at':idx,'total':len(source['records']),'stable_states':len(records),'births':births,'literal_calls':literal_calls,'seconds':round(time.time()-start,2)}),flush=True)
    # Independent full skeleton controls, including the inherited delayed failure.
    controls=[]
    for r in records[:12]+sorted(records,key=lambda r:-len(r['stable_ports']))[:5]:
        O=[i for i in range(9) if r['selected']>>i&1]
        for v in r['stable_ports']:
            h,u,w=ib.coefficients(r['rows'],v);co=[h[j]-sum(w[j][a][b] for a,b in combinations(O,2)) for j in range(1,10)]
            need(co==[1]+[0]*8,'independent full one-port skeleton');controls.append([r['code'],r['selected'],v])
    hostile=axis_values([90,28,9,16,36,79,30],72,2,3);need(hostile==[1,1,-23],'delayed-failure hostile')
    out={'scope':'Complete strong marked unit seeds of order9; each unselected port tested through the proved selected-endpoint run bound R=1+min(d+,d-). R successful positive values imply all-length stability. All marked types counted separately from their converse.','unit_input_sha256':hashlib.sha256((HERE/'unit_nine.json.gz').read_bytes()).hexdigest(),'strong_unit_orbits':source['strong_unit_orbits'],'stable_port_histogram':dict(hist),'stable_states':len(records),'stable_clone_births':births,'stable_clone_descendants':len(records)-births,'cyclic_stable_sets':sum(bool(r['cyclic_triples']) for r in records),'false_two_clone_screens':false_screens,'literal_calls':literal_calls,'independent_full_skeleton_controls':controls,'delayed_failure_control':hostile,'records':records}
    save('stable_nine.json.gz',out);print(json.dumps({k:v for k,v in out.items() if k not in ['records','false_two_clone_screens','independent_full_skeleton_controls']}),flush=True)
if __name__=='__main__':
    {'unit':phase1,'stable':phase2}[sys.argv[1]]()
