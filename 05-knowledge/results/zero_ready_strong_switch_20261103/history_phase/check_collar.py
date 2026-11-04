"""Independent inner-collar control, on the induced B14 (not full outer T)."""
from pathlib import Path
from itertools import combinations
import json,subprocess,tempfile
from validate import Engine,need,read,induced,strong,digest
HERE=Path(__file__).resolve().parent
PACKET=HERE.parent
NAMES=[1,2,4,6,7,8,14,21,22,31,32,33,40,54]
PORTS=[32,33,22]
def inner(rows,names):return induced(rows,[names.index(a) for a in NAMES])
def collar(T):
    r=NAMES.index(32);p=NAMES.index(40);C=[i for i in range(14) if i not in (r,p)]
    incoming=[NAMES[i] for i in range(14) if T[i]>>r&1]
    return dict(incoming32=incoming,prefix40_dominates_remainder=all(T[p]>>i&1 for i in C),remainder_strong=strong(induced(T,C)),zero_ready_collar=(incoming==[33] and all(T[p]>>i&1 for i in C) and strong(induced(T,C))))
def all_roots(E,T):
    q=E.query(T,list(range(14)));F=[z[3] for z in q['port_data']];h=q['H'];need(all(f>=h for f in F),'all rooted split injections E>=H')
    return dict(H=h,E2_by_name=dict(zip(map(str,NAMES),F)),defect_by_name={str(n):f-h for n,f in zip(NAMES,F)},zero_ready_names=[n for n,f in zip(NAMES,F) if f==h])
def support(edges):
    vertices=set(sum(edges,[]));adj={a:set() for a in vertices}
    for a,b in edges:adj[a].add(b);adj[b].add(a)
    left=set(vertices);components=[]
    while left:
        seen={min(left)}
        while True:
            nxt=seen|set().union(*(adj[a] for a in seen))
            if nxt==seen:break
            seen=nxt
        left-=seen;components.append(sorted(seen))
    degrees=sorted(map(len,adj.values()))
    need(degrees==[1,1,1,1,2] and sorted(map(len,components))==[2,3],'literal suffix support P3 disjoint K2')
    return dict(edges=edges,degrees=degrees,components=sorted(components,key=lambda x:(len(x),x)))
def main():
    f=PACKET/'literal/realized_generation2_seed4_depth3_all_prefixes.json';d=read(f)[1]
    parent=read(PACKET/'literal/realized_depth3.json')[4]
    need(d['source_rows']==parent['target_rows'] and d['source_names']==parent['target_names'],'exact first-generation +4 parent')
    records=[]
    with tempfile.TemporaryDirectory(prefix='inner-collar-audit-') as td:
        exe=Path(td)/'endpoint.exe';subprocess.run(['g++','-O2','-std=c++20',str(HERE/'endpoint_backward.cpp'),'-o',str(exe)],check=True);E=Engine(exe)
        for side in ('source','target'):
            T=inner(d[side+'_rows'],d[side+'_names']);q=E.query(T,[NAMES.index(a) for a in PORTS]);c=collar(T)
            joint=[]
            for i,j in combinations(range(3),2):joint.append(dict(ports=[PORTS[i],PORTS[j]],P2=q['P2'],port_data=[q['port_data'][i],q['port_data'][j]],ordered_W2=[q['ordered_W2'][a][b] for a in (i,j) for b in (i,j)]))
            records.append(dict(stage=side,outer_H=d[side+'_H'] if side=='source' else d['target_H'],carrier='induced outer strong suffix B14',names=NAMES,rows=T,ports=PORTS,actual_inports={str(p):[NAMES[i] for i in range(14) if T[i]>>NAMES.index(p)&1] for p in PORTS},H=q['H'],E32=q['port_data'][0][3],root32_defect=q['port_data'][0][3]-q['H'],all_ports=q,joint_pairs=joint,collar=c,all_roots=all_roots(E,T)))
    need(records[0]['H']==46073 and records[1]['H']==80319,'correct B14 Hamilton counts')
    need(records[0]['root32_defect']==0 and records[1]['root32_defect']>0,'zero-ready inner collar genuinely lost')
    need(records[0]['collar']['incoming32']==[33] and records[1]['collar']['incoming32']==[22,33],'second physical inport at22')
    need(records[0]['all_roots']['zero_ready_names']==[32] and not records[1]['all_roots']['zero_ready_names'],'no alternate zero-ready mark in final B14')
    controls=[];pending=[];paths=[PACKET/'literal'/name for name in ('realized_generation2_seed1_depth3_all_prefixes.json','realized_generation2_seed4_depth3_all_prefixes.json')]
    for path,indices in zip(paths,([3,5],[0,1])):
        data=read(path)
        for index in indices:
            z=data[index];A=inner(z['source_rows'],z['source_names']);B=inner(z['target_rows'],z['target_names']);edges=[e for e in z['edits'] if all(a in NAMES for a in e)]
            actual=[[NAMES[a],NAMES[b]] for a,b in combinations(range(14),2) if ((A[a]^B[a])>>b)&1]
            need({frozenset(e) for e in edges}=={frozenset(e) for e in actual},'literal suffix edits')
            before,after=collar(A),collar(B);expected=z['target_H']!=1537279
            need(before['zero_ready_collar'] and after['zero_ready_collar']==expected,'three preserve, one breaks exact inner collar')
            controls.append(dict(file=path.name,index=index,source_H=z['source_H'],target_H=z['target_H'],support=support(edges),before=before,after=after))
            pending.append(B)
    with tempfile.TemporaryDirectory(prefix='all-inner-roots-audit-') as td:
        exe=Path(td)/'endpoint.exe';subprocess.run(['g++','-O2','-std=c++20',str(HERE/'endpoint_backward.cpp'),'-o',str(exe)],check=True);E=Engine(exe)
        for row,T in zip(controls,pending):
            row['all_roots']=all_roots(E,T);need(row['all_roots']['zero_ready_names']==([] if row['target_H']==1537279 else [32]),'complete four-child zero-ready marking control')
    pins=[f,PACKET/'literal/realized_depth3.json',*paths,PACKET/'DEPTH.md',PACKET/'CLOCKS.md']
    out=dict(status='FINITE_EXACT_INDEPENDENT_BACKWARD_INNER_COLLAR_CONTROL',records=records,all_four_edit_supports=controls,source_pins={str(p.relative_to(PACKET)):digest(p) for p in pins},scope='H and E32 here belong to B14, not the outer512425/1537279 tournaments. All four suffix edit supports are P3 disjoint K2, but the support isomorphism type does not preserve inner collar incidence. No original field-gauge receipt is reinterpreted.')
    (HERE/'collar.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8',newline='\n')
    print('PASS: B14 collar32/in33/prefix40; exact defect '+str(records[0]['root32_defect'])+' -> '+str(records[1]['root32_defect'])+'; joint ports32,33,22; all four P3+K2 supports; all14 roots checked in each child, only1537279 has none zero-ready')
if __name__=='__main__':main()
