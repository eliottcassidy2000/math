#!/usr/bin/env python3
from pathlib import Path
import json,hashlib,subprocess,sys
import bank_audit as B
import kernel_audit as K
P=Path(__file__).resolve().parent

def main():
    odd=Path(sys.argv[1]).resolve() if len(sys.argv)>1 else P.parent/'session189-odd-transitive';anchor=Path(sys.argv[2]).resolve() if len(sys.argv)>2 else P.parent/'session189-collision-next'
    records=json.loads((odd/'physical_queries.json').read_text());controls=json.loads((odd/'controls.json').read_text());graphs=controls['graphs'];raw=(odd/'query_input.txt').read_text().splitlines();N=int(raw[0]);B.need(len(raw)==N+1,'query input length');parsed=[];requests=[]
    for line in raw[1:]:
        v=list(map(int,line.split()));n=v[0];B.need(len(v)==n+3,'query input record');rows=v[1:1+n];a,c=v[-2:];parsed.append((rows,a,c));q=B.ports(rows,a,c);requests.append(' '.join(map(str,[len(q),*q])))
    # Verify every named graph really has the declared module and local reversal.
    for case in controls['cases']:
        m=case['m'];source=graphs[case['graph_ids'][0]]['rows'];target=graphs[case['graph_ids'][1]]['rows'];coremask=(1<<m)-1;n=len(source)
        B.need(tuple(row&coremask for row in source[:m])==K.tt(m),'physical transitive source order')
        B.need(target==[row^(coremask^(1<<i))if i<m else row for i,row in enumerate(source)],'literal local module converse')
        B.need(case['exterior_rows']==[row>>m for row in source[m:]],'literal exterior rows')
        B.need(all((source[m+j]&coremask)==(coremask if case['into_mask']>>j&1 else 0)for j in range(n-m)),'original module cross incidences')
        if 'witness_groups'in case:
            X,Y=case['witness_groups'];core=list(range(m))+X+Y;rootmask=sum(1<<v for v in core);targetq=graphs[case['graph_ids'][2]]['rows'];want=list(source)
            B.need(all(source[a]>>x&1 for a in range(m)for x in X)and all(source[y]>>a&1 for y in Y for a in range(m)),'source cyclic witness incidences')
            for C in [core,list(range(m)),X,Y]:
                B.need(all(len({(source[v]>>i)&1 for i in C})==1 for v in range(n)if v not in C),'original selected source module')
                mask=sum(1<<v for v in C)
                for i in C:want[i]^=mask^(1<<i)
            B.need(want==targetq,'actual fully covered root-family quotient target')
            union=X+Y
            B.need(all(len({(source[v]>>i)&1 for i in union})==1 for v in range(n)if v not in union and v>=m),'ordered witness module after deleting A')
            B.need(all(source[x]>>y&1 for x in X for y in Y),'ordered witness X to Y')
            child_rows=[]
            for C in [X,Y]:child_rows.append(tuple(sum(1<<j for j,v in enumerate(C)if source[u]>>v&1)for u in C))
            profiles=[K.paths(r)[3].get((1<<len(r))-1,[])for r in child_rows];B.need(profiles[0]==profiles[1]and all(profiles[0][j]**2>=profiles[0][j-1]*profiles[0][j+1]for j in range(1,len(profiles[0])-1)),'complete equal log-concave witness profiles')
    seen=set();combinations={}
    for rec in records:
        case=controls['cases'][rec['case']];src=graphs[case['graph_ids'][0]]['rows'];m=case['m'];outside=case.get('ambient_names',[]);fixed=sum(1<<v for v in outside if src[v]&1);B.need(fixed==rec['outside_fixed_one_mask'],'physical fixed ambient mask');B.need(rec['group'] in case['witness_groups'],'exact physical witness group')
        for v in outside:B.need(all(((src[v]>>i)&1)==(src[v]&1)for i in range(m)),'ambient sees core child uniformly')
        for query in rec['queries']:
            f=query['flags'];a=fixed|(1 if f&1 else 0)|(2 if f&2 else 0);c=fixed|((1<<rec['r'])if f&4 else 0)|(sum(1<<i for i in rec['group'])if f&8 else 0);B.need((a,c)==(query['a_mask'],query['c_mask']),'physical finite-difference labels')
            for g in query['graph_records']:
                idx=g['id'];B.need(idx not in seen,'unique physical record');seen.add(idx);rows,aa,cc=parsed[idx];B.need((aa,cc)==(a,c)and rows==graphs[g['graph_id']]['rows'],'same physical graph and ports');B.need(query['finite_difference_sign']==(-1)**(4-f.bit_count()),'finite difference sign')
    B.need(seen==set(range(N)),'complete grouped query coverage');anchor_data=json.loads((anchor/'controls.json').read_text());expected=list(map(int,(odd/'query_output.txt').read_text().split()));anchor_records=[]
    for case in anchor_data['cases']:
        T=case['core_rows'];U=case['core_target_rows'];C=case['core_converse_rows'];B.need(tuple(C)==K.op(T),'literal whole core converse')
        for index,block in enumerate(case['child_blocks']):
            source_child=tuple(sum(1<<j for j,v in enumerate(block)if T[u]>>v&1)for u in block);target_child=tuple(sum(1<<j for j,v in enumerate(block)if U[u]>>v&1)for u in block)
            if not index:B.need(target_child==K.op(source_child),'literal distinguished child converse')
            sp=K.paths(source_child)[3].get((1<<len(block))-1,[]);tp=K.paths(target_child)[3].get((1<<len(block))-1,[]);B.need(sp==tp,'complete actual matching child profiles')
            for other in case['child_blocks']:
                if other!=block:B.need(all(((T[u]>>v)&1)+((U[u]>>v)&1)==1 for u in block for v in other),'quotient converse on same names')
        for query in case['queries']:
            for gid,kind in zip(case['graph_ids'],['source','target','whole_core_converse']):
                rows=anchor_data['graphs'][gid]['rows'];g=B.ports(rows,query['a_mask'],query['c_mask']);requests.append(' '.join(map(str,[len(g),*g])));expected.append(query[kind]);anchor_records.append({'case':case['name'],'kind':kind,'query':[query['a_mask'],query['c_mask']],'expected':query[kind]})
    subprocess.run(['g++','-O3','-std=c++17',str(P/'count.cpp'),'-o',str(P/'count')],check=True);output=subprocess.check_output([str(P/'count')],input='\n'.join(requests)+'\n',text=True);got=list(map(int,output.split()));B.need(got==expected,'all independent actual grouped/anchor queries')
    for rec in records:
        sums={}
        for query in rec['queries']:
            for g in query['graph_records']:sums[g['kind']]=sums.get(g['kind'],0)+query['finite_difference_sign']*got[g['id']]
        B.need(sums['source']-sums['local_converse']==sums['source']-sums['quotient_only'],'actual saturated quartic restoration');key=(rec['case'],tuple(rec['group']));combinations[key]=combinations.get(key,0)+rec['weight']*(sums['source']-sums['quotient_only'])
    summary=[{'case':case,'group':list(group),'contrast':v}for(case,group),v in combinations.items()];(P/'independent_queries.json').write_text(json.dumps({'odd_actual_counts':got[:N],'grouped_contrasts':summary,'anchor_actual_records':anchor_records},sort_keys=True,indent=2)+'\n')
    receipt={'status':'PASS','independent_actual_grouped_queries':N,'independent_anchor_query_counts':len(got)-N,'total_actual_queries':len(got),'grouped_results':summary,'sources':{f:hashlib.sha256((P/f).read_bytes()).hexdigest()for f in ['queries_audit.py','bank_audit.py','kernel_audit.py','count.cpp']},'artifacts':{'independent_queries.json':hashlib.sha256((P/'independent_queries.json').read_bytes()).hexdigest()},'source_data_pins':{label:{f:hashlib.sha256((src/f).read_bytes()).hexdigest()for f in files}for label,src,files in [('odd',odd,['physical_queries.json','controls.json','query_input.txt','query_output.txt']),('anchor',anchor,['controls.json'])]}}
    (P/'queries_receipt.json').write_text(json.dumps(receipt,sort_keys=True,indent=2)+'\n');print(json.dumps(receipt,sort_keys=True))
if __name__=='__main__':main()
