from pathlib import Path
import importlib.util,json
HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('port_boundary_counter',HERE/'check.py');e=importlib.util.module_from_spec(spec);spec.loader.exec_module(e)
contexts=ports=0;ready_by_size={};hostile=[]
for n in range(1,5):
    for B in e.labelled(n):
        source=next((v for v in range(n) if B[v].bit_count()==n-1),None)
        transitive=e.tables(B)[0][-1]==1
        for m in range(2,5):
            for l in range(m):
                T=e.carrier(B,m,l);resource_first=0 if l==0 else 1
                expected=[resource_first] if n==1 else ([m+1+source] if transitive else [])
                actual=[]
                for r in range(len(T)):
                    z=e.state(T,r);indegree=sum(bool(T[v]>>r&1) for v in range(len(T)))
                    if indegree==1 and z[2]-z[1]==1:actual.append(r)
                    ports+=1
                e.need(actual==expected,'complete unit-ready port classification')
                if source is not None and n>1:
                    r=m+1+source;z=e.state(T,r);D=e.tables(B)[0][((1<<n)-1)^(1<<source)]
                    e.need(z[2]-z[1]==D,'global B-source deficit')
                    if not transitive and m==2 and l==0:hostile.append(dict(B=B,root=r,state=z,deletion_source_H=D))
                contexts+=1;e.tables.cache_clear()
result=dict(status='PASS',contexts=contexts,ports=ports,scope='Every labelled B of orders1..4, every resource width2..4 and position; every physical vertex tested for indegree1 and q-s1',nontransitive_B_with_global_source_controls=hostile)
(HERE/('port_readiness.json' if __debug__ else 'port_readiness.optimized.json')).write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8',newline='\n')
print(json.dumps({k:v for k,v in result.items() if k!='nontransitive_B_with_global_source_controls'}))
