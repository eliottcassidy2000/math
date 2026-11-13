from pathlib import Path
import importlib.util,json
HERE=Path(__file__).resolve().parent
p=HERE.parent.parent/'changing_core_returns_20261113/preparation/check.py'
spec=importlib.util.spec_from_file_location('star_probe_counts',p);e=importlib.util.module_from_spec(spec);spec.loader.exec_module(e)
A=(518,388,456,1011,487,455,387,257,1,502);P=list(A)
for a,b in ((0,2),(0,7),(0,9)):P[a]^=1<<b;P[b]^=1<<a
P=tuple(P);T1691,readyroot=e.birth(P,0,True)
def induced(T,keep):return tuple(sum(1<<j for j,v in enumerate(keep) if T[u]>>v&1) for u in keep)
def probe(T,target):
    n=len(T);hits=[];choices=0;port_data=[]
    for v in range(n):
        old=[u for u in range(n) if u!=v];U=induced(T,old);_,start=e.tables(e.opposite(U));S=start[-1];ones=[old[i] for i,s in enumerate(S) if s==1];zeros=[old[i] for i,s in enumerate(S) if s==0]
        port_data.append(dict(port=v,deletion_H=sum(S),one_start_entries=ones,zero_start_entries=zeros))
        for w in ones:
            for mask in range(1<<len(zeros)):
                O={w}|{u for i,u in enumerate(zeros) if mask>>i&1};V=list(T);V[v]=sum(1<<u for u in O)
                for u in old:
                    if u in O:V[u]&=~(1<<v)
                    else:V[u]|=1<<v
                V=tuple(V);z=e.state(V,v);e.need(z[1]==1,'star classification start budget');choices+=1
                if z[0]+z[3]-1==target:
                    W,r=e.birth(V,v,True);actual=e.state(W,r);e.need(actual[0]==target and actual[2]-actual[1]==1,'literal unit-ready output')
                    hits.append(dict(port=v,outneighbors=sorted(O),star_flips=[u for u in old if ((T[v]^V[v])>>u)&1],parent_rows=V,parent_state=z,parent_strong=e.strong(V),output_rows=W,output_root=r,output_state=actual,output_strong=e.strong(W)))
                e.tables.cache_clear()
    return dict(source_H=e.state(T,0)[0],target=target,source_rows=T,ports=port_data,unique_start_star_choices=choices,hits=hits)
if __name__=='__main__':
    result=dict(status='PASS',cases=[probe(A,1691),probe(T1691,3311)])
    (HERE/('star_successor_probe.json' if __debug__ else 'star_successor_probe.optimized.json')).write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8',newline='\n')
    for x in result['cases']:print(json.dumps(dict(source=x['source_H'],target=x['target'],choices=x['unique_start_star_choices'],hits=[dict(port=h['port'],flips=h['star_flips'],state=h['parent_state'],strong=h['output_strong']) for h in x['hits']])))
