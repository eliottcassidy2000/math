"""One targeted next-order joint core/tail inverse, with actual graph output."""
from pathlib import Path
import sys,json,time,argparse
HERE=Path(__file__).resolve().parent
RESULTS=HERE.parent.parent
sys.path.insert(0,str(RESULTS/'changing_core_returns_20261113'/'returns'))
import inverse as old

def main(order):
    target=53084231;hits=[];t0=time.time();nodes=0;tested=0
    for rank in range(old.f.a.class_counts(order)[order]):
        state=old.f.compile_rank(order,rank);h,e,s,p,f,t=state;G=s+p
        ns,meta=old.inverse(target,G,h,e,G-e-f-1);nodes+=meta['nodes'];tested+=1
        if ns:
            row=dict(order=order,rank=rank,state=state,words=old.address(order,rank));n=min(ns)
            T,y,V,x,w=old.physical(row,n);new=old.f.advance(state,w)
            old.need(new[0]+new[4]-1==target,'actual full target compiler')
            old.need(old.f.a.strong(T),'actual source strong')
            old.need([v for v in range(len(T)) if T[v]>>y&1]==[x],'actual indegree-one output')
            core=old.literal_core(row)
            old.need(old.direct(core)==h and old.direct(core,0)==1,'independent literal core H/start')
            hits.append(dict(core=row,tail=n,word=''.join(map(str,w)),rows=T,root=y,predecessor=x,parent_state=new,state=(target,new[0]-1,new[0],target+1)))
            print('HIT',order,rank,n,'seconds',round(time.time()-t0,2),flush=True);break
        if rank and rank%500000==0:print('tested',rank,'seconds',round(time.time()-t0,2),flush=True)
    rec=dict(status='PASS',target=target,core_order=order,rooted_classes_tested=tested,complete_at_this_order=not hits,total_inverse_nodes=nodes,arbitrary_outer_tail=True,hits=hits)
    suffix='' if __debug__ else '.optimized'
    (HERE/('joint_probe'+suffix+'.json')).write_text(json.dumps(rec,indent=2)+'\n',encoding='utf-8',newline='\n')
    print('PASS target core-order probe',order,'hits',len(hits),'seconds',round(time.time()-t0,2),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--order',type=int,default=18);main(p.parse_args().order)
