"""Look past an outer-tail clock: all unique-start parents at bounded order."""
from pathlib import Path
import sys,time,json,gzip
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'berggren_target_return_20261111'/'reset'))
import funnel_six_state as f

def minimum_strong(n):
    q,r=divmod(n,3)
    return 3*5**(q-1) if r==0 else (5**q if r==1 else 9*5**(q-1))

def run(N=16):
    targets={405*(1<<k)+71 for k in range(1,7)}|{881,2645,7937};levels=[];hits=[];t0=time.time()
    for n in range(3,N+1):
        local=[]
        for rank in range(f.a.class_counts(n)[n]):
            state=f.compile_rank(n,rank);h=state[0]+state[4]-1
            if h in targets:local.append(dict(order=n,rank=rank,state=state,H=h))
        levels.append(dict(order=n,classes=f.a.class_counts(n)[n],hits=len(local)))
        hits+=local;print(n,len(local),'total hits',len(hits),'seconds',round(time.time()-t0,2),flush=True)
    tag='' if __debug__ else '.optimized'
    bounds=[]
    for h in sorted(targets):
        n=3
        while minimum_strong(n)<=(h+1)//2:n+=1
        bounds.append(dict(H=h,parent_H_bound=(h+1)//2,max_parent_order=n-1,excluded_next_order_minimum=minimum_strong(n),complete=N>=n-1,strong_ready_classes=sum(r['H']==h for r in hits)))
    (HERE/('next_target'+tag+'.json.gz')).write_bytes(gzip.compress(json.dumps(dict(scope=f'All strong unique-start rooted parents through{N}; exact ready child count; all-order status determined per target by the cited Busch bound',targets=sorted(targets),bounds=bounds,levels=levels,hits=hits),separators=(',',':')).encode(),mtime=0))
if __name__=='__main__':run(int(sys.argv[1]) if len(sys.argv)>1 else 16)
