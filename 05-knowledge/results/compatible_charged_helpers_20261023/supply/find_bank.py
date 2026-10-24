"""Deterministic discovery of fixed-cardinality surjective native helper banks."""
from pathlib import Path
from itertools import combinations
import json,random
HERE=Path(__file__).resolve().parent
def encoded(menu):
    a=menu['label'];return a[0],sum(x<<j for j,x in enumerate(a[1:]))
def section(menus,k):
    charges=[encoded(m) for m in menus];out={}
    def visit(start,left,total,parity,mask):
        if left==0:
            out.setdefault(((total&127)<<4)|parity,mask);return
        for i in range(start,len(menus)-left+1):
            x,p=charges[i];visit(i+1,left-1,total+x,parity^p,mask|(1<<i))
    visit(0,k,0,0,0)
    return out
def main():
    allmenus=json.loads((HERE/'uv_discovery.json').read_text())['menus']
    pool=[m for m in allmenus if m['best']['target_pinned_count']==2 and m['best']['common_clocks']]
    pool.sort(key=lambda m:(m['best']['cost'],m['source']))
    rng=random.Random(2026102302);found=[]
    for h,k,trials in [(17,8,100),(18,6,100),(22,5,100),(32,4,60),(len(pool),3,1)]:
        best=0
        for trial in range(trials):
            selected=pool[:h] if trial==0 else rng.sample(pool,h)
            selected.sort(key=lambda m:m['source']);sec=section(selected,k)
            best=max(best,len(sec))
            if len(sec)==2048:
                result=dict(helper_count=h,V_occupancy=k,atomic_support_bound=2+2*min(k,h-k),
                    source_ids=[m['source'] for m in selected],menus=selected,charge_section={str(q):v for q,v in sorted(sec.items())},
                    trial=trial,seed=2026102302,all_subset_count=sum(1 for _ in combinations(range(h),k)),
                    charge_encoding='16*(Xsum mod128) + XOR of parity bits Y1..Y4 in bit positions0..3')
                found.append(result);print('FOUND',h,k,'bound',result['atomic_support_bound'],'trial',trial,'sources',result['source_ids'],flush=True)
                break
        else:print('MISS',h,k,'maxcoverage',best,flush=True)
    (HERE/'banks_discovery.json').write_text(json.dumps(dict(pool_size=len(pool),banks=found),indent=2)+'\n',encoding='utf-8')
if __name__=='__main__':main()
