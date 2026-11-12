"""Search physical endpoint interfaces for both heads of the crossing ray.

Targets are -481+5841*2^(k+1), and9*(-107+649*2^(k+2)).
An interface of defectD and slopeA*2^j certifies the corresponding
affine ray only after its positive tail length becomes available.
"""
from pathlib import Path
from collections import Counter
import gzip,csv,json,time
import numpy as np
from numba import njit
from census_nine import read,save,need,tables,defect
from two_block_probe import decode
HERE=Path(__file__).resolve().parent

@njit(cache=True)
def matches(H,U,W):
    n=len(U);N=1<<n;slopes=np.zeros(N,dtype=np.int64);taxes=np.zeros(N,dtype=np.int64);found=[]
    for m in range(1,N):
        b=m&-m;v=0
        while (1<<v)!=b:v+=1
        rest=m^b;slopes[m]=slopes[rest]+U[v];tax=taxes[rest]
        for w in range(v+1,n):
            if rest>>w&1:tax+=W[v,w]
        taxes[m]=tax;d=H-tax;u=slopes[m]
        if d not in [-481,-107]:continue
        odd=u;j=0
        while odd>0 and odd%2==0:odd//=2;j+=1
        if (d==-481 and odd==5841) or (d==-107 and odd==649):found.append((m,d,u,j))
    return found

def graphs():
    with gzip.open(HERE.parent/'metagraph_aspect_sequences_20261111/records.jsonl.gz','rt') as f:
        for line in f:yield json.loads(line)['rows']
    with gzip.open(HERE.parent/'blue_quotient_20261020b/n9_classes.tsv.gz','rt') as f:
        for r in csv.DictReader(f,delimiter='\t'):yield decode(int(r['code']),9)
    for r in read('unit_ten.json.gz')['records']:yield r['rows']

def run():
    out=[];stats=Counter();start=time.time()
    for rows in graphs():
        n=len(rows);stats[n]+=1;h,u,w,_=tables(np.array(rows,dtype=np.int64))
        for mask,d,slope,j in matches(h,u,w):
            need(defect(h,w,[v for v in range(n) if mask>>v&1])==d,'literal target interface')
            out.append({'rows':rows,'selected':mask,'H':int(h),'defect':d,'U':slope,'dyadic_shift':j})
        if sum(stats.values())%50000==0:print(json.dumps({'sources':dict(stats),'interfaces':len(out),'seconds':round(time.time()-start,1)}),flush=True)
    need(stats==Counter({1:1,2:1,3:2,4:4,5:12,6:56,7:456,8:6880,9:191536,10:465023}),'complete declared source universe')
    result={'scope':'All ordinary tournaments through9 and every strong ordinary order10 source admitting at least one unit selection; all selected sets tested for two specified nonunit endpoint interface types. Not every order10 tournament.','source_counts':dict(stats),'interfaces':len(out),'records':out}
    save('head_repair.json.gz',result);print(json.dumps({k:v for k,v in result.items() if k!='records'}),flush=True)
if __name__=='__main__':run()
