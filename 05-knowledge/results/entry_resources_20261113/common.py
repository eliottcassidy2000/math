from pathlib import Path
import sys,json,gzip,hashlib
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
OLD=HERE.parent/'joint_stability_forbidden_20261112'
sys.path.insert(0,str(OLD))
from census_nine import tables,strong,defect
import numpy as np

def need(p,message):
    if not p:raise RuntimeError(message)

def save(name,data):
    tag='' if __debug__ else '.optimized'
    path=HERE/(name+tag+'.json.gz')
    path.write_bytes(gzip.compress((json.dumps(data,separators=(',',':'))+'\n').encode(),mtime=0))
    print(name,'saved',flush=True)

def load(path):
    p=Path(path)
    return json.loads(gzip.decompress(p.read_bytes()) if p.suffix=='.gz' else p.read_text(encoding='utf8'))

def labelled(n,code):
    rr=[0]*n;k=0
    for i in range(n):
        for j in range(i+1,n):
            u,v=(i,j) if code>>k&1 else (j,i)
            rr[u]|=1<<v;k+=1
    return rr

def selected_state(rows,O):
    h,u,w,_=tables(np.array(rows,dtype=np.int64))
    picks=[i for i in range(len(rows)) if O>>i&1]
    return int(h),int(defect(h,w,picks)),sum(int(u[i]) for i in picks)
