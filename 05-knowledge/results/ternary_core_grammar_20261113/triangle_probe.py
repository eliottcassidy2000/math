"""All cross matrices for triangle/strong-core amalgamation, at bounded core order."""
from pathlib import Path
import sys,json,gzip,time
import numpy as np
from numba import njit
HERE=Path(__file__).resolve().parent

def need(ok,msg):
    if not ok:raise RuntimeError(msg)

def save(name,data):
    tag='' if __debug__ else '.optimized'
    (HERE/(name+tag+'.json.gz')).write_bytes(gzip.compress(json.dumps(data,separators=(',',':')).encode(),mtime=0))

def bank(limit):
    p=HERE.parent/'multicore_metagraph_20261113/rooted_growth.json.gz'
    data=json.loads(gzip.decompress(p.read_bytes()))
    # The old final vertex is the named shared root. A strong tournament has
    # a Hamilton path to every vertex, so the inherited end bank includes all.
    from importlib.util import spec_from_file_location,module_from_spec
    path=HERE.parent/'multicore_metagraph_20261113/rooted_growth.py'
    spec=spec_from_file_location('inherited_growth',path);m=module_from_spec(spec);spec.loader.exec_module(m)
    out=[]
    for row in data['records']:
        n=row['n']
        if 3<=n<=limit and all(m.endpoints(tuple(row['rows']))):out.append(row)
    return out

@njit
def ends(T):
    n=len(T);d=np.zeros((1<<n,n),dtype=np.int64)
    for v in range(n):d[1<<v,v]=1
    for mask in range(1,1<<n):
        for v in range(n):
            if not(mask>>v&1):continue
            rest=mask^(1<<v)
            if rest:
                for u in range(n):
                    if T[u]>>v&1:d[mask,v]+=d[rest,u]
    return d[-1]

@njit
def attach(T,r,mask):
    n=len(T);U=np.zeros(n+2,dtype=np.int64);U[:n]=T
    # r -> a -> b -> r; bits encode a/b -> old nonroot vertices.
    a=n;b=n+1;U[r]|=1<<a;U[a]|=1<<b;U[b]|=1<<r;k=0
    for x in (a,b):
        for v in range(n):
            if v==r:continue
            if mask>>k&1:U[x]|=1<<v
            else:U[v]|=1<<x
            k+=1
    return U

@njit
def scan(T,r):
    total=1<<(2*(len(T)-1));values=np.zeros(total,dtype=np.int64)
    for mask in range(total):values[mask]=ends(attach(T,r,mask)).sum()
    return values

def run(limit=6):
    records=[];t0=time.time();roots=bank(limit)
    for i,rec in enumerate(roots):
        T=np.array(rec['rows'],dtype=np.int64);h=int(ends(T).sum());values=scan(T,len(T)-1)
        need(int(min(values))>=3*h,'Busch merge lower bound')
        counts={int(v):int(c) for v,c in zip(*np.unique(values,return_counts=True))}
        hits={str(d):next((j for j,x in enumerate(values) if x==3*h+d),None) for d in (0,2,4)}
        records.append(dict(n=len(T),rows=list(map(int,T)),core_H=h,ordinary_code=rec['ordinary_code'],minimum=int(min(values)),exits=hits,spectrum=counts))
        if i%25==0:print(i+1,'/',len(roots),'seconds',round(time.time()-t0,2),flush=True)
    groups={}
    for rec in records:
        key=(rec['n'],rec['ordinary_code'])
        g=groups.setdefault(key,dict(n=rec['n'],code=rec['ordinary_code'],H=rec['core_H'],minimum=rec['minimum'],exits=set()))
        g['minimum']=min(g['minimum'],rec['minimum'])
        g['exits'].update(int(d) for d,mask in rec['exits'].items() if mask is not None)
    classes=[]
    for g in groups.values():g['exits']=sorted(g['exits']);classes.append(g)
    failed=[g for g in classes if g['exits']!=[0,2,4]]
    result=dict(scope=f'Every rooted strong class of orders3..{limit}, every triangle cross orientation; source bank inherited and pinned',rooted_inputs=len(records),cross_matrices=sum(1<<(2*(r['n']-1)) for r in records),records=records,classes=classes,failed_classes=failed)
    save('triangle_probe',result)
    print('COUNTS',len(records),len(classes),'failed',len(failed),flush=True)
    print('FIRST FAILURES',failed[:12],flush=True)

if __name__=='__main__':run(int(sys.argv[1]) if len(sys.argv)>1 else 6)
