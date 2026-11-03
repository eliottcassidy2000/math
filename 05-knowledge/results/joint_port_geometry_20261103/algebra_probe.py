"""Exact finite-field word-algebra ranks for enlarged shared-port observers."""
import json
import numpy as np
import sympy as s
from return_probe import HERE,need

def algebra(ops,p=1009):
    n=next(iter(ops.values())).shape[0];d=n*n;E=np.zeros((d,d),dtype=np.int64);piv=[];basis=[];queue=[('',np.eye(n,dtype=np.int64))];head=0
    while head<len(queue):
        word,M=queue[head];head+=1;v=M.reshape(d).copy()
        r=len(piv)
        if r:v=(v-v[piv]@E[:r])%p
        nz=np.flatnonzero(v)
        if not len(nz):continue
        pivot=int(nz[0]);v=v*pow(int(v[pivot]),-1,p)%p
        if r:E[:r]=(E[:r]-E[:r,pivot,None]*v[None,:])%p
        E[r]=v;piv.append(pivot);basis.append(word)
        queue.extend((word+'/'+name,(M@A)%p) for name,A in ops.items())
        if len(piv)%50==0:print('algebra rank',len(piv),'dimension',n,flush=True)
    return {'prime':p,'dimension':n,'algebra_rank':len(piv),'word_basis':basis,'coordinate_pivots':piv}
def main():
    d=json.loads((HERE/'joint_probe.json').read_text());out=[]
    for c in d['cases']:
        if c['masks'] not in ([3],[0,1,2,3]):continue
        p=1009;ops={name:np.array([[int(s.Rational(x).p)*pow(int(s.Rational(x).q),-1,p)%p for x in row] for row in rows],dtype=np.int64) for name,rows in c['matrices'].items()}
        result=algebra(ops,p);result['masks']=c['masks'];out.append(result);print('RESULT',c['masks'],result['dimension'],result['algebra_rank'],flush=True)
    (HERE/'algebra_probe.json').write_text(json.dumps({'status':'FINITE_EXACT_MODULAR_ALGEBRA','cases':out},indent=2)+'\n',encoding='utf8',newline='\n')
if __name__=='__main__':main()
