"""Optional numerical metric discovery. Not required by the exact proof checker.
Requires numpy/scipy; output is an untrusted candidate until rationally checked.
"""
import argparse,json
import numpy as np
from scipy.optimize import linprog
from scipy.sparse import csr_matrix
from verify import active,block,mv,R

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--max-length',type=int,default=7);ap.add_argument('--rho2',type=float,default=.995);ap.add_argument('--output',default='candidate_metric.json');args=ap.parse_args()
    words=[w for w in active(args.max_length) if sum(a!=b for a,b in zip(w,w[1:]))<=3]
    mats=list(map(block,words));m=len(words);transitions=[]
    for u,U in enumerate(mats):
        for v,V in enumerate(mats):
            q=V[0][1];t=V[1][1]
            for r in R:
                a,w,d=mv(V,mv(U,r));rr=w/a;s=d/a
                transitions.append((u,v,np.array([[t-q*rr,rr],[-q*s,s]])))
    def form(z):return np.array([z[0]**2,2*z[0]*z[1],z[1]**2])
    rows=[];cols=[];data=[];rhs=[]
    def addrow(indices,coefficients,b):
        i=len(rhs);rows.extend([i]*len(indices));cols.extend(indices);data.extend(coefficients);rhs.append(b)
    def add(u,v,K,z):addrow(list(range(3*u,3*u+3))+list(range(3*v,3*v+3)),list(-args.rho2*form(z))+list(form(K@z)),-1e-6)
    for u,v,K in transitions:
        for z in [np.array([1,0]),np.array([0,1]),np.array([1,1]),np.array([1,-1])]:add(u,v,K,z)
    for step in range(100):
        C=csr_matrix((data,(rows,cols)),shape=(len(rhs),3*m))
        res=linprog(np.tile([1,0,1],m),A_ub=C,b_ub=rhs,bounds=[(1,None),(None,None),(1,None)]*m,method='highs')
        if not res.success:
            print('No feasible metric found in this numerical relaxation.',res.message);return
        Qs=[np.array([[a,b],[b,c]]) for a,b,c in res.x.reshape(m,3)];bad=0
        for u,v,K in transitions:
            values,vectors=np.linalg.eigh(args.rho2*Qs[u]-K.T@Qs[v]@K)
            if values[0]<-1e-8:add(u,v,K,vectors[:,0]);bad+=1
        for u,Q in enumerate(Qs):
            values,vectors=np.linalg.eigh(Q)
            if values[0]<1e-4:addrow(list(range(3*u,3*u+3)),list(-form(vectors[:,0])),-1);bad+=1
        print(step,bad,flush=True)
        if bad==0:
            with open(args.output,'w') as f:json.dump({'status':'NUMERICAL CANDIDATE, NOT CERTIFIED','words':words,'Q':res.x.tolist(),'rho2':args.rho2},f,indent=2)
            print('Candidate saved; rational rounding and exact PSD replay are still required.');return
    print('Iteration budget exhausted; no certificate claim.')
if __name__=='__main__':main()
