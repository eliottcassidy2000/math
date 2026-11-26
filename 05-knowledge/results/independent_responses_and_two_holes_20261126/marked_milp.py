from pathlib import Path
import sys
import numpy as np
from scipy.optimize import milp, LinearConstraint, Bounds
from scipy.sparse import coo_matrix

def parse(path):
    it = iter(map(int,Path(path).read_text().split()))
    n=next(it); rows=[]
    for h in range(n):
        a,k=next(it),next(it)
        rows.append((a,[tuple(next(it) for _ in range(4)) for j in range(k)]))
    assert next(it,None) is None and n==53
    mark=rows[0][1][0][1:]
    return [(h,r,y,ss) for h,(r,row) in enumerate(rows[1:]) for y,*ss in row
            if all(z!=w for z,w in zip(mark,ss))]

def run(path, cap, secs):
    b=parse(path); nn=len(b); assert nn==1741
    cov=[]; cols=[]; vals=[]; cost=[]; m=52+3*56
    bounds1=np.zeros(m+1); bounds2=np.ones(m+1)
    bounds1[:52]=1; bounds1[52:]=-np.inf; bounds2[m]=cap
    for j,(a,r,y,ss) in enumerate(b):
        for i in [a]+[52+56*k+s for k,s in enumerate(ss)]:
            cov.append(i); cols.append(j); vals.append(1)
        cost.append(int(r>=0 and r!=y))
    for k in range(168):
        cov.extend([52+k,m]); cols.extend([nn+k,nn+k]); vals.extend([-1,1]); cost.append(0)
    mat=coo_matrix((vals,(cov,cols)),shape=(m+1,nn+168)).tocsc()
    res=milp(np.array(cost,dtype=float),integrality=np.ones(nn+168),
             bounds=Bounds(0,1),constraints=LinearConstraint(mat,bounds1,bounds2),
             options={"time_limit":secs,"mip_rel_gap":0})
    print({"status":res.status,"message":res.message,"cost":res.fun,
           "nodes":getattr(res,"mip_node_count",None),"bound":getattr(res,"mip_dual_bound",None)},flush=True)
    if res.x is not None:
        pick=[b[j] for j in range(nn) if res.x[j]>0.5]
        assert len(pick)==52 and len({a for a,r,y,s in pick})==52
        z=np.zeros((3,56),dtype=int)
        for a,r,y,S in pick:
            for k,s in enumerate(S):z[k,s]+=1
        q=(z*(z-1)//2).sum(axis=1);ed=sum(r>=0 and r!=y for a,r,y,S in pick)
        assert q.sum()<=cap and abs(ed-res.fun)<1e-6
        print({"square_profile":q.tolist(),"edits":ed,"choices":[(a,y) for a,r,y,S in pick]},flush=True)

if __name__=="__main__":
    run(sys.argv[1],int(sys.argv[2]),float(sys.argv[3]))
