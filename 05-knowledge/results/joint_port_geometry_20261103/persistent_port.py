"""Three-live-port implementation of C3/C5 retaining the original spare port."""
from collections import Counter
import json
import sympy as s
from joint_probe import HERE,e,pivots,rank,close_labeled

def rename_vector(V,keep):
    ren={old:new for new,old in enumerate(keep)};out=Counter()
    for (sk,m),value in V.items():
        paths=[(ren.get(a,-1),ren.get(b,-1)) for a,b in sk]+[(-1,-1)]*m
        bare=paths.count((-1,-1));ss=tuple(sorted(p for p in paths if p!=(-1,-1)))
        if bare+sum(b==-1 for a,b in ss)<=2:out[ss,bare]+=value
    return out
def birth(V,k,mask,keep):
    out=Counter()
    for (sk,m),value in V.items():
        paths=list(sk)+[(-1,-1)]*m
        ins=[None]+[i for i,(a,b) in enumerate(paths) if b>=0 and mask>>b&1]
        outs=[None]+[i for i,(a,b) in enumerate(paths) if a==-1 or not(mask>>a&1)]
        for i in ins:
            for j in outs:
                if i is not None and i==j:continue
                ps=[x for t,x in enumerate(paths) if t not in (i,j)]
                ps.append((k if i is None else paths[i][0],k if j is None else paths[j][1]))
                bare=ps.count((-1,-1));ss=tuple(sorted(p for p in ps if p!=(-1,-1)))
                out[ss,bare]+=value
    return rename_vector(out,keep)
def macro(j):
    M=s.zeros(e.DIM)
    for row,key in enumerate(e.STATES):
        V=Counter({key:1});V=birth(V,2,1,(0,1,2))
        V=birth(V,3,4,(3,1,2))
        for _ in range(j-1):V=birth(V,3,4,(0,1,2))
        V=birth(V,3,1,(0,1,3));V=birth(V,3,4,(3,1,2))
        V=rename_vector(V,(2,1))
        for target,c in V.items():M[row,e.cap.INDEX[target]]=c
    return M
def main():
    N=e.macro(0);ops={'3p':macro(3),'5p':macro(5),'N':N,**{'J'+str(m):e.dense((0,m,2))*N for m in range(4)}}
    image=N[1:,:];B=image[list(pivots(image.T)),:]
    q=e.named_queries();Q,labels,depth=close_labeled([(k,q[k]) for k in ('h','s','E1','S2')],ops)
    cols=pivots(B*Q);R=Q[:,list(cols)];V=B*R
    reduced={k:V.gauss_jordan_solve(B*M*R)[0] for k,M in ops.items()}
    e.need(all(rank(s.Matrix.vstack(B,B*M))==B.rows for M in ops.values()),'persistent-spare sealed source domain is invariant')
    out={'status':'EXACT_OPERATOR_DISCOVERY','source_rank':B.rows,'observer_rank':V.cols,'depth':depth,
         'source_basis':e.data(B),'query_basis':e.data(R),'query_labels':[labels[i] for i in cols],
         'matrices':{k:e.data(M) for k,M in reduced.items()},'cap_matrices':{k:e.data(M) for k,M in ops.items()}}
    (HERE/'persistent_port.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf8',newline='\n')
    print('persistent spare: source rank',B.rows,'observer',V.cols,'depth',depth,flush=True)
if __name__=='__main__':main()
