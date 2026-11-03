"""The terminal-end filtration and its effect algebra on physical sources."""
import json
import sympy as s
from sympy.polys.matrices import DomainMatrix
import numpy as np
from joint_probe import HERE,e,pivots,rank
from algebra_probe import algebra

def rref(M):return DomainMatrix.from_Matrix(M).rref()[0].to_Matrix()
def nullrows(M):return DomainMatrix.from_Matrix(M).nullspace().to_Matrix()
def normalized(M):
    return s.Matrix.vstack(*[M[i,:]/s.gcd_list(list(M[i,:])) for i in range(M.rows)]) if M.rows else M
def main():
    packet=json.loads((HERE/'joint_probe.json').read_text());raw=s.Matrix(packet['physical_source_basis'])
    B=rref(raw);P=pivots(B);need=e.need
    need(B[:,list(P)]==s.eye(31),'coordinate pivots for physical source span')
    levels=[];prefix=[]
    for k in (2,1):
        low=[i for i,(sk,m) in enumerate(e.STATES) if m+sum(b==-1 for a,b in sk)<k]
        C=normalized(nullrows(B[:,low].T));levels.append(C.rows);prefix.append(C)
    stacked=s.Matrix.vstack(*prefix,s.eye(31));S=stacked[list(pivots(stacked.T)),:];Sinv=S.inv()
    need(levels==[14,26] and S.rows==31,'three terminal-end layers14,12,5')
    N=e.macro(0);ops={'3':e.macro(3),'5':e.macro(5),'N':N,**{'J'+str(m):e.dense((0,m,2))*N for m in range(4)}}
    A={k:S*(B*M[:,list(P)])*Sinv for k,M in ops.items()}
    for M in A.values():need(M[:14,14:]==s.zeros(14,17) and M[:26,26:]==s.zeros(26,5),'universal terminal filtration')
    diag={}
    for k,M in A.items():
        D=s.diag(M[:14,:14],M[14:26,14:26],M[26:,26:])
        diag[k]=np.array([[int(x.p)*pow(int(x.q),-1,1009)%1009 for x in row] for row in D.tolist()],dtype=np.int64)
    result=algebra(diag);print('DIAGONAL',result['algebra_rank'],flush=True)
    out={'status':'EXACT_FILTRATION_AND_MODULAR_DIAGONAL','layer_dimensions':[14,12,5],
         'cap_coordinate_pivots':P,'source_basis':e.data(B),'adapted_basis':e.data(S),
         'adapted_matrices':{k:e.data(M) for k,M in A.items()},'diagonal_algebra':result}
    (HERE/'filtration.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf8',newline='\n')
if __name__=='__main__':main()
