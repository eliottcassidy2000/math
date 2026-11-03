"""Exact rational splitting of the full shared-port return effect algebra."""
import json
import sympy as s
from sympy.polys.matrices import DomainMatrix
from return_probe import HERE,need
import effects as e

def sylvester(blocks):
    equations=[];rhs=[]
    for D,E,X in blocks:
        equations.append(s.kronecker_product(E,s.eye(D.rows))-s.kronecker_product(s.eye(E.rows),D.T))
        rhs.extend(list(X))
    A=s.Matrix.vstack(*equations);aug=A.row_join(s.Matrix(rhs))
    R,piv=DomainMatrix.from_Matrix(aug).rref();R=R.to_Matrix();n=A.cols
    need(tuple(piv)==tuple(range(n)),'unique exact common splitting')
    solution=R[:n,n];K=s.Matrix(blocks[0][1].rows,blocks[0][0].rows,list(solution))
    need(all(E*K-K*D==X for D,E,X in blocks),'every generator splits exactly')
    return K
def main():
    d=json.loads((HERE/'filtration.json').read_text());ops={k:s.Matrix(v) for k,v in d['adapted_matrices'].items()}
    K=sylvester([(M[:14,:14],M[14:,14:],M[14:,:14]) for M in ops.values()])
    T=s.eye(31);T[14:,:14]=K;U={k:T*M*T.inv() for k,M in ops.items()}
    L=sylvester([(M[14:26,14:26],M[26:,26:],M[26:,14:26]) for M in U.values()])
    V=s.eye(31);V[26:,14:26]=L;G=V*T;inv=G.inv();Z={k:G*M*inv for k,M in ops.items()}
    for M in Z.values():need(M==s.diag(M[:14,:14],M[14:26,14:26],M[26:,26:]),'full exact block diagonalization')
    need(d['diagonal_algebra']['algebra_rank']==14**2+12**2+5**2,'modular matrix-word certificate fills the rational diagonal upper bound')
    out={'status':'EXACT_RATIONAL_ALGEBRA_CERTIFICATE','simple_block_dimensions':[14,12,5],
         'algebra_dimension':365,'first_splitting':e.data(K),'second_splitting':e.data(L),
         'diagonalizer':e.data(G),'diagonal_matrices':{k:e.data(M) for k,M in Z.items()}}
    (HERE/'split_algebra.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf8',newline='\n')
    print('PASS: exact shared-port algebra M14(Q)+M12(Q)+M5(Q), eight invariant quotient types',flush=True)
if __name__=='__main__':main()
