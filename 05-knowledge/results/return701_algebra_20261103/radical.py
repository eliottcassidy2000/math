"""Exact nilpotence and the center of the semisimple quotient."""
import json
import sympy as s
from common import HERE,need,read,save,data,dm,rref,nullrows

def main():
    d=json.loads((HERE/'algebra_exact.json').read_text());R=s.Matrix(d['basis']);P=d['coordinate_pivots'];n=27
    K,JP=rref(s.Matrix(d['trace_kernel_coordinates']));F=[i for i in range(95) if i not in JP]
    Q=s.eye(95)[:,F]-s.eye(95)[:,list(JP)]*K[:,F]
    need(K*Q==s.zeros(69,26),'quotient projection annihilates exact trace kernel')
    J=s.Matrix(d['trace_kernel_matrices']);JD=[dm(s.Matrix(n,n,list(J[i,:]))) for i in range(J.rows)]
    W=s.eye(n);dims=[n];flags=[]
    while W.rows:
        WD=dm(W);stack=s.Matrix.vstack(*[(WD*A).to_Matrix() for A in JD]);W,_=rref(stack)
        dims.append(W.rows);flags.append(data(W));print('radical module filtration',dims,flush=True)
        need(len(dims)<=28 and dims[-1]<dims[-2],'strict nilpotent radical filtration')
    ops={k:s.Matrix(v) for k,v in d['generators'].items()};basis=[dm(s.Matrix(n,n,list(R[i,:]))) for i in range(95)]
    rgens={k:s.Matrix(C)[F,:]*Q for k,C in d['right_generators'].items()};lgens={}
    for k,M in ops.items():
        MD=dm(M);rows=s.Matrix([list((MD*basis[i]).to_Matrix()) for i in F]);coords=rows[:,P]
        need(dm(coords)*dm(R)==dm(rows),'left-generator algebra closure')
        lgens[k]=coords*Q
    A=s.Matrix.vstack(*[(rgens[k]-lgens[k]).T for k in ops]);Z,_=rref(nullrows(A))
    print('quotient center dimension',Z.rows,flush=True)
    # Right regular action of each chosen quotient basis representative.
    QB=[]
    for i in F:
        rows=s.Matrix([list((basis[j]*basis[i]).to_Matrix()) for j in F]);QB.append(rows[:,P]*Q)
    centers=[sum((Z[t,j]*QB[j] for j in range(26)),s.zeros(26)) for t in range(Z.rows)]
    x=s.Symbol('x')
    for seed in range(2,30):
        chosen=sum((seed**i*M for i,M in enumerate(centers)),s.zeros(26));cp=s.factor(chosen.charpoly(x).as_expr())
        if len(s.factor_list(cp)[1])==Z.rows:break
    need(len(s.factor_list(cp)[1])==7,'center separates all seven rational components')
    print('central characteristic polynomial',cp,flush=True)
    save('radical.json',{'status':'EXACT_NILPOTENT_TRACE_KERNEL','algebra_dimension':95,'radical_dimension':69,
      'module_radical_filtration_dimensions':dims,'module_radical_filtration_bases':flags,
      'quotient_dimension':26,'quotient_projection':data(Q),'quotient_representative_indices':F,
      'quotient_right_generators':{k:data(M) for k,M in rgens.items()},'quotient_basis_right_actions':[data(M) for M in QB],
      'center_basis_coordinates':data(Z),'center_right_actions':[data(M) for M in centers],
      'chosen_center_action':data(chosen),'center_seed':seed,'center_characteristic_polynomial':str(cp)})
if __name__=='__main__':main()
