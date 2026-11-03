"""Explicit split matrix factors and a faithful finite algebra presentation."""
import json
import sympy as s
from common import HERE,need,save,data,dm,rref

def main():
    a=json.loads((HERE/'algebra_exact.json').read_text());d=json.loads((HERE/'radical.json').read_text())
    R=s.Matrix(a['basis']);ops={k:s.Matrix(M) for k,M in a['generators'].items()};n=27
    flags=[s.Matrix(W) for W in d['module_radical_filtration_bases'][:-1]]
    stacked=s.Matrix.vstack(*reversed(flags),s.eye(n));_,ps=rref(stacked.T);S=stacked[list(ps),:];SI=S.inv()
    sizes=[2,5,14,6];ends=[0,2,7,21,27]
    def diagonal(M):
        A=S*M*SI
        for end in ends[1:-1]:need(A[:end,end:]==s.zeros(end,n-end),'radical module flag invariant')
        return [A[ends[i]:ends[i+1],ends[i]:ends[i+1]] for i in range(4)]
    blockops={k:diagonal(M) for k,M in ops.items()}
    Z=s.Matrix(d['center_basis_coordinates']);F=d['quotient_representative_indices'];seed=d['center_seed']
    z=sum((seed**i*Z[i,:] for i in range(Z.rows)),s.zeros(1,26))
    center=s.Matrix(n,n,list(z*R[F,:]));CZ=diagonal(center)
    x=s.Symbol('x');regular=s.Matrix(d['chosen_center_action']);factors=s.factor_list(regular.charpoly(x).as_expr())[1]
    lambdas=sorted([-s.Poly(f,x).nth(0)/s.Poly(f,x).nth(1) for f,m in factors]);need(len(lambdas)==7,'seven center eigenvalues')
    components=[];extractors=[]
    for lam in lambdas:
        candidates=[];ranks=[]
        for j,C in enumerate(CZ):
            P=s.eye(C.rows)
            for other in lambdas:
                if other!=lam:P=P*(C-other*s.eye(C.rows))/(lam-other)
            need(P*P==P,'central primitive projection on every semisimple layer')
            W,wp=rref(P);ranks.append(W.rows)
            if W.rows:candidates.append((W.rows,j,W,wp))
        dim,j,W,wp=min(candidates,key=lambda q:q[:2]);extractors.append((j,W,wp))
        small={k:(W*blocks[j])[:,list(wp)] for k,blocks in blockops.items()}
        need(all(W*blockops[k][j]==M*W for k,M in small.items()),'simple factor intertwines every generator')
        components.append({'center_eigenvalue':int(lam),'dimension':dim,'ranks_in_layers':ranks,'chosen_layer':j,
                           'layer_basis':data(W),'layer_pivots':wp,'generator_matrices':{k:data(M) for k,M in small.items()}})
    image=[]
    for i in range(95):
        blocks=diagonal(s.Matrix(n,n,list(R[i,:])));row=[]
        for j,W,wp in extractors:row.extend(list((W*blocks[j])[:,list(wp)]))
        image.append(row)
    Phi=s.Matrix(image);_,pp=rref(Phi);K=s.Matrix(a['trace_kernel_coordinates'])
    need(Phi.cols==26==len(pp) and K*Phi==s.zeros(69,26),'full split semisimple quotient with exactly the computed kernel')
    need(sorted(c['dimension'] for c in components)==[1,1,1,1,2,3,3],'complete rational matrix-factor list')
    save('semisimple.json',{'status':'EXACT_RATIONAL_ALGEBRA_CLASSIFICATION','dimension':95,'radical_dimension':69,
      'radical_nilpotence_index':4,'simple_dimensions':[c['dimension'] for c in components],
      'radical_layer_dimensions':sizes,'adapted_basis':data(S),'components':components,'quotient_map':data(Phi),
      'presentation':'The95 explicit basis matrices and seven right-generator actions in algebra_exact.json completely specify the rational algebra; generators span it by the listed inherited words.'})
    print('PASS A/J = Q^4 + M2(Q) + M3(Q)^2, J^4=0 != J^3; A dimension95',flush=True)
if __name__=='__main__':main()
