"""Rational closure and trace radical for the persistent27 operation algebra."""
import json
import sympy as s
from common import HERE,PREV,need,read,save,data,dm,rref,nullrows

def main():
    packet=read('persistent_port.json');n=27
    S,SP=rref(s.Matrix(packet['source_basis']));QC=S*s.Matrix(packet['query_basis'])
    ops={k:S*s.Matrix(M)[:,list(SP)] for k,M in packet['cap_matrices'].items()}
    need(QC.det()!=0 and all(M*QC==QC*s.Matrix(packet['matrices'][k]) for k,M in ops.items()),'exact conjugacy from physical cap coordinates to previous query coordinates')
    words=read('persistent_certificate.json')['modular_algebra']['word_basis'];matrices=[]
    for word in words:
        M=s.eye(n)
        for letter in word.split('/'):
            if letter:M=M*ops[letter]
        matrices.append(M)
    B=s.Matrix([list(M) for M in matrices]);R,P=rref(B)
    need(R.rows==95,'exact independent word basis95')
    print('rational basis95 pivots',len(P),flush=True)
    RD=dm(R);basis=[s.Matrix(n,n,list(R[i,:])) for i in range(95)];basisD=[dm(M) for M in basis]
    coords={};OD={k:dm(M) for k,M in ops.items()}
    for k,M in OD.items():
        rows=s.Matrix([list((A*M).to_Matrix()) for A in basisD]);C=rows[:,list(P)]
        need(dm(C)*RD==dm(rows),'rational right closure '+k);coords[k]=C
        print('closure',k,flush=True)
    transpose=[j*n+i for i in range(n) for j in range(n)]
    G=(RD*dm(R[:,transpose].T)).to_Matrix();GR,GP=rref(G)
    K,_=rref(nullrows(G));J=(dm(K)*RD).to_Matrix()
    need(K.rows==95-len(GP),'trace kernel dimension')
    print('trace rank',len(GP),'radical candidate',K.rows,flush=True)
    save('algebra_exact.json',{'status':'EXACT_RATIONAL_CLOSURE','dimension':95,'matrix_dimension':27,
      'coordinate_pivots':P,'basis':data(R),'right_generators':{k:data(C) for k,C in coords.items()},
      'source_basis':data(S),'source_pivots':SP,'query_conjugacy':data(QC),'generators':{k:data(M) for k,M in ops.items()},
      'trace_rank':len(GP),'trace_kernel_coordinates':data(K),'trace_kernel_matrices':data(J),
      'trace_gram':data(G),'radical_status':'candidate until nilpotence and semisimple quotient verified'})
if __name__=='__main__':main()
