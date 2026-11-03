"""Finite exact algebra and phase discovery for the joint return observer."""
import json
from itertools import combinations
import sympy as s
from effects import HERE,need,data

def flat(M):return s.Matrix(list(M))
def reduce_basis(pairs):
    Q=s.Matrix.hstack(*[flat(M) for word,M in pairs]);piv=Q.rref()[1]
    return [pairs[i] for i in piv]
def power(M,n,m):
    out=s.eye(M.rows)
    while n:
        if n&1:out=(out*M).applyfunc(lambda x:x%m)
        M=(M*M).applyfunc(lambda x:x%m);n//=2
    return out

def mod_rank(M,p):
    a=[[int(x)%p for x in row] for row in M.tolist()];row=0
    for col in range(M.cols):
        pivot=next((i for i in range(row,M.rows) if a[i][col]),None)
        if pivot is None:continue
        a[row],a[pivot]=a[pivot],a[row];c=pow(a[row][col],-1,p);a[row]=[(x*c)%p for x in a[row]]
        for i in range(M.rows):
            if i!=row:
                c=a[i][col];a[i]=[(x-c*y)%p for x,y in zip(a[i],a[row])]
        row+=1
        if row==M.rows:break
    return row

def main():
    d=json.loads((HERE/'effects.json').read_text());ops={k:s.Matrix(v) for k,v in d['natural_matrices'].items()}
    B=[('',s.eye(6))]
    while True:
        C=reduce_basis(B+[(word+k,M*A) for word,M in B for k,A in ops.items()])
        if len(C)==len(B):break
        B=C
    need(len(B)==14,'full direct sum M2+Q+M3 algebra')
    witnesses=s.Matrix.hstack(*[flat(M) for word,M in B]);units=[]
    for block in ((0,1),(2,),(3,4,5)):
        for i in block:
            for j in block:
                E=s.zeros(6);E[i,j]=1;c=witnesses.gauss_jordan_solve(flat(E))[0]
                need(sum((v*M for v,(_,M) in zip(c,B)),s.zeros(6))==E,'constructive block matrix unit')
                units.append({'i':i,'j':j,'coefficients':[str(x) for x in c]})
    periods=[]
    for k,n in (('3',12),('5',48)):
        A=ops[k]
        need(power(A,n,7)==s.eye(6) and all(power(A,j,7)!=s.eye(6) for j in range(1,n)),'exact joint phase period')
        lead=(power(A,n,49)-s.eye(6))/7
        need(all(x.q==1 for x in lead) and any(x%7 for x in lead),'nonzero first 7-adic lift')
        periods.append({'letter':k,'order_mod7':n,'first_lift_mod7':data(lead.applyfunc(lambda x:x%7))})
    cubes=[]
    for mask in range(8):
        selected=[k for i,k in enumerate(('3','5','N')) if mask>>i&1]
        # Existing count sector remains closed; mixed words determine all observations.
        basis=[s.eye(6)[:,i] for i in (0,1,2,3)]
        while True:
            Q=s.Matrix.hstack(*basis,*[ops[k]*q for k in selected for q in basis]);p=Q.rref()[1]
            nxt=[Q[:,i] for i in p]
            if len(nxt)==len(basis):break
            basis=nxt
        cubes.append({'operations':selected,'future_joint_dimension':len(nxt)})
    rank7=mod_rank(witnesses,7);need(rank7==14,'same full block algebra over F7')
    block=ops['3'][3:6,3:6];I=s.eye(3)
    line=(2*(block*block+2*I)).applyfunc(lambda x:x%7);plane=(I-line).applyfunc(lambda x:x%7)
    need((line*line-line).applyfunc(lambda x:x%7)==s.zeros(3) and mod_rank(line,7)==1,'spectral projector onto eigenvalue3')
    J=power(block,3,7)
    need((plane*(J*J+I)).applyfunc(lambda x:x%7)==s.zeros(3),'C3 cubed is a quarter-turn on the quadratic plane')
    v=next(plane[i,:] for i in range(3) if any(plane[i,:]));w=(v*J).applyfunc(lambda x:x%7)
    intertwiner=s.Matrix.vstack(v,w)
    need(mod_rank(intertwiner,7)==2 and (intertwiner*J-s.Matrix([[0,1],[-1,0]])*intertwiner).applyfunc(lambda x:x%7)==s.zeros(2,3),'explicit C4 intertwiner')
    out={'status':'FINITE_EXACT_ALGEBRA_CERTIFICATE','dimension':len(B),'dimension_mod7':rank7,'algebra_words':[w for w,M in B],
         'matrix_units':units,'phase_periods':periods,'operation_cube':cubes,
         'characteristic_polynomials':{k:str(s.factor(M.charpoly().as_expr())) for k,M in ops.items()},
         'quarter_turn':{'scalar_projector':data(line),'plane_projector':data(plane),'J':data(J),'intertwiner':data(intertwiner)}}
    (HERE/'algebra.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf8',newline='\n')
    print('PASS: full14-dimensional effect algebra, joint periods12/48 and operation cube',cubes)
if __name__=='__main__':main()
