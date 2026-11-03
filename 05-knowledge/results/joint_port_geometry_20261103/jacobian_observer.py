"""Exact mixed-response recovery and the boundary term of finite Weyl truncation."""
import json
import sympy as s
from return_probe import HERE,need

def main():
    x,y=s.symbols('x y');examples=[]
    for P,B in ((x,y),(x+y**2,y),(x+(y+x**2)**2,y+x**2)):
        D=lambda f:s.expand(s.diff(P,x)*s.diff(f,y)-s.diff(P,y)*s.diff(f,x))
        k=D(B);need(k==1,'actual supplied polynomial mate')
        Q=x**3+2*x*y+y**4;R1=D(Q);R2=D(B*Q)
        need(s.expand(R2-B*R1)==Q,'mixed response exactly reconstructs the source')
        need(s.expand(D(R2)-B*D(R1)-2*R1)==0,'exact image compatibility')
        U=x*y+1;V=y**3-x;source=s.expand(V-B*U);defect=s.expand(D(V)-B*D(U)-2*U)
        need(s.expand(D(source)-defect-U)==0 and s.expand(B*(D(source)-defect)+source-V)==0,'split source/defect normal form')
        # A nonunit mixed multiplier still separates, but carries divisibility.
        f=B**2;a=D(f);N=s.expand(D(f*Q)-f*D(Q))
        need(s.expand(N-a*Q)==0,'nonunit image divisibility is retained')
        examples.append({'P':str(P),'B':str(B),'Q':str(Q),'D_Q':str(R1),'D_BQ':str(R2),'nonunit_coefficient':str(a)})
    truncations=[]
    for d in range(1,13):
        D=s.zeros(d+1);M=s.zeros(d+1)
        for j in range(1,d+1):D[j-1,j]=j
        for j in range(d):M[j+1,j]=1
        E=s.zeros(d+1);E[d,d]=1
        need(D*M-M*D==s.eye(d+1)-(d+1)*E,'exact finite-dimensional boundary commutator')
        truncations.append({'degree_cap':d,'boundary_rank':1,'boundary_coefficient':d+1})
    (HERE/'jacobian_observer.json').write_text(json.dumps({'status':'FINITE_EXACT_SYMBOLIC_CERTIFICATE','examples':examples,
         'truncation_controls':truncations,'scope':'mate supplied; this does not prove mate existence or polynomial invertibility'},indent=2)+'\n',encoding='utf8',newline='\n')
    print('PASS: exact mixed-response reconstruction, image defect and finite-rank boundary corrections',flush=True)
if __name__=='__main__':main()
