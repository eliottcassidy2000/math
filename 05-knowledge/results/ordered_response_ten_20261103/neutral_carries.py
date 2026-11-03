"""All-carry formula, cover dominance, and actual source-preserving native programs."""
from pathlib import Path
from math import factorial,comb
import sys,json,sympy as s
from bank import subset_h,path_poly,strong,gate,op,need
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
PRIOR=ROOT/'05-knowledge/results/apex_wedge_constructor_20261103'
sys.path.insert(0,str(PRIOR))
from mechanisms import chart,source,gap_h,reset
from profile import response
from controllers import native
r=s.Symbol('r')

def carry(P):
    n=P+1;j=0
    while n%2==0:j+=1;n//=2
    return j

def formula(P):
    k=P.bit_length();j=carry(P)
    a=k+1 if P==(1<<k)-1 else 2*k-j
    return s.expand(-r**a*(r*r-1)*((r+1)**(j-1)-(r-1)**(j-1)))

def falling_coeffs(P):
    if P==0:return [0]
    deg=s.degree(P,r);v=[int(P.subs(r,i)) for i in range(deg+1)]
    return [sum((-1)**(k-i)*comb(k,i)*v[i] for i in range(k+1))//factorial(k) for k in range(deg+1)]

def replay(program):
    A=tuple(program['source']);n=len(A);mask=program['birth_out']
    B=tuple(row|((1<<n) if not(mask>>i&1) else 0) for i,row in enumerate(A))+(mask,)
    for P in program['paths']:B=native(B,tuple(P))
    endpoint=tuple(row&((1<<n)-1) for row in B[:n])
    need(endpoint==tuple(program['target']),'independent native replay with named ports')

def main():
    formulas=[]
    for P in range(1,256,2):
        w,u=chart(P,P+1);D=formula(P)
        need(s.expand(response(w,u,r)-D)==0,'uniform binary carry response')
        need(gap_h(w,u)==gap_h(w,u,False),'scalar H-neutrality')
        need((D==0)==(carry(P)==1),'complete profile-neutral classification')
        coeff=falling_coeffs(-D)
        need(all(v>=0 for v in coeff) and (D==0 or coeff[2]>0),'path-cover dominance')
        need(s.expand(sum(v*s.prod(r-j for j in range(k)) for k,v in enumerate(coeff))+D)==0,'exact falling basis reconstruction')
        program=reset(w,u);replay(program)
        formulas.append({'P':P,'R':P+1,'carry':carry(P),'order':len(program['source']),
                         'H':gap_h(w,u),'F_plus_minus':str(s.factor(D)),'cover_difference':coeff})
    examples=[]
    for P in (1,3,5,7,11):
        w,u=chart(P,P+1);plus=source(w,u);minus=source(w,u,False);D=formula(P)
        HP=subset_h(plus);HM=subset_h(minus)
        qp=sum(v*HP[(len(HP)-1)^i] for i,v in enumerate(HP));qm=sum(v*HM[(len(HM)-1)^i] for i,v in enumerate(HM))
        need(strong(plus) and strong(minus) and HP[-1]==HM[-1]==gap_h(w,u),'literal strong neutral endpoints')
        need(qp-qm==D.subs(r,2),'ordered two-bin convolution')
        record={'P':P,'program':reset(w,u),'H':HP[-1],'F2_plus':qp,'F2_minus':qm}
        if P in (1,3,5):
            FP,_,cp=path_poly(plus);FM,_,cm=path_poly(minus)
            need(s.expand(FP-FM-D)==0 and all(b>=a for a,b in zip(cp,cm)),'all literal cover counts')
            record.update({'profile_plus':str(FP),'profile_minus':str(FM)})
        if P in (1,3,7):
            after_plus=subset_h(gate(plus))[-1];after_minus=subset_h(gate(minus))[-1]
            need(after_minus-after_plus==-D.subs(r,2),'literal triangle attachment improvement')
            record.update({'G_plus_H':after_plus,'G_minus_H':after_minus})
        examples.append(record)
    # Positive profile-operator propagation, with a nonzero scalar effect after every further gate.
    propagation=[];g=-formula(3)
    for depth in range(1,6):
        g=op(g);coeff=falling_coeffs(g)
        need(all(v>=0 for v in coeff) and g.subs(r,1)>0,'iterated triangle dominance')
        propagation.append({'depth':depth,'H_improvement':int(g.subs(r,1)),'profile_difference':str(g)})
    out={'status':'PROVED_SCOPED','odd_P_universe':[1,255],'symbolic_and_native_checks':len(formulas),
         'formulas':formulas,'literal_examples':examples,'iterated_G_improvements':propagation}
    (HERE/'neutral_carries.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf8',newline='\n')
    print('PASS: 128 carry formulas and named native resets; H497 neutral repair gains48 after G')
if __name__=='__main__':main()
