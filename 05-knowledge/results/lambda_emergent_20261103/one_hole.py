"""Named gap automata: all-order one-hole ready counts and full golden data."""
from itertools import product
import json
import sympy as s
from effects import HERE,old,junction,need,data
from controls import values

def model(bins,forbidden=()):
    gaps=[];fixed=[]
    for k,ports in enumerate(bins):
        for j in range(len(ports)+1):
            pair=(ports[j-1] if j else None,ports[j] if j<len(ports) else None)
            (fixed if (k,j) in forbidden else gaps).append(pair)
    states=list(product(range(3),repeat=len(gaps)));idx={v:i for i,v in enumerate(states)}
    def transfer(bits):
        A=s.zeros(len(states))
        for i,state in enumerate(states):
            for j,(left,right) in enumerate(gaps):
                if state[j]==0 and left is not None and not bits[left]:continue
                new=list(state);new[j]=2 if right is None or not bits[right] else 1
                A[i,idx[tuple(new)]]+=1
        return A
    def empty_ok(left,right):return left is None or right is None or (left,right)==(0,1)
    final=s.Matrix([int(all(x==2 if x else empty_ok(*pair) for x,pair in zip(state,gaps)) and all(empty_ok(*pair) for pair in fixed)) for state in states])
    initial=s.zeros(1,len(states));initial[0,idx[(0,)*len(gaps)]]=1
    A=transfer((1,1));B=transfer((1,0));C=transfer((0,0))
    residual=s.eye(A.rows)
    for j in range(len(gaps)+1):residual=residual*(A-j*s.eye(A.rows))
    need(residual==s.zeros(A.rows),'universal commuting-idempotent annihilator')
    return initial,A,B,C,final,len(gaps)

def compute(M,p,q):
    initial,A,B,C,final,_=M
    return int((initial*(A**p)*B*(A**q)*C*final)[0])
def family(p,q):return (1,)*(p+q+1)+(0,),(1,)*p+(0,)+(1,)*q+(0,)
def formula(p,q):return 2**q*(2*3**p+6*2**p-3)-1-int(q==0)
def main():
    models={'h':[model(((0,1),)),model(((1,0),))],
            'S2':[model(((0,1),()),((0,0),)),model(((0,),(1,)),((0,0),))],
            'E1':[model(((1,0),),((0,2),))]}
    bases=(3,2,1);grid=[(p,q) for p in (1,2,3) for q in (1,2,3)]
    design=s.Matrix([[a**p*b**q for a in bases for b in bases] for p,q in grid])
    records={};checks=[]
    for name,ms in models.items():
        value=lambda p,q:sum(compute(M,p,q) for M in ms)
        coef=design.inv()*s.Matrix([value(p,q) for p,q in grid])
        boundary=s.Matrix([[a**p for a in bases] for p in (1,2,3)]).inv()*s.Matrix([value(p,0) for p in (1,2,3)])
        records[name]={'bases':bases,'coefficients_q_positive':data(s.Matrix(3,3,list(coef))),'coefficients_q_zero':[str(c) for c in boundary],
                       'gap_counts':[M[-1] for M in ms]}
        for p in range(1,6):
            for q in range(5):
                predicted=sum(coef[i*3+j]*a**p*b**q for i,a in enumerate(bases) for j,b in enumerate(bases)) if q else sum(c*a**p for c,a in zip(boundary,bases))
                need(value(p,q)==predicted,'annihilator formula initial and extra controls')
                if name=='h':need(predicted==formula(p,q),'closed ready Hamilton formula')
                if p+q<=5:
                    w,u=family(p,q);T=old.source(w,u);v,z=values(T,len(w));actual={'h':v[0],'S2':v[3],'E1':v[2]}[name]
                    need(predicted==actual and z[3]-z[0]==1,'independent endpoint DP and readiness')
                    checks.append([name,p,q,int(predicted)])
    # Actual two-polarity three-port updates commute up to fresh-label swap.
    updates=[]
    def birth(T,ports,polarity):
        a,b,ell=ports;full=(1<<len(T))-1
        incoming=(full^(1<<ell)) if polarity else ((1<<a)|(1<<b))
        return junction.ear(T,incoming)
    for p,q in ((1,0),(2,1),(3,2)):
        w,u=family(p,q);T=old.source(w,u);m=len(w);ports=(m,m+1,m-1);n=len(T)
        X=birth(birth(T,ports,0),ports,1);Y=birth(birth(T,ports,1),ports,0)
        order=list(range(n))+[n+1,n]
        need(X==junction.sub(Y,tuple(order)),'p/q updates commute after exchanging fresh names')
        need(old.hp(X)==formula(p+1,q+1),'literal commuting updates realize formula')
        updates.append({'p':p,'q':q,'ports':ports,'source_rows':T,'output_H':old.hp(X)})
    lanes=[]
    for p,q in ((1,3),(2,1),(2,3),(3,2),(3,3)):
        h=formula(p,q);lanes.append({'p':p,'q':q,'H':h,'H_factorization':{str(a):b for a,b in s.factorint(h).items()},'NR':3*h+2,'NR_is_prime':bool(s.isprime(3*h+2))})
    out={'status':'EXACT_OPERATOR_CERTIFICATE','formula':'2^q*(2*3^p+6*2^p-3)-1-[q=0]',
         'domain':'p>=1,q>=0','full_golden_formulas':records,'independent_small_checks':checks,'literal_commuting_updates':updates,'lanes':lanes}
    (HERE/'one_hole.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf8',newline='\n')
    print('PASS: all-order annihilator certificates, independent gap/endpoint counts, literal commuting updates; formulas',records)
if __name__=='__main__':main()
