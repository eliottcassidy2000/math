"""Exact controls for the unrestricted real trace-square-zero cubic theorem."""
from pathlib import Path
from hashlib import sha256
from itertools import combinations,product
from collections import Counter
from fractions import Fraction as Q
import json,sys,random

sys.stdout.reconfigure(newline='\n')
gates=Counter()
def check(ok,label):
    gates[label]+=1
    if not ok:raise RuntimeError('always-active gate failed: '+label)
def mul(A,B):
    return [[sum(a*b for a,b in zip(row,col)) for col in zip(*B)] for row in A]
def transpose(A):return [list(c) for c in zip(*A)]
def tr(A):return sum(A[i][i] for i in range(len(A)))
def energy(A):return sum(a*a for row in A for a in row)
def triple_sums(A):
    C=T=0
    for i,j,k in combinations(range(len(A)),3):
        cyc=A[i][j]*A[j][k]*A[k][i]+A[i][k]*A[k][j]*A[j][i]
        total=(A[i][j]+A[j][i])*(A[i][k]+A[k][i])*(A[j][k]+A[k][j])
        C+=cyc;T+=total-cyc
    return C,T
def oriented_check(A,label,do_matrix=True):
    C,T=triple_sums(A);E=energy(A);F=3*C-T
    check(T>=0,label+' nonnegative transitive penalty')
    check(27*C*C<=E**3,label+' unrestricted cyclic cubic bound')
    check(F<=3*C,label+' production versus cyclic trace')
    check(F<0 or 3*F*F<=E**3,label+' actual positive production norm bound')
    if do_matrix:
        A2=mul(A,A)
        check(tr(A2)==0,label+' literal trace square zero')
        check(tr(mul(A2,A))==3*C,label+' literal trace cube versus cycles')
        At=transpose(A)
        N=[[u+v for u,v in zip(ar,br)] for ar,br in zip(A,At)]
        L=[[u-v for u,v in zip(ar,br)] for ar,br in zip(A,At)]
        check(tr(mul(N,mul(L,L)))==2*F,label+' independent shear matrix contraction')
    return (E,C,T,F)
def main():
    totals=Counter();digest=sha256()
    for n in range(1,5):
        edges=list(combinations(range(n),2))
        for word in product(range(5),repeat=len(edges)):
            A=[[0]*n for _ in range(n)]
            for (i,j),z in zip(edges,word):
                if z in (1,2):A[i][j]=z
                if z in (3,4):A[j][i]=z-2
            vals=oriented_check(A,'exhaustive small')
            digest.update(repr((n,word,vals)).encode())
            totals[str(n)]+=1
    check(dict(totals)=={'1':1,'2':5,'3':125,'4':15625},'complete declared finite universe')
    rng=random.Random(150909)
    for n in range(5,13):
        for _ in range(100):
            A=[[0]*n for _ in range(n)]
            for i,j in combinations(range(n),2):
                w=rng.randrange(6)
                if rng.randrange(2):A[i][j]=w
                else:A[j][i]=w
            vals=oriented_check(A,'larger fixed-seed scout')
            digest.update(repr((n,A,vals)).encode())
    equality=[]
    # The common-edge multi-return equality family includes zero-edge supports.
    for q in range(1,6):
        A=[[0]*(q*q+2) for _ in range(q*q+2)]
        A[0][1]=q
        for j in range(2,len(A)):A[1][j]=A[j][0]=1
        E,C,T,F=oriented_check(A,'multi-return equality')
        check(T==0 and 27*C*C==E**3 and 3*F*F==E**3,'nontriangle equality is preserved')
        check(mul(A,transpose(A))==mul(transpose(A),A),'equality matrix normality')
        # M^3 = q^2 times the orthogonal projection scaled by q; trace and
        # Frobenius tests above are independent of this finite rank witness.
        import sympy as s
        check(s.Matrix(A).rank()==3,'equality rank three')
        equality.append(dict(returns=q*q,E=E,C=C,F=F))
    # Similarity preserves trace moments while changing the energy gauge.
    A=[[Q(0),Q(2),Q(0)],[Q(0),Q(0),Q(3)],[Q(1,6),Q(0),Q(0)]]
    E,C,T,F=oriented_check(A,'nonnormal similarity')
    check(C==1 and 27*C*C<E**3,'nonnormal similarity is strict')
    check(mul(A,transpose(A))!=mul(transpose(A),A),'nonnormal defect is actual')
    # Rational orthogonal similarity goes outside the oriented cone but keeps
    # the stronger REAL matrix theorem and its equality exact.
    P=[[Q(0),Q(1),Q(0)],[Q(0),Q(0),Q(1)],[Q(1),Q(0),Q(0)]]
    outside=0
    for vv in [(1,2,3),(1,1,2),(2,-1,3),(3,2,-2)]:
        norm=sum(z*z for z in vv)
        O=[[Q(int(i==j))-Q(2*vv[i]*vv[j],norm) for j in range(3)] for i in range(3)]
        check(mul(O,transpose(O))==[[Q(int(i==j)) for j in range(3)] for i in range(3)],'rational orthogonal transport')
        A=mul(mul(O,P),transpose(O));A2=mul(A,A);cube=tr(mul(A2,A));E=energy(A)
        check(tr(A2)==0 and cube==3,'real outside-cone trace moments')
        check(3*cube*cube==E**3,'real outside-cone sharp equality')
        check(mul(A,transpose(A))==mul(transpose(A),A),'real outside-cone equality normality')
        outside+=int(any(A[i][i]!=0 for i in range(3)) or any(z<0 for row in A for z in row))
    check(outside==4,'all four real controls genuinely outside oriented cone')
    # Exact scalar critical-point arithmetic uses w^2=1/6 and t=2w.
    import sympy as s
    w=s.symbols('w',positive=True)
    f=(1-2*w*w)**s.Rational(3,2)-2*w**3+3*w
    check(s.simplify(s.diff(f,w)-3*s.sqrt(1-2*w*w)*(s.sqrt(1-2*w*w)-2*w))==0,'scalar derivative identity')
    check(s.simplify(f.subs(w,1/s.sqrt(6))-2*s.sqrt(s.Rational(2,3)))==0,'sharp scalar maximum value')
    check(s.simplify(f.subs(w,0))==1 and s.simplify(f.subs(w,1/s.sqrt(2)))==s.sqrt(2),'strict scalar endpoint values')
    cert=dict(scope='Unrestricted real trace-square-zero cubic theorem; finite controls support analytic proof.',
              source_sha256=sha256(Path(__file__).read_bytes()).hexdigest(),
              exhaustive_counts=dict(totals),larger_scout=800,semantic_sha256=digest.hexdigest(),
              equality_cases=equality,gates=dict(sorted(gates.items())),total_gates=sum(gates.values()))
    dest=Path(__file__).with_name(Path(__file__).stem+'_certificate.json')
    if Path(__file__).parent.name=='04-computation':dest=Path(__file__).parents[1]/'05-knowledge/results'/dest.name
    dest.write_text(json.dumps(cert,indent=2,sort_keys=True)+'\n',encoding='utf-8',newline='\n')
    print('EXHAUSTIVE_COUNTS',json.dumps(dict(totals),sort_keys=True),'LARGER_SCOUT',800)
    print('SEMANTIC_SHA256',digest.hexdigest())
    print('EQUALITY_CASES',json.dumps(equality,sort_keys=True))
    print('TOTAL_ALWAYS_ACTIVE_GATES',sum(gates.values()))
    print('CERTIFICATE_SHA256',sha256(dest.read_bytes()).hexdigest())
if __name__=='__main__':main()
