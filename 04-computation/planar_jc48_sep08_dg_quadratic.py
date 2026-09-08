#!/usr/bin/env python3
"""Exact DG filtration controls; the inherited quadratic theorem is analytic."""
import hashlib
import json
import sympy as s

x,t,r,b=s.symbols('x t r b')
gates=0
def need(ok,label):
    global gates
    gates+=1
    if not ok:
        raise RuntimeError(label)

def chart(f):
    return s.expand(f.subs({x:1/r,t:-r**2-r**4*b},simultaneous=True))

rows=[]
for n in range(4):
    basis=[x**a*t**(n-j)*(1+x*x*t)**j for a in range(2*n+1) for j in range(n+1)]
    for a in range(2*n+1):
        for j in range(n+1):
            f=x**a*t**(n-j)*(1+x*x*t)**j
            need(s.expand(chart(f)-(-1)**n*r**(2*n-a)*b**j*(1+r*r*b)**(n-j))==0,'actual chart identity')
    # Independent coefficient kernel in the complete necessary degree box.
    # Regularity of the highest b power and descending elimination bound
    # every source coefficient by 4n; this larger box includes all candidates.
    mons=[x**a*t**k for k in range(n+1) for a in range(4*n+1)]
    columns=[]; keys=set()
    for mon in mons:
        d={}
        for term in s.Add.make_args(chart(mon)):
            powers=term.as_powers_dict(); ir=int(powers.get(r,0)); ib=int(powers.get(b,0))
            if ir<0:
                d[(ir,ib)]=d.get((ir,ib),0)+term/(r**ir*b**ib)
        columns.append(d); keys.update(d)
    matrix=s.Matrix([[d.get(key,0) for d in columns] for key in sorted(keys)]) if keys else s.zeros(0,len(mons))
    dimension=len(mons)-matrix.rank()
    need(dimension==(2*n+1)*(n+1),'independent Laurent kernel dimension')
    coefficients=s.Matrix([[s.Poly(s.expand(f),x,t).coeff_monomial(mon) for f in basis] for mon in mons])
    need(coefficients.rank()==dimension,'basis independent and exhaustive in necessary box')
    need(matrix*coefficients==s.zeros(matrix.rows,dimension),'basis in actual kernel')
    rows.append([n,len(mons),matrix.rank(),dimension])

aa=s.symbols('a0:9'); bb=s.symbols('b0:5'); c0=s.symbols('c0')
A=sum(aa[i]*x**i for i in range(9))
B=2*aa[8]*x**6+2*aa[7]*x**5+sum(bb[i]*x**i for i in range(5))
C=aa[8]*x**4+aa[7]*x**3+(bb[4]-aa[6])*x*x+(bb[3]-aa[5])*x+c0
F=A*t*t+B*t+C
Fc=chart(F)
need(all(term.as_powers_dict().get(r,0)>=0 for term in s.Add.make_args(Fc)),'full free family globally regular')
need(s.expand(Fc.subs(r,0)-(aa[8]*b*b+(2*aa[6]-bb[4])*b+aa[4]-bb[2]+c0))==0,'full boundary restriction')
need(s.expand(s.diff(t,x)*s.diff(-x,t)-s.diff(t,t)*s.diff(-x,x))==1,'positive affine plane mate')
need(chart(-x)==-1/r,'negative global mate control')

# Literal multiplication splitting of the bidegree boxes (n=1 with n=2).
for a in range(7):
    for j in range(4):
        a1=min(a,2); a2=a-a1; j1=min(j,1); j2=j-j1
        need(0<=a2<=4 and 0<=j2<=2,'box splitting admissible')
        f1=x**a1*t**(1-j1)*(1+x*x*t)**j1
        f2=x**a2*t**(2-j2)*(1+x*x*t)**j2
        need(s.expand(f1*f2-x**a*t**(3-j)*(1+x*x*t)**j)==0,'actual product filtration')

payload={'kernel_rows':rows,'L2_parameters':15,'boundary':'a8*b^2+(2*a6-b4)*b+a4-b2+c0','gates':gates}
digest=hashlib.sha256(json.dumps(payload,sort_keys=True,separators=(',',':')).encode()).hexdigest()
print('DG complete source filtration: (n, universe, rank, dimension) =',rows)
print('Quadratic layer: 15 free parameters; boundary and pole controls PASS')
print('Quadratic Keller exclusion is RECOVERED from THM-2071, not inferred from finite tests')
print('Always-active exact gates:',gates)
print('Semantic SHA256:',digest)
print('RESULT: PASS')
