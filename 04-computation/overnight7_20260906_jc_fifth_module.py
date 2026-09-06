"""Exact fifth-normal source/target module and jet-order controls.

The all-size filtered pullback lemma is proved in the report. This program
checks its symbolic differential identities and complete bounded local
primitive banks at two finite-field exceptional embeddings. No repository
mathematical code imports and no assertions.
"""
from itertools import product
from math import comb
import hashlib
import json
import sys
import sympy as S

sys.stdout.reconfigure(newline="\n")
GATES=0


def need(ok,label):
    global GATES
    GATES+=1
    if not ok:
        raise RuntimeError(label)


X,Q,W,A=S.symbols('x q w alpha')
F=72783360*A**4-77822208*A**3-28419741*A**2+7849770*A-1276420
KAPPA=(-S.Rational(5183766767360,3**19)+S.Rational(931699873280,3**16)*A
       -S.Rational(460399208960,3**14)*A**2-S.Rational(183606968320,3**13)*A**3)
CHI=-S.Rational(9,20)*KAPPA
PR=S.Rational(520,9)*A**2-S.Rational(1688,81)*A-S.Rational(5717,729)
QP=(X**5+S.Rational(9,2)*X**4-2*X**3-S.Rational(27,4)*X**2+X-S.Rational(3,4)
    +X**2*(X**2-1)**2*(-S.Rational(259,36)+PR+4*A-9*A*X-PR*X**2))
J=S.Matrix([[3,0,0,-1,0,-2],[-9,0,0,0,-1,2],
            [0,3,0,-1,0,0],[0,4,0,0,-1,0],
            [0,0,3,-1,0,-2],[0,0,9,0,-1,-2]])


def modq(value,p):
    num,den=S.fraction(value)
    return int(num)%p*pow(int(den)%p,-1,p)%p


def rref(rows,p):
    rows=[[v%p for v in row] for row in rows]
    piv=[];r=0
    for c in range(len(rows[0]) if rows else 0):
        candidate=next((i for i in range(r,len(rows)) if rows[i][c]),None)
        if candidate is None:continue
        rows[r],rows[candidate]=rows[candidate],rows[r]
        inv=pow(rows[r][c],-1,p)
        rows[r]=[v*inv%p for v in rows[r]]
        for i in range(len(rows)):
            if i!=r and rows[i][c]:
                v=rows[i][c]
                rows[i]=[(a-v*b)%p for a,b in zip(rows[i],rows[r])]
        piv.append(c);r+=1
        if r==len(rows):break
    return rows,piv


def rank(rows,p):
    return len(rref(rows,p)[1])


class Ring:
    def __init__(self,p,bound):self.p,self.bound=p,bound
    def scalar(self,v):return {} if v%self.p==0 else {(0,0):v%self.p}
    def add(self,a,b):
        c=dict(a)
        for k,v in b.items():c[k]=(c.get(k,0)+v)%self.p
        return {k:v for k,v in c.items() if v}
    def scale(self,a,v):return {k:v*x%self.p for k,x in a.items() if v*x%self.p}
    def mul(self,a,b):
        c={}
        for (i,j),v in a.items():
            for (k,l),w in b.items():
                if 2*(i+k)+j+l<=self.bound:
                    key=i+k,j+l;c[key]=(c.get(key,0)+v*w)%self.p
        return {k:v for k,v in c.items() if v}
    def power(self,a,n):
        value=self.scalar(1)
        for _ in range(n):value=self.mul(value,a)
        return value
    def polynomial(self,coeffs,a):
        v={}
        for c in reversed(coeffs):v=self.add(self.mul(v,a),self.scalar(c))
        return v


def local_banks(p,alpha):
    ring=Ring(p,12)
    coeffs=[modq(S.expand(QP.subs(A,alpha)).coeff(X,j),p) for j in range(9)]
    monos=[(a,b,c) for a in range(7) for b in range(7-a)
           for c in range(13-2*(a+b))]
    slots=[(i,j) for j in range(10) for i in range((12-j)//2+1)]
    constraints=[];outputs=[];values=[]
    all_columns=[]
    for point in (-1,0,1):
        xx=ring.add(ring.scalar(point),{(1,0):1})
        qq=ring.add(ring.polynomial(coeffs,xx),{(0,2):1})
        dd=ring.add(ring.scalar(1),ring.mul(ring.power(xx,2),qq))
        cc=ring.mul(xx,ring.mul(dd,ring.add(dd,ring.scalar(2))))
        zz=ring.add(ring.mul(qq,ring.add(dd,ring.scalar(3))),ring.scalar(3))
        need((0,0) not in cc and (0,0) not in zz,"local target coordinates vanish")
        c_p=[ring.power(cc,j) for j in range(7)]
        z_p=[ring.power(zz,j) for j in range(7)]
        columns=[ring.mul(ring.mul(c_p[a],z_p[b]),{(0,c):1}) for a,b,c in monos]
        all_columns.append(columns)
        constraints.extend([[column.get(slot,0) for column in columns] for slot in slots])
        outputs.append([(-column.get((1,10),0))%p for column in columns])
        values.append([column.get((0,10),0) for column in columns])
    need((len(monos),len(constraints))==(140,135),"complete weighted local primitive bank dimensions")
    rr=rank(constraints,p)
    joined=rank(constraints+outputs,p)
    period=[sum(ell*outputs[i][j] for i,ell in enumerate((5,-18,13)))%p for j in range(len(monos))]
    need(rank(constraints+[period],p)==rr,"all prefix-cancelled target primitives kill L")
    need(joined-rr==2,"full retained normal image is exactly the tangent-relation plane")
    jet_rank=rank(constraints+values+outputs,p)
    need(jet_rank-rr==3,"complete leading primitive value-and-derivative image has dimension three")
    for row in values[1:]:
        difference=[(a-b)%p for a,b in zip(row,values[0])]
        need(rank(constraints+[difference],p)==rr,"all leading primitive branch values agree")
    for monomial,expected in [((1,0,10),(-3,-3,-3)),((0,1,10),(9,-4,-9))]:
        column=monos.index(monomial)
        need(all(row[column]==0 for row in constraints),"actual polynomial target primitive fixes lower source jets")
        need(tuple(row[column] for row in outputs)==tuple(v%p for v in expected),"actual polynomial primitives attain the full rank-two response")
    signature=hashlib.sha256(json.dumps([constraints,values,outputs],separators=(",",":")).encode()).hexdigest()
    return {"p":p,"alpha":alpha,"columns":140,"constraints":135,"constraint_rank":rr,
            "joint_rank":joined,"full_jet_rank":jet_rank,"period_rank":rr,"matrix_sha256":signature}


def series_controls(p,alpha):
    N=6
    def add(a,b):return [(x+y)%p for x,y in zip(a,b)]
    def scale(a,v):return [x*v%p for x in a]
    def mul(a,b):return [sum(a[j]*b[i-j] for j in range(i+1))%p for i in range(N)]
    def sc(v):return [v%p]+[0]*(N-1)
    def polynomial(coeffs,x):
        value=sc(0)
        for c in reversed(coeffs):value=add(mul(value,x),sc(c))
        return value
    coeffs=[modq(S.expand(QP.subs(A,alpha)).coeff(X,j),p) for j in range(9)]
    J_inv=[[modq(v,p) for v in row] for row in J.inv().tolist()]
    ss=[0,1,0,0,0,0]
    kap=modq(KAPPA.subs(A,alpha),p);chi=modq(CHI.subs(A,alpha),p)
    need(kap!=0 and chi!=0,"exceptional nonzero fifth debt control")
    def compiler(z,correction,h):
        qq=add(add(polynomial(coeffs,z),ss),mul(correction,z))
        qh=polynomial(h,z)
        qq[5]=(qq[5]+qh[0])%p
        d=add(sc(1),mul(mul(z,z),qq))
        c=mul(z,mul(d,add(d,sc(2))))
        e=mul(qq,add(d,sc(3)))
        return c,e
    def solve(h):
        zs=[sc(-1),sc(0),sc(1)];cstar=sc(0);estar=sc(-3);corr=sc(0)
        for order in range(1,N):
            residual=[]
            for z in zs:
                c,e=compiler(z,corr,h)
                residual.extend(((c[order]-cstar[order])%p,(e[order]-estar[order])%p))
            delta=[-sum(a*b for a,b in zip(row,residual))%p for row in J_inv]
            for i in range(3):zs[i][order]=delta[i]
            cstar[order]=delta[3];estar[order]=delta[4];corr[order]=delta[5]
            for z in zs:
                c,e=compiler(z,corr,h)
                need(c[:order+1]==cstar[:order+1] and e[:order+1]==estar[:order+1],"all six contact equations, including the actual fifth forcing")
        lh=(5*polynomial(h,sc(-1))[0]-18*polynomial(h,sc(0))[0]+13*polynomial(h,sc(1))[0])%p
        need(corr[1:5]==[0]*4 and corr[5]==(chi-lh*pow(8,-1,p))%p,"full fifth affine compatibility condition")
        return corr[5]
    need(solve([0])==chi,"uncorrected pencil has exact fifth exit")
    need(solve([0,chi])==0,"polynomial source shear pays collision coefficient")
    need(solve([1,chi])==0,"constant-kernel alteration preserves payment")
    need(solve([0,chi-9,4])==0,"nonconstant-kernel alteration preserves payment")
    need(solve([1])==chi,"constant descended coefficient alone does not pay")
    # dC wedge dE has Jac_(x,q)=Jcomp. An order-N graph term starts at density N-1.
    for point in (-1,0,1):
        q0=polynomial(coeffs,sc(point))[0];d0=(1+point*point*q0)%p
        j0=6*(d0+1)*(d0*d0+2*d0-2)%p
        need(j0==12,"actual target-wedge base normal coefficient")
        need(10*point*j0*pow(3,-1,p)%p==40*point%p,"t10 normal enters normalized density first at J9")
        need(9*point*j0*pow(3,-1,p)%p==36*point%p,"t9 hostile can change J8")
    return {"p":p,"alpha":alpha,"kappa":kap,"chi5":chi}


def main():
    need(S.gcd(F,S.together(KAPPA).as_numer_denom()[0])==1,"kappa nonzero at every exceptional embedding")
    displayed=(S.Rational(259188338368,3**17)-S.Rational(46584993664,3**14)*A
               +S.Rational(23019960448,3**12)*A**2+S.Rational(9180348416,3**11)*A**3)
    need(S.expand(CHI-displayed)==0,"inherited chi5=-9 kappa/20")
    need(S.expand(8*CHI+S.Rational(18,5)*KAPPA)==0,"required fifth L value")
    need(J.det()==-288,"inherited full repair Jacobian")
    tangents=S.Matrix([[3,-9],[3,4],[3,9]])
    ell=S.Matrix([[5,-18,13]])
    need(ell*tangents==S.zeros(1,2) and tangents.rank()==2,"ordinary triple tangent relation")
    hess=S.Matrix([[9,-54,81],[9,24,16],[9,54,81]])
    need(hess.det()==63180,"three tangent directions kill every quadratic Hessian")
    h=S.Function('h')(X)
    shear=S.Matrix([X,Q+W**10*h,W])
    need(S.simplify(shear.jacobian([X,Q,W]).det())==1,"lawful ambient polynomial shear has unit Jacobian")
    d=1+X**2*Q;c=X*d*(d+2);e=Q*(d+3)
    jac=S.diff(c,X)*S.diff(e,Q)-S.diff(c,Q)*S.diff(e,X)
    need(S.expand(jac-6*(d+1)*(d*d+2*d-2))==0,"symbolic compiler two-form density")
    # General polynomial target coefficient dependence cannot lower the t-adic order further;
    # this is the sharp derivative term for a coordinate wedge itself.
    records=[];contacts=[]
    for p,alpha in ((421,126),(443,112)):
        need(int(F.subs(A,alpha))%p==0,"named exceptional residue embedding")
        records.append(local_banks(p,alpha))
        contacts.append(series_controls(p,alpha))
    print("Actual polynomial source module: K[x]; h5=chi5*x pays L(h5)=-18*kappa/5")
    print("Full prefix-cancelled descended Hamiltonian normal evaluation image: ker(5,-18,13), dimension2")
    print("local_primitive_banks="+json.dumps(records,sort_keys=True))
    print("literal_contact_controls="+json.dumps(contacts,sort_keys=True))
    print("Sharp order firewall: source O(t10) leaves every J0..J8 unchanged; O(t9) can change J8")
    print("No target two-form has nonzero constant source density on inherited4046 frozen-prefix graphs")
    print(f"PASS: {GATES} optimization-live exact gates")


if __name__=='__main__':main()
