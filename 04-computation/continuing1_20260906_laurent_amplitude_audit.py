"""Independent raw-source and exact algebra audit of the positive-root repair."""
from fractions import Fraction as F
from math import factorial, comb
import sys
sys.stdout.reconfigure(newline='\n')
GATES=0

def need(ok,why):
    global GATES
    GATES+=1
    if not ok: raise ArithmeticError(why)

def add(a,b):
    c=[F(0)]*max(len(a),len(b))
    for i,x in enumerate(a):c[i]+=x
    for i,x in enumerate(b):c[i]+=x
    return c

def mul(a,b):
    c=[F(0)]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):c[i+j]+=x*y
    return c

def reduce(a):
    a=list(map(F,a))+[F(0)]*max(0,3-len(a))
    for i in range(len(a)-1,2,-1):
        c=a[i];a[i]=0
        a[i-1]+=28*c;a[i-2]-=14*c;a[i-3]+=c/3
    return a[:3]

def evaluate(a,x):
    v=F(0)
    for c in reversed(a):v=v*x+c
    return v

def imul(a,b):
    v=[x*y for x in a for y in b]
    return min(v),max(v)

def iadd(a,b):return a[0]+b[0],a[1]+b[1]
def ineg(a):return -a[1],-a[0]
def ihorner(a,I):
    v=(F(0),F(0))
    for c in reversed(a):v=iadd(imul(v,I),(c,c))
    return v

P=[-F(1,3),F(14),F(-28),F(1)]
ZPOLY=[-F(1,3),F(0),F(14),F(0),F(-28),F(0),F(1)]
ABC=[F(142115562175391338833022911,115962939903341750549938130),
     F(137850584919079100136401223,115962939903341750549938130),
     -F(825111792094668079242879,23192587980668350109987626)]

def main():
    # Reconstruct first and second CT fibers by literal charge, not a row formula.
    first={};doubled={}
    for n,target in ((11,first),(22,doubled)):
        for x in range(n+1):
            for y in range(n-x+1):
                z=n-x-y
                if -21*x+y+12*z:continue
                target[z]=factorial(n)//(factorial(x)*factorial(y)*factorial(z))
    need(first=={1:110,3:4620,5:9240,7:330},'raw first fiber')
    need(all(z%2==0 for z in doubled),'second gamma exponents')
    T=[F(0)]*8
    for z,c in doubled.items():T[z//2]=c*(-1 if (z//2-1)%2 else 1)
    need(T[0]==-22,'lower carry retained')
    need(reduce([-F(110,330),F(4620,330),-F(9240,330),F(1)])==[0,0,0],'exact first normalization')

    # Different coefficient route: differentiate complete H coefficients first,
    # then sum all ordered pairs of source degrees adding to14.
    B=(1,-20,28,-10,1)
    rows=[]
    for r in range(7):
        J=[F(0)]*8
        for d in range(r,15-r):
            e=14-d
            for i in range(5):
                for j in range(5):
                    if not (0<=d-2*i<=11 and 0<=e-2*j<=11):continue
                    v=B[i]*B[j]*comb(11,d-2*i)*comb(11,e-2*j)
                    v*=factorial(d)//factorial(d-r)* (factorial(e)//factorial(e-r))
                    J[7-i-j]+=v
        rows.append(reduce(J))
    need(reduce(mul(rows[0],ABC))==reduce(T),'T=R J0 exact raw quotient')
    need(ABC[0]>0 and ABC[1]>0 and ABC[2]<0,'signed ratio coefficients')
    print('RAW', 'first',first,'second',doubled,'ratio_checked',tuple(map(str,ABC)))

    # Independent root isolation directly in z, not square roots of producer boxes.
    brackets=[(158318,158320),(695519,695522),(5243203,5243205)]
    roots=[]
    for a,b in brackets:
        lo,hi=F(a,10**6),F(b,10**6)
        need(evaluate(ZPOLY,lo)*evaluate(ZPOLY,hi)<0,'sextic sign bracket')
        sign=evaluate(ZPOLY,lo)>0
        for _ in range(48):
            mid=(lo+hi)/2
            if (evaluate(ZPOLY,mid)>0)==sign:lo=mid
            else:hi=mid
        roots.append((lo,hi))
    need(all(roots[i][1]<roots[i+1][0] for i in range(2)),'three positive roots exhaust degree')
    SI=[imul(I,I) for I in roots]
    for r,row in enumerate(rows):
        for I in SI:need(ihorner(row,I)[1]<0,'all seven windows negative')
    for I in SI:need(ihorner(ABC,I)[0]>0,'original signed ratio positive at roots')
    y=(-32559418467680575845,-813366686033280291,-19065454274144095)
    numerator=[y[2]-28*y[1]+14*y[0],y[1]-28*y[0],y[0]]
    need(numerator==[-433076656792870357777,910850350409022843369,-32559418467680575845],
         'independent Lagrange numerator')
    bounds={0:(-19914989775031856088736,-6371590793283941536065),
            6:(-88591556071110678488277235694,-31956183471936284505432288886)}
    for r in (0,6):
        total=(F(0),F(0))
        for I,ZI in zip(SI,roots):
            derivative=ihorner((14,-56,3),I)
            need(not derivative[0]<=0<=derivative[1],'simple roots for Lagrange weights')
            reciprocal=(1/derivative[1],1/derivative[0])
            value=imul(imul(imul(ihorner(rows[r],I),ihorner(numerator,I)),reciprocal),ZI)
            total=iadd(total,value)
        need(bounds[r][0]<=total[0]<=total[1]<=bounds[r][1]<0,'published half-power enclosure')
        print('HALF_POWER',r,'integer_bounds',total[0].numerator//total[0].denominator,
              -((-total[1].numerator)//total[1].denominator))

    e1=(F(0),F(0));e2=(F(0),F(0));e3=(F(1),F(1))
    for I in roots:e1=iadd(e1,I);e3=imul(e3,I)
    for i,j in ((0,1),(0,2),(1,2)):e2=iadd(e2,imul(roots[i],roots[j]))
    a,b,c=ABC
    A=iadd((a,a),imul((c,c),imul(e1,e3)))
    BB=imul((c,c),iadd(e3,ineg(imul(e1,e2))))
    C=iadd((b,b),imul((c,c),iadd(imul(e1,e1),ineg(e2))))
    for name,I in zip(('A','B','C'),(A,BB,C)):
        need(I[0]>0,'positive real algebraic coefficient')
        scale=10**6
        lo=I[0]*scale;hi=I[1]*scale
        print('AMPLITUDE',name,'millionth_bounds',lo.numerator//lo.denominator,-((-hi.numerator)//hi.denominator))
    need(imul(e3,e3)[0]<=F(1,3)<=imul(e3,e3)[1],'E3 positive branch of1/sqrt3')

    # Universal reduction checked as a polynomial identity in formal E variables
    # through exact symbolic sparse multiplication (independent 3-variable ring).
    def smul(a,b):
        out={}
        for p,x in a.items():
            for q,y in b.items():
                k=tuple(i+j for i,j in zip(p,q));out[k]=out.get(k,F(0))+x*y
        return {k:v for k,v in out.items() if v}
    def sadd(*args):
        out={}
        for a in args:
            for k,v in a.items():out[k]=out.get(k,F(0))+v
        return {k:v for k,v in out.items() if v}
    zero=(0,0,0);one={zero:F(1)}
    E=[{tuple(int(i==j) for i in range(3)):F(1)} for j in range(3)]
    coeff=[smul(E[0],E[2]),sadd(E[2],{k:-v for k,v in smul(E[0],E[1]).items()}),sadd(smul(E[0],E[0]),{k:-v for k,v in E[1].items()})]
    # Product (z+E1)(z^3-E1 z^2+E2 z-E3).
    qplus=[{k:-v for k,v in E[2].items()},E[1],{k:-v for k,v in E[0].items()},one]
    product_coeff=[{} for _ in range(5)]
    for i,x in enumerate((E[0],one)):
        for j,y in enumerate(qplus):product_coeff[i+j]=sadd(product_coeff[i+j],smul(x,y))
    for j in range(3):need(sadd(coeff[j],product_coeff[j])=={},'z4 formal reduction low coefficient')
    need(product_coeff[3]=={} and product_coeff[4]==one,'z4 formal high coefficients')

    # Irreducibility and parity obstruction arithmetic.
    need({(x**3)%7 for x in range(7)}=={0,1,6},'cubic residues mod7')
    need(all((3*x**3-1)%7 for x in range(7)),'3p irreducible mod7')
    need(-P[0]==F(1,3),'cubic norm of s')
    # v3(1/3)=-1 is odd: a rational square has even valuation at every prime.
    need(F(1,3).numerator%3!=0 and F(1,3).denominator==3,'nonsquare rational norm')
    print('FIELD', 'p irreducible mod7; norm1/3 nonsquare; p(z^2) irreducible degree6')
    print('BOUNDARY', 'positive real-algebraic branch only; nonnegative real full-six and rational positive-branch identities excluded')
    print('PASS',GATES,'independent always-active gates')

if __name__=='__main__':main()
