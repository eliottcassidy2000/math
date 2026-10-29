"""Second-phase certificates. Imports frozen first-phase exact arithmetic."""
from fractions import Fraction as F
from check import add, mul, fruit, primitive, require

G=(F(-4),F(28));T=(F(56),F(728));H=(F(-196),F(196));U=(F(0),F(0))

def add2(P,Q):return add(P,Q,421,43904)
def mul2(n,P):
    if n<0:return mul2(-n,(P[0],-P[1]))
    R=None
    while n:
        if n&1:R=add2(R,P)
        P=add2(P,P);n//=2
    return R

def phi(P):
    if P is None or P[0]==4:return None
    x,y=P
    return x*(x-56)**2/(x-4)**2,(x-56)*(x*x+44*x+224)*y/(x-4)**3

def source(n,e,j=0):
    if n%3:raise ValueError('nonzero rational isogeny-image residue')
    return add(mul(n//3,G),mul(e+2*j,T))

def check():
    counts=[]
    for p in [11,17]:
        require((16*43904**2*(421**2-4*43904))%p!=0,'good reduction')
        counts.append(1+sum(sum(y*y%p==(x**3+421*x*x+43904*x)%p for y in range(p)) for x in range(p)))
    require(counts==[12,18],'torsion count bound')
    # Exact identity for the 3-division polynomial.
    from collections import defaultdict
    coeff=defaultdict(int)
    for i,a in enumerate([784,3]):
        for j,b in enumerate([-2458624,9408,300,1]):coeff[i+j]+=a*b
    require([coeff[i] for i in range(5)]==[-43904**2,0,6*43904,4*421,3],'division polynomial')
    require(all((x**3+300*x*x+9408*x-2458624)%19 for x in range(19)),'irreducible cubic mod19')
    x=F(-784,3)
    require(x**3+421*x*x+43904*x==F(-15366400,27)<0,'linear root no rational y')
    require(mul2(3,H)==phi(G)==(F(-225),F(210)),'primitive-generator bridge')
    require(phi(T)==U,'torsion surjection')
    ncert=0
    for m in range(1,21):
        for e in range(2):
            Q=add2(mul2(3*m,H),mul2(e,U))
            triples=[]
            for j in range(3):
                P=source(3*m,e,j)
                require(phi(P)==Q,'named source lift')
                for sign in [1,-1]:
                    R=P if sign==1 else (P[0],-P[1])
                    t=primitive(fruit(R));triples.append(t)
                    a,b,c=t
                    require((a+b)*(a+c)*(b+c)!=0,'exceptional guard')
                    require(F(a,b+c)+F(b,a+c)+F(c,a+b)==4,'original expression')
            require(len(set(triples))==6,'six distinct ordered lifts')
            require(len({tuple(sorted(t)) for t in triples})==1,'single unordered triple')
            q=-Q[0]/196-1
            QU=add2(Q,U);q2=-QU[0]/196-1
            require(q2==(1-7*q)/(7*(q+1)),'reciprocity on quotient')
            if q>F(1,6) and 1+29*q-196*q*q>0:require(q2<0,'disjoint positive partners')
            ncert+=1
    for n in [1,2,4,5]:
        try:source(n,0)
        except ValueError:pass
        else:raise ArithmeticError('nonlift accepted')
    require((16*43904**2*(421**2-4*43904))%3==2,'good reduction at3')
    points={(x,y) for x in range(3) for y in range(3) if (y*y-x**3-421*x*x-43904*x)%3==0}
    require(points=={(0,0),(1,1),(1,2),(2,1),(2,2)},'finite field six')
    def reduce3(P):
        if P is None or P[0].denominator%3==0:return None
        return tuple(int(t.numerator*pow(t.denominator,-1,3)%3) for t in P)
    require(len({reduce3(mul2(n,H)) for n in range(6)})==6,'Hbar generates')
    nguards=0
    for n in range(-40,41):
        for e in range(2):
            Q=add2(mul2(n,H),mul2(e,U))
            if Q is None:continue
            q=-Q[0]/196-1
            guard=q.denominator%3==0 or (q.numerator+q.denominator)%3==0
            require(guard==(n%3==0),'constant-size ternary classifier')
            require(guard==(reduce3(Q) in [None,(0,0)]),'projective finite reduction')
            nguards+=1
    print('PASS:',nguards,'arbitrary-index mod3 guards, including denominator poles')
    print('PASS: Eprime torsion certificate C2; 3H=phi(G);',ncert,'unordered families with six exact labels; reciprocity guard')
if __name__=='__main__':check()
