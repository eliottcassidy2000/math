"""Exact four-atom positive-branch repair for actual support (-27,1,15).

All numerical optimization belongs to the separate exploratory scout.
This producer uses literal charge counts, complete carrier coefficients,
rational root intervals, integer square roots, and exact determinant bounds.
"""
from fractions import Fraction as F
from itertools import permutations
from math import comb,factorial,isqrt
from pathlib import Path
import argparse
import json
import sys

sys.stdout.reconfigure(newline='\n')
GATES=0


def need(ok,label):
    global GATES
    GATES+=1
    if not ok:raise ArithmeticError(label)


def trim(row):
    row=list(map(F,row))
    while len(row)>1 and row[-1]==0:row.pop()
    return row


def add(*rows):
    return trim([sum(row[i] if i<len(row) else 0 for row in rows)
                 for i in range(max(map(len,rows)))])


def scale(row,c):return trim([c*v for v in row])


def cv(a,b):
    out=[F(0)]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):out[i+j]+=x*y
    return trim(out)


def rem(a,p):
    a=trim(a)
    while len(a)>=len(p) and any(a):
        c=a[-1]/p[-1];shift=len(a)-len(p)
        for i,v in enumerate(p):a[i+shift]-=c*v
        a=trim(a)
    return a


def choose(n,k):return comb(n,k) if 0<=k<=n else 0
def ff(n,r):return factorial(n)//factorial(n-r) if 0<=r<=n else 0


def ia(*rows):return sum(row[0] for row in rows),sum(row[1] for row in rows)
def minus(row):return -row[1],-row[0]


def im(A,B):
    vals=[a*b for a in A for b in B]
    return min(vals),max(vals)


def divide(A,B):
    if B[0]<=0<=B[1]:return None
    return im(A,(1/B[1],1/B[0]))


def horner(row,I):
    out=(F(0),F(0))
    for c in reversed(row):out=ia(im(out,I),(c,c))
    return out


def at(row,x):return horner(row,(x,x))[0]


def root_interval(I):
    M=10**25
    a=isqrt(I[0].numerator*M*M//I[0].denominator)
    b=isqrt(I[1].numerator*M*M//I[1].denominator)+1
    lo,hi=F(a,M),F(b,M)
    need(0<lo and lo*lo<=I[0]<=I[1]<=hi*hi,'literal rational square-root enclosure')
    return lo,hi


PERMS=[(p,(-1)**sum(p[i]>p[j] for i in range(4) for j in range(i+1,4)))
       for p in permutations(range(4))]


def determinant(matrix):
    terms=[]
    for p,sign in PERMS:
        term=(F(1),F(1))
        for i in range(4):term=im(term,matrix[i][p[i]])
        terms.append(term if sign==1 else minus(term))
    return ia(*terms)


def actual_rows():
    h,x,q,k,g,n=4,1,5,9,14,24
    beta={j:comb(13-2*j,1+j) for j in range(-1,5)}
    kernel=[(-1)**(q-i)*beta[h-i] for i in range(q+1)]
    need(kernel==[-1,35,-84,55,-13,1],'full genuine beta coefficients including both linear/quadratic anchors')
    H=[trim([sum(choose(g,j-2*i)*kernel[i] if q-i==e else 0 for i in range(q+1))
             for e in range(q+1)]) for j in range(n+1)]
    first=[F((-1)**j*factorial(g),factorial(1+j)*factorial(12-3*j)*factorial(1+2*j)) for j in range(5)]
    P=scale(first,1/first[-1])
    need(P==[F(1,11),-10,84,-60,1],'actual monic quartic')
    need(all((charge-1)%14==0 for charge in (-27,1,15)) and -27+12+15==0,
         'first support mass is exactly fourteen')
    # The mod-three monic reduction is s^4+2s+2. A quartic without
    # linear or irreducible quadratic factors is irreducible.
    modrow=[2,2,0,0,1]
    need(all(sum(c*a**i for i,c in enumerate(modrow))%3 for a in range(3)),
         'mod-three quartic has no linear factor')
    irreducible_quadratics=[]
    for c in range(3):
        for b in range(3):
            qpoly=[c,b,1]
            if any((a*a+b*a+c)%3==0 for a in range(3)):continue
            irreducible_quadratics.append(qpoly)
            residue=list(modrow)
            for j in range(4,1,-1):
                lead=residue[j]%3
                for pos,v in enumerate(qpoly):residue[j-2+pos]=(residue[j-2+pos]-lead*v)%3
            need(any(residue[:2]),'mod-three quartic has no irreducible quadratic factor')
    need(len(irreducible_quadratics)==3,'complete mod-three irreducible quadratic universe')
    need(isqrt(11)**2!=11,'root norm one-eleventh is not a rational square')
    T=[F((-1)**(j%2)*factorial(28),factorial(2+j)*factorial(24-3*j)*factorial(2+2*j)) for j in range(-1,9)]
    need(H[k]==scale([0]+first,-1),'original first-row condition without changing the Euler scale')
    windows=[]
    for r in range(k):
        raw=add(*[scale(cv(H[j],H[2*k-j]),ff(j,r)*ff(2*k-j,r)) for j in range(r,2*k-r+1)])
        need(raw[0]==0,'actual square carries one common positive phase power')
        windows.append(trim(raw[1:]))
    W=[F((-1)**(j%2)*comb(28,2+2*j)*sum(beta[a]*beta.get(j-a,0) for a in beta)) for j in range(-1,9)]
    need(W==windows[0],'independent complete beta-path convolution reproduces alpha-completed square')
    for mass,z in ((14,1),(28,2)):
        literal={}
        for a in range(mass+1):
            for b in range(mass-a+1):
                c=mass-a-b
                if -27*a+b+15*c:continue
                need((c-z)%2==0,'literal source-channel parity')
                literal[(c-z)//2]=factorial(mass)//(factorial(a)*factorial(b)*factorial(c))
        expected={j:int((-1)**j*first[j]) for j in range(5)} if mass==14 else {j:int((-1)**(j%2)*T[j+1]) for j in range(-1,9)}
        need(literal==expected,'literal constant terms retain every first/doubled coefficient and the lower carry')
    return P,T,windows


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--certificate')
    args=parser.parse_args()
    P,T,windows=actual_rows()
    R=[rem(row,P) for row in windows]
    target=rem(T,P)
    labels=[(0,0),(1,7),(1,8),(3,2)]
    normalizers=[F(ff(9,r)**2) for d,r in labels]
    Is=[(F(8419,849544),F(11993,1210189)),(F(259526,2155711),F(291249,2419213)),
        (F(11376199,8744207),F(1124341,864214)),(F(120947883,2065060),F(91906768,1569213))]
    need(all(Is[i][1]<Is[i+1][0] for i in range(3)),'four disjoint positive root brackets')
    for I in Is:need(at(P,I[0])*at(P,I[1])<0,'opposite endpoint signs exhaust degree four')
    for refinement in range(100):
        Z=[root_interval(I) for I in Is]
        matrix=[];rhs=[];valid=True
        for I,zI in zip(Is,Z):
            base=horner(R[0],I)
            if base[1]>=0:valid=False;break
            row=[]
            for (d,r),norm in zip(labels,normalizers):
                if r==d==0:
                    row.append((F(1),F(1)));continue
                ratio=divide(horner(R[r],I),im((norm,norm),base))
                if ratio is None:valid=False;break
                power=(F(1),F(1))
                for _ in range(d):power=im(power,zI)
                row.append(im(power,ratio))
            if not valid:break
            matrix.append(row);rhs.append(divide(horner(target,I),base))
        if valid:
            det=determinant(matrix)
            numerators=[]
            for j in range(4):
                replacement=[[rhs[i] if a==j else matrix[i][a] for a in range(4)] for i in range(4)]
                numerators.append(determinant(replacement))
            if det[1]<0 and all(v[1]<0 for v in numerators):
                coefficients=[divide(v,det) for v in numerators]
                if all(c[0]>0 and c[1]-c[0]<F(1,10**7) for c in coefficients):break
        refined=[]
        for a,b in Is:
            m=(a+b)/2
            need(at(P,m)!=0,'nonrational dyadic midpoint')
            refined.append((a,m) if at(P,a)*at(P,m)<0 else (m,b))
        Is=refined
    else:raise ArithmeticError('determinant enclosure budget exhausted')
    need(det[1]<0,'exact interpolation determinant is nonzero')
    for D,c in zip(numerators,coefficients):
        need(D[1]<0 and c[0]>0,'strict positive algebraic Cramer coefficient')
    for I in Is:
        need(at(P,I[0])*at(P,I[1])<0,'retained root bracket')
        for r in range(9):need(horner(R[r],I)[1]<0,'every original window negative at every first root')
    print('ACTUAL support=(-27,1,15) firstmass=14 doubledmass=28; all source coefficients and carries retained')
    print('P=s^4-60s^3+84s^2-10s+1/11; H degree24 selected coefficient9')
    print('ATOMS (z_power,derivative_order)=',labels,'normalizers=',','.join(map(str,normalizers)))
    print('REFINEMENTS',refinement)
    for i,I in enumerate(Is):print('ROOT',i+1,','.join(map(str,I)))
    print('DETERMINANT_SIGN negative; all four Cramer numerator signs negative')
    for j,c in enumerate(coefficients):
        promised=[(984869,984870),(421733,421734),(4973466,4973467),(5162,5163)][j]
        need(F(promised[0],10**6)<=c[0]<=c[1]<=F(promised[1],10**6),'published rational coefficient enclosure')
        lo,hi=c[0]*10**6,c[1]*10**6
        print('COEFFICIENT',j,'millionth_enclosure=',lo.numerator//lo.denominator,-((-hi.numerator)//hi.denominator))
    print('EXACT_FORM T=beta0 J0+beta1 z J7/(9)_7^2+beta2 z J8/(9)_8^2+beta3 z^3 J2/(9)_2^2 on the four positive branches')
    print('FIELD p(s) irreducible mod3; Norm(s)=1/11 nonsquare, hence p(z^2) irreducible over Q')
    print('SCOPE positive real algebraic coefficients; not an all-eight-root identity; no universal h or rational-coefficient claim')
    print('PASS',GATES,'always-active gates; no floating point or optimization in this producer')
    if args.certificate:
        bank={'P':P,'full_target':T,'full_windows':windows,'residues':R,'target':target,
              'labels':labels,'normalizers':normalizers,'intervals':Is,'square_roots':Z,
              'matrix':matrix,'rhs':rhs,'determinant':det,'numerators':numerators,'coefficients':coefficients}
        def encode(v):
            if isinstance(v,F):return str(v)
            if isinstance(v,dict):return {k:encode(w) for k,w in v.items()}
            if isinstance(v,(list,tuple)):return [encode(w) for w in v]
            return v
        Path(args.certificate).write_bytes((json.dumps(encode(bank),sort_keys=True,indent=2)+'\n').encode())


if __name__=='__main__':main()
