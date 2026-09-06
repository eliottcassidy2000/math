"""Independent determinant/cofactor and physical-entry audit; no producer imports."""
from fractions import Fraction as F
from functools import reduce
from itertools import combinations, product
from math import gcd
import sys

sys.stdout.reconfigure(newline="\n")
Q=91**6
GATES=0


def need(ok,message):
    global GATES
    GATES+=1
    if not ok:
        raise ArithmeticError(message)


def connected(n,edges):
    reached={0}
    while True:
        more={b for a,b in edges if a in reached}|{a for a,b in edges if b in reached}
        new=reached|more
        if new==reached:
            return len(new)==n
        reached=new


def determinant(matrix):
    """Integer Bareiss elimination, independent of tree product formulas."""
    A=[list(row) for row in matrix]
    n=len(A)
    sign=1
    previous=1
    for k in range(n-1):
        if A[k][k]==0:
            pivot=next((i for i in range(k+1,n) if A[i][k]),None)
            if pivot is None:
                return 0
            A[k],A[pivot]=A[pivot],A[k]
            sign=-sign
        pivot=A[k][k]
        for i in range(k+1,n):
            for j in range(k+1,n):
                value=A[i][j]*pivot-A[i][k]*A[k][j]
                if value%previous:
                    raise ArithmeticError('Bareiss exact divisibility')
                A[i][j]=value//previous
        for i in range(k+1,n):
            A[i][k]=0
        previous=pivot
    return sign*A[-1][-1]


def minors(n,edges,weights):
    rows=[]
    for (u,v),(a,b) in zip(edges,weights):
        row=[0]*n
        row[u],row[v]=a,-b
        rows.append(row)
    return [abs(determinant([[a for j,a in enumerate(row) if j!=v] for row in rows])) for v in range(n)]


def check_tree(n,C,edges,weights):
    cof=minors(n,edges,weights)
    ell=sum(sum(v in e for e in edges)==1 for v in range(n))
    need(min(cof)>0,"positive determinant cofactors")
    for (u,v),(a,b) in zip(edges,weights):
        need(a*cof[u]==b*cof[v],"determinants form literal positive kernel")
    m=min(cof)
    need(m*(n-1)**(n-1)<=C**(n-1)*(n-2)**(n-2),"universal minimum cofactor bound")
    need(m**ell*ell**(ell*(n-1))<=C**(ell*(n-1))*(ell-1)**((ell-1)*(n-1)),"topology-dependent leaf bound")
    primitive=[v//reduce(gcd,cof) for v in cof]
    need(min(primitive)<=m,"primitive normalization direction")


def tree_controls():
    count=0
    census={}
    for n,C in ((3,6),(4,5),(5,3)):
        candidates=list(combinations(range(n),2))
        weights=[(a,b) for a in range(1,C) for b in range(1,C-a+1)]
        trees=[edges for edges in combinations(candidates,n-1) if connected(n,edges)]
        need(len(trees)==n**(n-2),"independent labelled tree enumeration")
        for edges in trees:
            for ws in product(weights,repeat=n-1):
                check_tree(n,C,edges,ws)
                count+=1
        census[n]=len(trees)
    for n in range(3,11):
        edges=tuple((0,j) for j in range(1,n))
        cof=minors(n,edges,((1,n-2),)*(n-1))
        need(min(cof)==(n-2)**(n-2) and max(cof)==(n-2)**(n-1),"sharp cofactor stars from independent minors")
        need(min(v//reduce(gcd,cof) for v in cof)==1,"sharpness does not transfer to primitive minimum")
    cof=minors(5,((0,1),(0,2),(0,3),(0,4)),((1,3),)*4)
    need(min(cof)==27>2**4,"naive half-sum product hostile")
    cof=minors(3,((0,1),(0,2)),((1,100),)*2)
    need(max(cof)>F(101**2,4)>=min(cof),"minimum cannot be replaced by maximum")
    return count,census


def allowed_sums():
    prime=[True]*357
    prime[0]=prime[1]=False
    for p in range(2,19):
        if prime[p]:
            for k in range(p*p,357,p):
                prime[k]=False
    vals={1}
    for p in range(2,357):
        if prime[p] and p%3==2:
            vals={v*p**e for v in vals for e in range(3) if v*p**e<=356}
    return vals


def label_components(row,sums):
    edges=[(i,j) for i,j in combinations(range(len(row)),2) if (row[i]+row[j])//gcd(row[i],row[j]) in sums]
    unseen=set(range(len(row)))
    out=[]
    while unseen:
        reached={min(unseen)}
        while True:
            new=reached|{j for i,j in edges if i in reached}|{i for i,j in edges if j in reached}
            if new==reached:
                break
            reached=new
        unseen-=reached
        out.append(frozenset(reached))
    return out


def box_member(a,b,x,limit=Q):
    # Complement-coordinate implementation: x=Q(a+b)-(au+bv), 0<=u,v<=2Q.
    M=limit*(a+b)-x
    u=(pow(a,-1,b)*M)%b
    v=(M-a*u)//b
    lo=max(-(u//b),-((-(v-2*limit))//a))
    hi=min((2*limit-u)//b,v//a)
    return lo<=hi


def no_crossing(A,B,Y):
    d=gcd(A,B)
    a,b=sorted((A//d,B//d))
    need(b<=Q,"physical internal height")
    delta=gcd(d,Y)
    c,x=d//delta,Y//delta
    return c>Q or not box_member(a,b,x)


def norm(x):
    return min(x%1,(-x)%1)


def actual_controls():
    larger=(6,10,15,20,24,30,60,18,2,3,40,45)
    shapes=((1,),(2,3),(2,3,6),(2,3,4,6),(2,3,4,6,9),(2,3,4,6,9,12))
    eta=(F(1,2),F(1,5),F(1,4),F(1,5),F(1,5),F(1,5))
    sums=allowed_sums()
    reports=[]
    for a,V in enumerate(shapes,1):
        b=13-a
        U=tuple(sorted(larger[:b]))
        t=120*Q+1
        row=tuple(t*v for v in V)+U
        need(reduce(gcd,V)==reduce(gcd,U)==reduce(gcd,row)==1 and 1 not in U,"actual primitive nonunit larger shape")
        need(len(set(row))==13 and sum(row)<=Q**2,"full physical box")
        need(set(label_components(row,sums))=={frozenset(range(a)),frozenset(range(a,13))},"independent actual graph")
        need(t>2*Q*max(U),"nonzero smaller partial sum dominates every larger partial sum")
        count=0
        for inside,outside in ((row[:a],row[a:]),(row[a:],row[:a])):
            for A,B in combinations(inside,2):
                for Y in outside:
                    need(no_crossing(A,B,Y),"independent complement-box mixed support test")
                    count+=1
        need(count==11*a*b//2,"complete orientation universe")
        L=max(U)
        u=min((u for u in U if u<L),key=lambda u:gcd(u,L))
        D=gcd(u,L)
        need(all(gcd(z,L)>1 for z in U if z<L),"all maximum pairs lack a coprime partner")
        delta=gcd(D,t*min(V))
        A,B=u//D,L//D
        R=Q*(A+B)-(A-1)*(B-1)
        need(D//delta<=Q and a*delta*R>=7*(b+1)*L*min(V),"native theorem hypotheses at actual row")
        # Independent grid selector: nearest integer to t*zeta-eta, so its
        # image lies within half a grid spacing of the chosen safe phase.
        zeta=F(1,7)
        center=t*zeta-eta[a-1]
        j=(center+F(1,2)).numerator//(center+F(1,2)).denominator
        x=(eta[a-1]+j)/t
        clearance=min(norm(v*x) for v in row)
        need(clearance>F(1,14),"independent literal full-row strict safe phase")
        reports.append((a,b,D,count,str(x),str(clearance)))
    return reports


def main():
    count,census=tree_controls()
    bounds=(1,177,31684,6684150,1694040507,468424857663)
    caps=(1,2,4,9,30,90)
    expected=(6240321451,76388115,698294,4854,26,0)
    for a,(M,C,cut) in enumerate(zip(bounds,caps,expected),1):
        b=13-a
        need(min(Q//C,a*Q//(7*(b+1)*M))==cut,"exact uniform endpoint-gcd cutoffs")
        if cut:
            need(C*cut<=Q and 7*(b+1)*cut*M<=a*Q,"actual inequality direction")
    # Exact divisibility refinement, independent finite positive-integer bank.
    div_cases=0
    for R in range(1,31):
        for delta in range(1,13):
            for v in range(1,13):
                z=gcd(delta,v)
                d,e=delta//z,v//z
                low=d*(R//e+1)
                candidates=[t for t in range(1,low+1) if t*v%delta==0 and t*v>delta*R]
                need(candidates==[low],"first scale satisfying retained divisibility and strict radius")
                div_cases+=1
    box_cases=0
    for limit in range(2,8):
        for a in range(1,limit):
            for b in range(a+1,limit+1):
                if gcd(a,b)!=1:
                    continue
                literal={a*r+b*s for r in range(-limit,limit+1) for s in range(-limit,limit+1)}
                for x in range(-limit*(a+b)-1,limit*(a+b)+2):
                    need(box_member(a,b,x,limit)==(x in literal),"independent complement box versus literal signed coefficients")
                    box_cases+=1
    reports=actual_controls()
    print('AUDIT PASS: leaf cofactor inequality and unit-free maximum-endpoint entry')
    print('INDEPENDENT TREE UNIVERSE:',count,'weighted labelled trees; topology census',census)
    print('MECHANISM: exact Bareiss maximal minors, independent subset tree enumeration')
    print('DIVISIBILITY_REFINEMENT_CONTROLS:',div_cases,'SIGNED_BOX_CONTROLS:',box_cases)
    for report in reports:
        print('ACTUAL_CONTROL:',report)
    print('CUTOFFS:',expected)
    print('SCOPE: sharp real cofactors; sufficient actual entry; no automatic balanced closure')
    print('ACTIVE GATES:',GATES)


if __name__=='__main__':
    main()
