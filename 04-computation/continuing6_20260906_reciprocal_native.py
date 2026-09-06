"""Reciprocal output swap 0<->1: native packets, characters and templates.

No repository imports. All gate checks remain active under python -O.
"""
from collections import defaultdict
from fractions import Fraction as F
from itertools import combinations
from math import gcd, isqrt
import hashlib
import json
from pathlib import Path
import sys

sys.stdout.reconfigure(newline="\n")
GATES=0


def gate(ok,label):
    global GATES
    GATES+=1
    if not ok:
        raise ArithmeticError(label)


def prime(p):
    return p>=2 and all(p%d for d in range(2,isqrt(p)+1))


def det(P,R,S):
    return (R[0]-P[0])*(S[1]-P[1])-(S[0]-P[0])*(R[1]-P[1])


def invs(p):
    values=[0]*p
    values[1]=1
    for x in range(2,p):
        values[x]=(-(p//x)*values[p%x])%p
    for x in range(1,p):
        gate(x*values[x]%p==1,"literal inverse recurrence")
    return values


def board(p):
    f=invs(p)
    return [(0,1),(1,0)]+[(x,f[x]) for x in range(2,p)]


def anchored_lines(p,B):
    """Actual integer primitive slopes, without modular quadratic tests."""
    groups=defaultdict(list)
    for x,y in B[2:]:
        d=gcd(y-1,x)
        groups[( (y-1)//d,x//d)].append((x,y))
    answer={}
    for (a,b),points in groups.items():
        gate(len(points)<=2,"conic bounds anchored line multiplicity")
        if len(points)==2:
            answer[a,b]=tuple(sorted(points))
    return answer


def packets(p):
    """Complete primitive integer packet, not a search through board points."""
    answer={}
    records=[]
    for h in range(5,p):
        if (p-1)%h:
            continue
        a=(p-1)//h
        for k in range(2,(h-1)//2+1):
            l=h-k
            for c in range(1,k):
                num=h+c*p
                if num%(k*l):
                    continue
                b=num//(k*l)
                if gcd(a,b)!=1 or b*l>p-1:
                    continue
                points=((b*k,a*k+1),(b*l,a*l+1))
                gate((a,b) not in answer,"primitive packet has one owner")
                answer[a,b]=points
                records.append([a,b,h,k,l,c])
                gate(1<=b<2*a and 5*a<=p-1,"strict slope and short-cofactor bounds")
    return answer,records


def character_lines(p):
    """Quadratic discriminant plus the standard-root interval sidecar."""
    squares={s*s%p:s for s in range(1,p)}
    answer={}
    for a in range(1,(p-1)//5+1):
        if (p-1)%a:
            continue
        h=(p-1)//a
        for b in range(1,2*a):
            if gcd(a,b)!=1:
                continue
            delta=b*(b+4*a)%p
            if not delta or delta not in squares:
                continue
            s=squares[delta]
            den=pow(2*a*b,-1,p)
            roots=sorted({(-b+s)*den%p,(-b-s)*den%p})
            gate(len(roots)==2,"nonzero discriminant gives distinct roots")
            H=min(h-1,(p-1)//b)
            if not all(1<=r<=H for r in roots):
                continue
            k,l=roots
            gate(k+l==h,"bounded roots force the exact unwrapped sum")
            answer[a,b]=((b*k,a*k+1),(b*l,a*l+1))
    return answer


records=[]
# Complete small hostile/minimality universe; the larger cases are prescribed
# arithmetic templates or the gauge/extra-line controls, not a prime census.
small=[p for p in range(3,44,2) if prime(p)]
controls=small+[61,97,113,197]
for p in controls:
    B=board(p)
    actual=anchored_lines(p,B)
    arithmetic,pack=packets(p)
    character=character_lines(p)
    gate(actual==arithmetic==character,"three independent exact line representations")
    predicted=set()
    for R,S in actual.values():
        predicted.add(tuple(sorted(((0,1),R,S))))
        predicted.add(tuple(sorted(((1,0),(R[1],R[0]),(S[1],S[0])))))
    literal=set()
    for triple in combinations(B,3):
        if det(*triple)==0:
            literal.add(tuple(sorted(triple)))
    gate(literal==predicted,"complete original integer-triangle census")
    gate(len(literal)==2*len(actual),"all actual defects are disjoint transpose pairs")
    chi=0 if p==5 else (1 if pow(5,(p-1)//2,p)==1 else -1)
    if chi==1:
        gate((1,1) in actual,"discriminant-five obstruction is a native line")
    records.append({"p":p,"chi5":chi,"triples":len(literal),"packets":pack})

by_p={r["p"]:r for r in records}
gate(any(r["chi5"]==-1 and r["triples"] for r in records),"negative-character hostile exists")
gate(min(r["p"] for r in records if r["chi5"]==-1 and r["triples"])==37,"minimal negative-five-character hostile, all smaller primes checked")
gate(all(by_p[p]["triples"]==0 for p in (3,5,7,13,17,23)),"retained genuine positive controls")
gate(by_p[11]["triples"]==2 and by_p[37]["triples"]==2,"two minimal native obstructions")
gate(by_p[113]["triples"]==4,"a construction's two forced triples need not be all triples")
gate([12,11,5,2,3,1] in by_p[61]["packets"],"short-cofactor lower bound five is attained")

# Exact hostile to dropping the lift interval, in the eligible slope range.
p,a,b=13,2,1
h=(p-1)//a
H=min(h-1,(p-1)//b)
roots=[k for k in range(p) if (a*b*k*k+b*k-1)%p==0]
gate(roots==[7,12] and H==5 and b*(b+4*a)==9,"square discriminant with absent native roots")
gate(det((0,1),(7,2),(12,12))==65,"modular zero is a nonzero native determinant")
gate(by_p[13]["triples"]==0,"character-only test falsely rejects a successful repair")

# Repeated root at p5 creates no pair of distinct conic points.
gate([x for x in range(5) if (x*x+x-1)%5==0]==[2],"ramified golden-ratio root")


def mul_linear(u,v):
    return [u[0]*v[0],u[0]*v[1]+u[1]*v[0],u[1]*v[1]]


def linear_nonnegative(u):
    return u[0]>=0 and u[1]>=0


templates=[(37,360,18,8,10,6),(43,70,7,2,5,1),(97,120,8,3,5,1)]
template_records=[]
for residue,modulus,h,k,l,c in templates:
    # p=residue+modulus*r. Identity and box checks are coefficientwise for all
    # integer r>=0, and do not depend on which progression values are prime.
    A=[F(residue-1,h),F(modulus,h)]
    D=[F(h+c*residue,k*l),F(c*modulus,k*l)]
    gate(all(t.denominator==1 for t in A+D),"integral affine construction coefficients")
    gate(h==k+l and 2<=k<l and 1<=c<k,"general construction packet assumptions")
    gate(gcd(residue,modulus)==1,"Dirichlet progression is primitive")
    gate(modulus%5==0 and residue%5 in (2,3),"entire prime progression has chi5 negative")
    points=[]
    quotients=[]
    for j in (k,l):
        X=[j*z for z in D]
        Y=[j*A[0]+1,j*A[1]]
        gate(linear_nonnegative([X[0]-2,X[1]]) and linear_nonnegative([Y[0]-2,Y[1]]),"coordinates stay above removed conic point")
        gate(linear_nonnegative([residue-1-X[0],modulus-X[1]]) and linear_nonnegative([residue-1-Y[0],modulus-Y[1]]),"all native coordinates stay in the integer box")
        product=mul_linear(X,Y)
        product[0]-=1
        q1=product[2]/modulus
        q0=product[0]/residue
        gate(product==mul_linear([residue,modulus],[q0,q1]),"exact all-parameter reciprocal product identity")
        gate(q0.denominator==q1.denominator==1,"integral reciprocal-product quotient")
        gate(mul_linear(D,[Y[0]-1,Y[1]])==mul_linear(A,X),"exact all-parameter anchored line identity")
        points.append([[int(t) for t in X],[int(t) for t in Y]])
        quotients.append([int(q0),int(q1)])
    template_records.append({"residue":residue,"modulus":modulus,"h":h,"k":k,"l":l,"c":c,"a":[int(t) for t in A],"b":[int(t) for t in D],"points":points,"product_quotients":quotients})

# The r=1 first h7-template has a nonprimitive displayed slope. Primitive
# normalization changes (h,k,l,c) and must precede the bijective count.
gate(gcd(16,12)==4,"p113 construction has nonprimitive displayed slope")
gate([4,3,28,8,20,4] in by_p[113]["packets"],"primitive packet restores the factor-four gauge")

payload={"scope":"completed inverse permutation, only output swap0 with1, actual integer triples","controls":records,"templates":template_records}
raw=json.dumps(payload,sort_keys=True,indent=2)+"\n"
base=Path(__file__).resolve()
outdir=base.parent.parent/"05-knowledge"/"results" if base.parent.name=="04-computation" else base.parent
outdir.mkdir(parents=True,exist_ok=True)
(outdir/(base.stem+"_certificate.json")).write_bytes(raw.encode())
print("PROOF CANDIDATE: exact primitive native-line packet and root-interval criterion; independent audit pending.")
print("Controls: all odd primes<=43 plus61,97,113,197; three representations and literal full triangle sets agree.")
for r in records:
    print("CONTROL",json.dumps(r,separators=(",",":")))
for r in template_records:
    print("TEMPLATE",json.dumps(r,separators=(",",":")))
print("Native lift hostile p13:(a,b)=(2,1),Delta9,roots7/12,H5,determinant65; actual repair succeeds.")
print("All three primitive progressions have chi5=-1 and two forced native triples; infinity uses CITED Dirichlet.")
print("The original one-defect graph is repaired to zero defects or changed to an even positive number.")
print("No all-prime repair, general zero-swap classification, or2p-point construction is inferred.")
print("CERTIFICATE_SHA256",hashlib.sha256(raw.encode()).hexdigest())
print("PASS",GATES,"always-active exact gates.")
