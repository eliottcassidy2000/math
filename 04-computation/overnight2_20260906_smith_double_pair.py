#!/usr/bin/env python3
"""Universal finite-symbolic certificate for the p=3 double-pair cluster.

Two independent polynomial determinant algorithms verify every one of the
923 nonempty minors. Termwise affine dominance and thirteen unit-leading
witness minors prove all positive scale/depth values, with arbitrary units.
Independent integer DVR and literal-minor controls check the specialization.
No project imports, floating point, or optimization-disabled assertions.
"""
from fractions import Fraction
from hashlib import sha256
from itertools import combinations, permutations, product
from math import comb
import json
import sys

sys.stdout.reconfigure(newline="\n")
GATES=0


def require(condition,label):
    global GATES
    GATES+=1
    if not condition:
        raise RuntimeError(label)


def add(a, b, sign=1):
    c = dict(a)
    for ij, value in b.items():
        c[ij] = c.get(ij, 0) + sign*value
        if not c[ij]:
            del c[ij]
    return c


def mul(a, b):
    c = {}
    for (i,j), x in a.items():
        for (k,l), y in b.items():
            ij = (i+k,j+l)
            c[ij] = c.get(ij,0)+x*y
    return {ij:v for ij,v in c.items() if v}


def vp(value):
    if not value:
        raise ValueError("nonzero valuation required")
    result = 0
    while value % 3 == 0:
        value //= 3
        result += 1
    return result


def fraction_vp(value):
    value=Fraction(value)
    return vp(abs(value.numerator))-vp(value.denominator)


def residual():
    # After clearing the value and first derivative at zero, degrees 2..7.
    rows = []
    for node in ("one","A","oneplusB"):
        for order in (0,1):
            row=[]
            for degree in range(2,8):
                q=degree-order
                factor=degree if order else 1
                if node=="one": p={(0,0):factor}
                elif node=="A": p={(q,0):factor}
                else: p={(0,j):factor*comb(q,j) for j in range(q+1)}
                row.append(p)
            rows.append(row)
    return rows


def polynomial_minors():
    matrix=residual()
    previous={((),()):{(0,0):1}}
    all_ranks=[]
    for k in range(1,7):
        current={}
        for rows in combinations(range(6),k):
            for cols in combinations(range(6),k):
                polynomial={}
                for j,col in enumerate(cols):
                    prior=previous[rows[:-1],cols[:j]+cols[j+1:]]
                    polynomial=add(polynomial,mul(matrix[rows[-1]][col],prior),(-1)**(k-1+j))
                current[rows,cols]=polynomial
        all_ranks.append(current)
        previous=current
    return all_ranks


def independent_minor(rows,cols):
    # Literal Leibniz expansion; no matrix polynomial or Laplace cache.
    answer={}
    for perm in permutations(cols):
        inversions=sum(perm[i]>perm[j] for i in range(len(perm))
                       for j in range(i+1,len(perm)))
        coefficient=(-1)**inversions
        adegree=bdegree=0
        for row,col in zip(rows,perm):
            degree=col+2
            order=row%2
            if order:
                coefficient*=degree
            if row//2==1:
                adegree+=degree-order
            elif row//2==2:
                bdegree+=degree-order
        for j in range(bdegree+1):
            ij=(adegree,j)
            answer[ij]=answer.get(ij,0)+coefficient*comb(bdegree,j)
    return {ij:v for ij,v in answer.items() if v}


def translated(form):
    w,i,j,c=form
    return w,i,j,w+i+j+c


def dominates(a,b):
    return all(x<=y for x,y in zip(translated(a),translated(b)))


def pareto(forms):
    return sorted(a for a in forms if not any(b!=a and dominates(b,a) for b in forms))


FRONTIERS=(
    ((1,0,0,0),),
    ((3,0,1,1),(3,1,0,1),(4,0,0,0)),
    ((6,1,1,1),(7,0,1,0),(7,1,0,0)),
    ((11,1,1,0),(12,0,4,0),(12,4,0,0)),
    ((17,1,4,0),(17,4,1,0)),
    ((24,4,4,0),),
)


def certificate(all_ranks):
    manifest=[]
    term_count=0
    for k,minors in enumerate(all_ranks,1):
        terms={}
        robust={}
        for (rows,cols),poly in minors.items():
            require(poly==independent_minor(rows,cols),"independent polynomial determinant")
            manifest.append([rows,cols,sorted((i,j,v) for (i,j),v in poly.items())])
            w=sum(c+2 for c in cols)-sum(r%2 for r in rows)
            for (i,j),coefficient in poly.items():
                term_count+=1
                form=(w,i,j,vp(coefficient))
                terms.setdefault(form,(rows,cols,i,j,coefficient))
                require(any(dominates(f,form) for f in FRONTIERS[k-1]),
                        "every polynomial term obeys a universal lower bound")
                # Its unique minimum survives every positive integer depth.
                if all((a>=i and b>=j and
                        (a-i)+(b-j)+vp(v)-vp(coefficient)>0)
                       for (a,b),v in poly.items() if (a,b)!=(i,j)):
                    robust.setdefault(form,(rows,cols,i,j,coefficient))
        frontier=pareto(terms)
        robust_frontier=pareto(robust)
        require(tuple(frontier)==FRONTIERS[k-1],"literal lower frontier")
        require(tuple(robust_frontier)==FRONTIERS[k-1],"every frontier has robust witness")
        print("RANK",k,"MINORS",len(minors),"TERMS",len(terms),flush=True)
        print("LOWER",frontier,flush=True)
        print("ROBUST",robust_frontier,flush=True)
        for form in robust_frontier:
            print("WITNESS",form,robust[form],flush=True)
    digest=sha256(json.dumps(manifest,separators=(",",":"),sort_keys=True).encode()).hexdigest()
    print("SYMBOLIC_MINORS",len(manifest),"NONZERO_COEFFICIENTS",term_count)
    print("SYMBOLIC_SHA256",digest)
    return digest


def formula(e,d,r):
    s,t=sorted((d,r))
    divisors=(0,e,min(4*e,3*e+s+1),min(7*e+s,6*e+s+t+1),
              min(11*e+s+t,12*e+4*s),17*e+4*s+t,24*e+4*s+4*t)
    return (0,0)+tuple(b-a for a,b in zip(divisors,divisors[1:]))


def observer(nodes):
    return [[x**q if not order else (q*x**(q-1) if q else 0)
             for q in range(2*len(nodes))] for x in nodes for order in (0,1)]


def rational_smith(matrix):
    matrix=[[Fraction(v) for v in row] for row in matrix]
    answer=[]
    for k in range(len(matrix)):
        i,j=min(((i,j) for i in range(k,len(matrix)) for j in range(k,len(matrix))
                 if matrix[i][j]),key=lambda ij:fraction_vp(matrix[ij[0]][ij[1]]))
        matrix[k],matrix[i]=matrix[i],matrix[k]
        for row in matrix:
            row[k],row[j]=row[j],row[k]
        pivot=matrix[k][k]
        answer.append(fraction_vp(pivot))
        for i in range(k+1,len(matrix)):
            multiplier=matrix[i][k]/pivot
            require(not multiplier or fraction_vp(multiplier)>=0,"integral DVR elimination")
            matrix[i]=[a-multiplier*b for a,b in zip(matrix[i],matrix[k])]
    require(answer==sorted(answer),"rational invariants ordered")
    return tuple(answer)


def modular_smith(nodes):
    determinant_v=4*sum(vp(abs(y-x)) for i,x in enumerate(nodes) for y in nodes[i+1:])
    precision=determinant_v+1
    modulus=3**precision
    matrix=[[x%modulus for x in row] for row in observer(nodes)]
    answer=[]
    def val(x):
        return vp(x) if x else precision
    for k in range(len(matrix)):
        i,j=min(((i,j) for i in range(k,len(matrix)) for j in range(k,len(matrix))),
                key=lambda ij:val(matrix[ij[0]][ij[1]]))
        matrix[k],matrix[i]=matrix[i],matrix[k]
        for row in matrix:
            row[k],row[j]=row[j],row[k]
        v=val(matrix[k][k])
        require(v<precision,"precision exceeds every pivot valuation")
        answer.append(v)
        power=3**v
        reduced=modulus//power
        inverse=pow(matrix[k][k]//power,-1,reduced)
        for i in range(k+1,len(matrix)):
            require(matrix[i][k]%power==0,"modular elimination divisibility")
            multiplier=(matrix[i][k]//power)*inverse%reduced
            matrix[i]=[(a-multiplier*b)%modulus for a,b in zip(matrix[i],matrix[k])]
            require(matrix[i][k]==0,"modular elimination cleared")
    require(sum(answer)==determinant_v,"Vandermonde determinant mass")
    return tuple(answer)


def bareiss(matrix):
    matrix=[list(row) for row in matrix]
    sign=previous=1
    for k in range(len(matrix)-1):
        i=next((i for i in range(k,len(matrix)) if matrix[i][k]),None)
        if i is None:
            return 0
        if i!=k:
            matrix[k],matrix[i]=matrix[i],matrix[k]
            sign*=-1
        pivot=matrix[k][k]
        for i in range(k+1,len(matrix)):
            for j in range(k+1,len(matrix)):
                numerator=pivot*matrix[i][j]-matrix[i][k]*matrix[k][j]
                require(numerator%previous==0,"Bareiss exact division")
                matrix[i][j]=numerator//previous
            matrix[i][k]=0
        previous=pivot
    return sign*matrix[-1][-1]


def literal_all_minors(nodes):
    # Full original 8x8 observer: does not assume the two cleared unit pivots.
    matrix=observer(nodes)
    values=[0]
    for k in range(1,9):
        least=None
        for rows in combinations(range(8),k):
            for cols in combinations(range(8),k):
                determinant=bareiss([[matrix[i][j] for j in cols] for i in rows])
                if determinant:
                    v=vp(abs(determinant))
                    least=v if least is None else min(least,v)
        require(least is not None,"nonzero minor each rank")
        values.append(least)
    return tuple(b-a for a,b in zip(values,values[1:]))


def numeric_audit():
    manifest=[]
    units=(1,2,4,5,7,8,10,11,-1,-2)
    for e,d,r in product(range(1,5),repeat=3):
        for a,b in product(units,repeat=2):
            nodes=tuple(3**e*z for z in (0,1,3**d*a,1+3**r*b))
            result=modular_smith(nodes)
            require(result==formula(e,d,r),"all unit lifts agree with universal formula")
            manifest.append([e,d,r,a,b,result])
    for e,d,r in product(range(1,5),range(1,11),range(1,11)):
        nodes=tuple(17-2*3**e*z for z in (0,1,-3**d,1+2*3**r))
        result=rational_smith(observer(nodes))
        require(result==formula(e,d,r),"independent rational deep and translated controls")
    literal_cases=((1,1,1,1,1),(1,1,5,-2,4),(4,1,1,2,-1),(2,3,2,7,-2))
    for e,d,r,a,b in literal_cases:
        nodes=tuple(3**e*z for z in (0,1,3**d*a,1+3**r*b))
        require(literal_all_minors(nodes)==formula(e,d,r),"every original eight-by-eight minor")
    # Optional outer-depth-zero boundary is independently CRT, not used above.
    for d,r in product(range(1,7),repeat=2):
        nodes=(0,1,3**d,1+2*3**r)
        expected=tuple(sorted((0,0,d,3*d,0,0,r,3*r)))
        require(modular_smith(nodes)==expected==formula(0,d,r),"unit-separated CRT boundary")
    shallow=tuple(sorted((0,0,1,3,4,9,9,22)))
    require(formula(1,1,5)!=shallow,"naive sorted shallow continuation rejected")
    require(sum(formula(1,1,3))==sum(formula(1,2,2)),"same determinant hostile mass")
    require(formula(1,1,3)!=formula(1,2,2),"determinant alone loses cluster split")
    digest=sha256(json.dumps(manifest,separators=(",",":")).encode()).hexdigest()
    print("NUMERIC_UNIVERSE unit_lifts",len(manifest),"independent_rational",400,
          "full_literal_all_minor_cases",len(literal_cases),"CRT_boundary",36)
    print("HOSTILE_shallow",formula(1,1,5),"naive",shallow)
    print("HOSTILE_scale",formula(4,1,1))
    print("HOSTILE_same_determinant",formula(1,1,3),formula(1,2,2))
    print("NUMERIC_SHA256",digest)


def main():
    all_ranks=polynomial_minors()
    certificate(all_ranks)
    numeric_audit()
    print("PASS universal symbolic certificate and independent exact specializations")
    print("GATES",GATES)


if __name__=="__main__":
    main()
