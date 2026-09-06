#!/usr/bin/env python3
"""Exact unit-sensitive ternary (3,1) cluster certificate.
Independent polynomial determinant algorithms, all weighted minor minima,
actual attaining witnesses and exact integer controls; no project imports.
Arithmetic helpers inherited from the double-pair certificate.
"""
from fractions import Fraction
from hashlib import sha256
from itertools import combinations, permutations, product
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
    for node in ("one","A","C"):
        for order in (0,1):
            row=[]
            for degree in range(2,8):
                q=degree-order
                factor=degree if order else 1
                if node=="one": p={(0,0):factor}
                elif node=="A": p={(q,0):factor}
                else: p={(0,q):factor}
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
        ij=(adegree,bdegree)
        answer[ij]=answer.get(ij,0)+coefficient
    return {ij:v for ij,v in answer.items() if v}


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


LOW=(2,8,15,27,42,64)
EXPECTED=(
    [(((1,),(0,)),{(0,0):2})],
    [(((0,1),(0,1)),{(0,0):1}),(((1,3),(0,1)),{(1,0):1})],
    [(((0,1,3),(0,1,2)),{(1,0):2})],
    [(((0,1,3,5),(0,1,2,3)),{(2,1):1,(3,1):1})],
    [(((0,1,2,3,5),(0,1,2,3,4)),{(6,1):2})],
    [(((0,1,2,3,4,5),(0,1,2,3,4,5)),{(8,4):1})],
)


def weight(rows,cols):
    return sum(c+2 for c in cols)-sum(r%2 for r in rows)


def normalized_residue(poly,w,level):
    residue={}
    for (i,j),c in poly.items():
        value=2*w+i+2*j+vp(c)
        require(value>=level,"all weighted polynomial coefficients above bound")
        if value==level:
            residue[i,j]=(c//3**vp(c))%3
    return residue


def symbolic_certificate():
    ranks=polynomial_minors()
    manifest=[]
    coefficients=0
    for k,minors in enumerate(ranks):
        survivors=[]
        for (rows,cols),poly in minors.items():
            require(poly==independent_minor(rows,cols),"independent coefficient dictionaries")
            coefficients+=len(poly)
            manifest.append([rows,cols,sorted((i,j,c) for (i,j),c in poly.items())])
            residue=normalized_residue(poly,weight(rows,cols),LOW[k])
            if residue:
                survivors.append(((rows,cols),residue))
        require(survivors==EXPECTED[k],"complete normalized minor residue table")
        print("RANK",k+1,"MINORS",len(minors),"FLOOR",LOW[k],"SURVIVORS",survivors)
    rows,cols=(0,1,2,3),(0,1,2,3)
    ceiling=normalized_residue(ranks[3][rows,cols],weight(rows,cols),28)
    require(ceiling=={(4,0):1},"unit witness caps exceptional valuation at 28")
    print("EXCEPTIONAL_CEILING",28,"ROWS",rows,"DEGREES",tuple(c+2 for c in cols),"RESIDUE",ceiling)
    for a,b in product((1,2),repeat=2):
        critical=sum(c*a**i*b**j for (i,j),c in EXPECTED[3][0][1].items())%3
        require((critical==0)==(a==2),"critical residue criterion for all units")
    digest=sha256(json.dumps(manifest,separators=(",",":")).encode()).hexdigest()
    print("SYMBOLIC_UNIVERSE",len(manifest),"MINORS",coefficients,"COEFFICIENTS")
    print("SYMBOLIC_SHA256",digest)


def spectrum(a):
    kappa=int(a%3==2)
    return (0,0,2,6,7,12+kappa,15-kappa,22)


def canonical_tree(nodes):
    matrix=[[vp(abs(a-b)) if a!=b else 0 for b in nodes] for a in nodes]
    return min(tuple(matrix[order[i]][order[j]] for i in range(4) for j in range(i+1,4))
               for order in permutations(range(4)))


def numeric_audit():
    units=tuple(x for x in range(-40,41) if x%3)
    manifest=[]
    for a,b in product(units,repeat=2):
        result=modular_smith((0,9,27*a,81*b))
        require(result==spectrum(a),"arbitrary unit specializations")
        manifest.append([a,b,result])
    independent_units=(-10,-5,-2,-1,1,2,4,5,7,8,10,11)
    for a,b in product(independent_units,repeat=2):
        nodes=tuple(17-2*x for x in (0,9,27*a,81*b))
        require(rational_smith(observer(nodes))==spectrum(a),"independent rational translated unit controls")
    controls=((1,1),(2,1),(-1,-2),(8,7))
    for a,b in controls:
        require(literal_all_minors((0,9,27*a,81*b))==spectrum(a),"every original eight-by-eight minor")
    first,second=(0,9,27,81),(0,9,54,81)
    require(canonical_tree(first)==canonical_tree(second),"identical full metric trees")
    require(spectrum(1)!=spectrum(2),"metric-only implication refuted")
    require(sum(spectrum(1))==sum(spectrum(2))==64,"determinant mass retained")
    deep_cancellations=[]
    for a,expected in ((2,1),(8,2),(17,3)):
        value=vp(85*a-14)
        require(value==expected,"critical minor can cancel deeper")
        require(modular_smith((0,9,27*a,81))==spectrum(2),"competing minor caps deeper cancellation")
        deep_cancellations.append((a,value,spectrum(a)))
    # Every translated integer configuration of diameter <81 and positive
    # common depth has this form. Common depth zero is handled analytically
    # by integral CRT and the independently proved <=3-node metric law.
    trees={}
    configurations=0
    for tail in combinations(range(3,81,3),3):
        nodes=(0,)+tail
        tree=canonical_tree(nodes)
        result=modular_smith(nodes)
        require(tree not in trees or trees[tree]==result,"no smaller-diameter positive-depth counterexample")
        trees[tree]=result
        configurations+=1
    require(configurations==2600,"complete positive-depth diameter universe")
    print("NUMERIC_UNIVERSE",len(manifest),"UNIT_PAIRS",144,"RATIONAL_CONTROLS",len(controls),"ALL_MINOR_CONTROLS")
    print("MINIMAL_DIAMETER_CHECK",configurations,"CONFIGURATIONS",len(trees),"TREES","NO_COUNTEREXAMPLE_BELOW",81)
    print("HOSTILE",first,spectrum(1),second,spectrum(2))
    print("HOSTILE_TREE",canonical_tree(first))
    print("CAPPED_DEEP_CANCELLATIONS",deep_cancellations)
    print("NUMERIC_SHA256",sha256(json.dumps(manifest,separators=(",",":")).encode()).hexdigest())


def main():
    symbolic_certificate()
    numeric_audit()
    print("PASS universal one-residue criterion and exact minimal-diameter controls")
    print("GATES",GATES)


if __name__=="__main__":
    main()


