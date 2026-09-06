"""Independent SymPy-domain minor and integer-Smith audit; no producer imports."""
from itertools import combinations, permutations
import sys
import sympy as s
from sympy.matrices.normalforms import smith_normal_form
from sympy.polys.matrices import DomainMatrix

sys.stdout.reconfigure(newline="\n")
checks=0


def need(ok,label):
    global checks
    checks+=1
    if not ok:
        raise RuntimeError(label)


def valuation(x):
    x=abs(int(x))
    if not x:
        raise ValueError('zero valuation')
    v=0
    while x%3==0:
        x//=3
        v+=1
    return v


def observer(nodes):
    return s.Matrix([[x**j if d==0 else (j*x**(j-1) if j else 0)
                      for j in range(2*len(nodes))] for x in nodes for d in (0,1)])


def spectrum(nodes):
    d=smith_normal_form(observer(nodes),domain=s.ZZ)
    return tuple(valuation(d[i,i]) for i in range(d.rows))


def main():
    A,C=s.symbols('A C')
    ring=s.ZZ.poly_ring(A,C)
    matrix=s.Matrix([[u**j if order==0 else j*u**(j-1)
                      for j in range(2,8)] for u in (1,A,C) for order in (0,1)])
    floors=(2,8,15,27,42,64)
    survivors=[]
    semantic=[]
    count=coefficients=0
    for size in range(1,7):
        found=[]
        for rows in combinations(range(6),size):
            for cols in combinations(range(6),size):
                # SymPy's domain determinant is a third algebra implementation.
                d=DomainMatrix.from_Matrix(matrix.extract(rows,cols)).convert_to(ring).det()
                terms={ij:int(v) for ij,v in d.to_dict().items()}
                semantic.append([rows,cols,sorted((i,j,c) for (i,j),c in terms.items())])
                count+=1;coefficients+=len(terms)
                weight=2*(sum(j+2 for j in cols)-sum(i%2 for i in rows))
                residue={}
                for (i,j),c in terms.items():
                    v=valuation(c)
                    need(weight+i+2*j+v>=floors[size-1],'coefficient floor')
                    if weight+i+2*j+v==floors[size-1]:
                        residue[i,j]=(c//3**v)%3
                if residue:
                    found.append((rows,tuple(j+2 for j in cols),residue))
        survivors.append(found)
        print('INDEPENDENT_RANK',size,'FLOOR',floors[size-1],'SURVIVORS',found)
    expected=[
        [((1,),(2,),{(0,0):2})],
        [((0,1),(2,3),{(0,0):1}),((1,3),(2,3),{(1,0):1})],
        [((0,1,3),(2,3,4),{(1,0):2})],
        [((0,1,3,5),(2,3,4,5),{(2,1):1,(3,1):1})],
        [((0,1,2,3,5),(2,3,4,5,6),{(6,1):2})],
        [((0,1,2,3,4,5),(2,3,4,5,6,7),{(8,4):1})],
    ]
    need(survivors==expected,'complete residue table')
    need(count==923 and coefficients==5587,'complete minor universe')
    ceiling=s.factor(matrix.extract((0,1,2,3),(0,1,2,3)).det())
    need(ceiling==A**4*(A-1)**4,'uniform ceiling witness')
    critical=s.factor(matrix.extract((0,1,3,5),(0,1,2,3)).det())
    need(s.expand(critical+2*A*C*(A-1)*(C-1)*(A-C)*(10*A*C-5*A-5*C+3))==0,
         'critical cancellation factor')
    print('MINORS',count,'COEFFICIENTS',coefficients,'SYMBOLIC_RECONSTRUCTION_PASS')
    controls=0
    for a in (-11,-4,1,2,8,17):
        for b in (-7,-1,1,5):
            nodes=(0,9,27*a,81*b)
            kappa=int(a%3==2)
            actual=spectrum(nodes)
            need(actual==(0,0,2,6,7,12+kappa,15-kappa,22),'integer Smith full-family control')
            X=s.symbols('X')
            F=s.prod(X-v for v in nodes)
            derivative=s.diff(F,X)
            second=s.diff(F,X,2)
            hermite=[]
            for v in nodes:
                S=valuation(derivative.subs(X,v))
                d2=second.subs(X,v)
                hermite.append(max(2*S,3*S-valuation(d2)) if d2 else 2*S)
            need(max(hermite)==22,'incoming Hermite largest-factor identity')
            controls+=1
    trees={}
    perms=list(permutations(range(4)))
    for tail in combinations(range(3,81,3),3):
        nodes=(0,)+tail
        metrics=[[0 if i==j else valuation(nodes[i]-nodes[j]) for j in range(4)] for i in range(4)]
        key=min(tuple(metrics[p[i]][p[j]] for i in range(4) for j in range(i+1,4)) for p in perms)
        actual=spectrum(nodes)
        if key in trees:
            need(trees[key]==actual,'no metric collision below diameter81')
        else:
            trees[key]=actual
    need(len(trees)==11,'eleven metric trees')
    x=spectrum((0,9,27,81));y=spectrum((0,9,54,81))
    need(sum(x)==sum(y)==64 and x[-1]==y[-1]==22,'same determinant and largest invariant')
    need(sum(min(13,v) for v in x)==53 and sum(min(13,v) for v in y)==54,'different actual kernels')
    print('INDEPENDENT_INTEGER_SMITH',controls,'UNIT_CONTROLS',2600,'DIAMETER_CONTROLS',len(trees),'TREES')
    print('PASS full symbolic residue proof, minimal integer diameter81, kernel distinction')
    print('CHECKS',checks)


if __name__=='__main__':
    main()
