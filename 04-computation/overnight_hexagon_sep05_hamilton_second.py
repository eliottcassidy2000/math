"""Exact controls for the n>=16 second Hamilton weight theorem candidate."""
from itertools import combinations, permutations
from math import comb, factorial
import sympy as s

gates=0


def check(good,label):
    global gates
    if not good:
        raise RuntimeError(label)
    gates+=1


def identity(left,right,label):
    check(s.cancel(left-right)==0,label)


def signed_matching(n,m):
    numer=sum((-4)**j*comb(m,j)*factorial(n-j-1) for j in range(m+1))
    check(numer%2==0,"matching signed-sum integrality")
    return numer//2


def pair(u,v):
    return tuple(sorted((u,v)))


motifs={
    "P3":{(0,1),(1,2)},
    "P4":{(0,1),(1,2),(2,3)},
    "C3":{(0,1),(1,2),(0,2)},
    "C4":{(0,1),(1,2),(2,3),(0,3)},
    "C5":{(0,1),(1,2),(2,3),(3,4),(0,4)},
}


def motif_formula(n,name):
    e,a,c,d=(factorial(n-j) for j in (2,3,4,5))
    return {"P3":2*e-2*a,"P4":3*e-8*a+4*c,"C3":3*e-6*a,
            "C4":4*e-16*a+16*c,"C5":5*e-30*a+60*c-40*d}[name]


total_cycles=0
for n in range(3,10):
    weights=[0]*(n//2+1)
    finite={name:0 for name in motifs} if n>=6 else {}
    total=0
    for p in permutations(range(1,n)):
        if p[0]>p[-1]:
            continue
        word=(0,)+p
        edges={pair(word[j],word[(j+1)%n]) for j in range(n)}
        total+=1
        parity=0
        for m in range(n//2+1):
            if m:
                parity^=(2*m-2,2*m-1) in edges
            weights[m]+=parity
        for name in finite:
            finite[name]+=len(edges&motifs[name])%2
    check(total==factorial(n-1)//2,"literal Hamilton universe")
    for m,w in enumerate(weights):
        check(total-2*w==signed_matching(n,m),"literal matching versus exact integral expansion")
    for name,w in finite.items():
        check(w==motif_formula(n,name),"literal finite-motif weight")
    if n==6:
        check(finite["C5"]==20<factorial(4),"negative C5 plus apex hostile")
    total_cycles+=total
    print(f"n={n} cycles={total} matching_weights={weights} motif_controls={len(finite)}")

for n in range(4,101):
    for m in range(n//2):
        check(signed_matching(n,m)-signed_matching(n,m+1)==4*signed_matching(n-1,m),
              "exact contraction recurrence")
for n in range(3,101):
    for m in range(n//2+1):
        if n-m-1>=4:
            check(signed_matching(n,m)>0,"scoped integral positivity")
check(signed_matching(3,1)==-1,"negative signed-sum hostile n3")
check(signed_matching(5,1)==0,"zero signed-sum hostile n5")
check(signed_matching(6,3)==-4,"negative perfect-matching signed-sum hostile n6")

n,t,q=s.symbols("n t q")
F=q*((n-2)*(n-3)-2*q)
target=2*(n-4)**2*(n-3)
lo=n**3-26*n**2+193*n-444
hi=(n-4)*(n**3-22*n**2+124*n-200)/8
identity(F.subs(q,3*(n-5))-target,lo,"R2 lower endpoint")
identity(F.subs(q,(n-2)**2/4)-target,hi,"R2 upper real endpoint")
for poly in (lo,hi,n*n-9*n+22,n*n-11*n+32,3*n**3-53*n*n+320*n-664):
    check(all(c>0 for c in s.Poly(poly.subs(n,16+t),t).all_coeffs()),
          "all-height positive shifted coefficients")
E=(n-2)*(n-3)*(n-4)
A=(n-3)*(n-4)
C=n-4
identity(3*E-8*A+4*C-(2*E-4*A),C*(n*n-9*n+22),"P4 excess")
identity(4*E-16*A+16*C-(2*E-4*A),2*C*(n*n-11*n+32),"C4 excess")
identity(5*E-30*A+60*C-40-(2*E-4*A),3*n**3-53*n*n+320*n-664,"C5 excess")
x=s.symbols("x")
identity(s.diff(s.log((4+x)/(4-x))-x/2,x), x*x/(2*(16-x*x)),
         "integral reflected-density logarithmic derivative")


def r_values(n,edges):
    rows=[set() for _ in range(n)]
    for u,v in edges:
        rows[u].add(v);rows[v].add(u)
    return [len((rows[u]^rows[v])-{u,v}) for u in range(n) for v in range(u+1,n)]


def switch(n,edges,selected):
    return {pair(u,v) for u in range(n) for v in range(u+1,n)
            if ((pair(u,v) in edges)^((u in selected)!=(v in selected)))}


n=16
families=list(motifs.values())+[{(2*j,2*j+1) for j in range(m)} for m in range(9)]
for edges in families:
    check(max(r_values(n,edges))<=2,"all structural survivors obey R2")
    for selected in ({0},{1,3},{0,2,4,7},{0,1,2,3,4,5,6,7}):
        vals=r_values(n,switch(n,edges,selected))
        check(set(vals)<={0,1,2,n-4,n-3,n-2},"switching preserves projective R2")
for name,edges in {
        "P5":{(0,1),(1,2),(2,3),(3,4)},
        "C6":{(0,1),(1,2),(2,3),(3,4),(4,5),(0,5)},
        "P3_plus_edge":{(0,1),(1,2),(3,4)},
        "three_star":{(0,1),(0,2),(0,3)}}.items():
    check(max(r_values(n,edges))>2,"structural hostile "+name)
for a in range(3,25):
    for b in range(3,25):
        if (a-2)*(b-2)<=4:
            check(a+b<=9,"low/high partition capacity")
check(2*(13-5)>12,"degree-four common-neighbour capacity")
check(13-4>6,"degree-three common-neighbour capacity")

# Enumerate only the claimed labelled equality family, not 2^105 characters.
index={edge:j for j,edge in enumerate(combinations(range(1,n),2))}
full=(1<<len(index))-1


def root_gauged_mask(edges):
    raw=set(edges)
    root_negative={v for u,v in raw if u==0}
    out=0
    for edge,j in index.items():
        u,v=edge
        if (edge in raw)^((u in root_negative)!=(v in root_negative)):
            out|=1<<j
    return out


labels=set()
for a,b,c,d in combinations(range(n),4):
    for edges in ({(a,b),(c,d)},{(a,c),(b,d)},{(a,d),(b,c)}):
        labels.add(root_gauged_mask(edges))
check(len(labels)==3*comb(n,4),"labelled equality family has no switching collisions")
check(not(labels&{full^h for h in labels}),"global-negative equality family distinct")
print(f"literal_Hamilton_cycles={total_cycles} recurrence_and_sign_universe=3..100")
print("signed_sum_hostiles=W3(1):-1,W5(1):0,W6(3):-4")
print(f"n16_equality_family={len(labels)} even_global_signed_family={2*len(labels)}")
print(f"explicit_failure_gates={gates}")
print("RESULT=PASS")
