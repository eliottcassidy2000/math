"""Exact rational Fourier audit of periodic gradient minors and shear mixtures.

Universe: four explicitly declared three-dimensional trigonometric fields;
all 2x2 minors by Laurent convolution; independent 4^3 rational grid values;
the matched-energy three-shear mixture; eight exact sharpness parameters.
No numerical sampling, optimizer, repository imports, or erased assertions.
"""
from fractions import Fraction as F
from itertools import product
import hashlib
import json

ZERO = (0,0,0)
GATES = 0


def check(ok, label):
    global GATES
    GATES += 1
    if not ok:
        raise RuntimeError(label)


def addc(a,b):
    return (a[0]+b[0],a[1]+b[1])


def mulc(a,b):
    return (a[0]*b[0]-a[1]*b[1],a[0]*b[1]+a[1]*b[0])


def scale(P,t):
    return {k:mulc(v,(F(t),F(0))) for k,v in P.items() if t and v != (0,0)}


def add(*polys):
    out={}
    for P in polys:
        for k,v in P.items():
            out[k]=addc(out.get(k,(F(0),F(0))),v)
    return {k:v for k,v in out.items() if v != (0,0)}


def mul(P,Q):
    out={}
    for k,v in P.items():
        for ell,w in Q.items():
            key=tuple(x+y for x,y in zip(k,ell))
            out[key]=addc(out.get(key,(F(0),F(0))),mulc(v,w))
    return {k:v for k,v in out.items() if v != (0,0)}


def diff(P,j):
    return {k:mulc(v,(F(0),F(k[j]))) for k,v in P.items() if k[j]}


def constant(t):
    return {ZERO:(F(t),F(0))} if t else {}


def avg(P):
    re,im=P.get(ZERO,(F(0),F(0)))
    check(im == 0,"real average")
    return re


def velocity(waves):
    v=[{} for _ in range(3)]
    for k,a in waves:
        check(sum(x*y for x,y in zip(k,a)) == 0,"divergence-free amplitude")
        minus=tuple(-x for x in k)
        for i in range(3):
            # sin(k.x)=(exp(i k.x)-exp(-i k.x))/(2i).
            sine={k:(F(0),-F(a[i],2)),minus:(F(0),F(a[i],2))}
            v[i]=add(v[i],sine)
    return v


def minor(A,rows,cols,algebraic=False):
    i,j=rows
    k,l=cols
    if algebraic:
        return add(mul(A[i][k],A[j][l]),scale(mul(A[i][l],A[j][k]),-1))
    return A[i][k]*A[j][l]-A[i][l]*A[j][k]


def cof(A,algebraic=False):
    out=[]
    for i in range(3):
        row=[]
        for j in range(3):
            value=minor(A,[s for s in range(3) if s!=i],
                        [s for s in range(3) if s!=j],algebraic)
            row.append(scale(value,(-1)**(i+j)) if algebraic else (-1)**(i+j)*value)
        out.append(row)
    return out


def norm2(A):
    return sum(v*v for row in A for v in row)


def cos_quarter(q):
    return (1,0,-1,0)[q%4]


def direct_matrix(M,waves,point):
    # Independent differentiation and exact rational trigonometric values.
    return [[M[i][j]+sum(a[i]*k[j]*cos_quarter(sum(k[s]*point[s] for s in range(3)))
                        for k,a in waves)
             for j in range(3)] for i in range(3)]


def audit_field(M,waves,name):
    v=velocity(waves)
    check(add(*(diff(v[i],i) for i in range(3))) == {},"Fourier divergence")
    A=[[add(constant(M[i][j]),diff(v[i],j)) for j in range(3)] for i in range(3)]
    check(all(avg(A[i][j]) == M[i][j] for i in range(3) for j in range(3)),"gradient mean")
    check(add(*(A[i][i] for i in range(3))) == {},"trace-free field")
    C=cof(A,True)
    CM=cof(M)
    check([[avg(P) for P in row] for row in C] == CM,"Laurent cofactor mean")
    e=avg(add(*(mul(P,P) for row in A for P in row)))
    check(e == norm2(M)+sum(F(sum(x*x for x in a)*sum(x*x for x in k),2)
                            for k,a in waves),"independent Fourier orthogonality energy")
    # Every coordinate frequency in a minor or squared entry has magnitude
    # at most 2. Hence the 4-grid averages have no nonconstant alias.
    for row in C:
        for P in row:
            check(all(max(map(abs,k)) <= 2 for k in P),"cubature frequency bound")
    grid_cof=[[F(0) for _ in range(3)] for _ in range(3)]
    grid_energy=F(0)
    for point in product(range(4),repeat=3):
        B=direct_matrix(M,waves,point)
        CB=cof(B)
        check(sum(B[i][i] for i in range(3)) == 0,"point trace")
        grid_energy+=norm2(B)/F(64)
        for i in range(3):
            for j in range(3):
                grid_cof[i][j]+=CB[i][j]/F(64)
    check(grid_cof == CM,"independent exact-grid minors")
    check(grid_energy == e,"independent exact-grid energy")
    return {"name":name,"energy":str(e),"cofactor_mean":[list(map(str,row)) for row in CM],
            "cofactor_mean_norm2":str(norm2(CM)),
            "certified_mean_shear_distance2_lower":str(norm2(CM)/e) if e else "0"}, A,C


def elementary(i,j,t=1):
    return [[F(t) if (s,u)==(i,j) else F(0) for u in range(3)] for s in range(3)]


def main():
    M=[[F(0),F(1),F(0)],[F(0),F(0),F(1)],[F(1),F(0),F(0)]]
    waves=[((1,1,0),(1,-1,0)),((0,1,1),(0,1,-1)),((1,0,1),(1,0,-1))]
    result={"fields":[]}
    info,A,C=audit_field(M,waves,"matched_energy_cyclic_mean")
    result["fields"].append(info)
    # The diagonal minor has two nonzero quadratic averages which cancel.
    aa=avg(mul(A[0][0],A[1][1]))
    bb=avg(mul(A[0][1],A[1][0]))
    check(aa == bb == F(-1,2),"nontrivial quadratic cancellation")
    result["quadratic_minor_cancellation"]=[str(aa),str(bb),str(aa-bb)]
    # A two-grid is not a valid energy cubature: squared Fourier frequencies alias.
    energy2=sum(norm2(direct_matrix(M,waves,point)) for point in product((0,2),repeat=3))/F(8)
    check(energy2 == 15 and energy2 != F(info["energy"]),"coarse-grid alias hostile")
    result["coarse_two_grid_energy_hostile"] = str(energy2)
    Z=[[F(0) for _ in range(3)] for _ in range(3)]
    for mean,wavebank,name in [(Z,waves,"zero_mean_nonconstant"),
                               (M,[],"constant_full_rank"),
                               (elementary(0,1),[((0,1,0),(1,0,0))],"compatible_pure_shear_wave")]:
        data,_,C1=audit_field(mean,wavebank,name)
        result["fields"].append(data)
        if name == "compatible_pure_shear_wave":
            check(all(P == {} for row in C1 for P in row),"pointwise rank-one positive control")
    mixture=[elementary(0,1,3),elementary(1,2,3),elementary(2,0,3)]
    mixed_mean=[[sum(T[i][j] for T in mixture)/F(3) for j in range(3)] for i in range(3)]
    mixed_cof=[[sum(cof(T)[i][j] for T in mixture)/F(3) for j in range(3)] for i in range(3)]
    mixed_energy=sum(norm2(T) for T in mixture)/F(3)
    check(mixed_mean == M,"mixture matching mean")
    check(mixed_energy == 9,"mixture matching energy")
    check(mixed_cof == Z and norm2(cof(M)) == 3,"mixture violates gradient minors")
    check(all(sum(T[i][i] for i in range(3)) == 0 and cof(T)==Z for T in mixture),
          "mixture pure-shear support")
    result["mixture"]={"energy":str(mixed_energy),"mean_shear_distance2":"0",
                        "mean_cofactor_norm2":str(norm2(mixed_cof)),
                        "cofactor_of_mean_norm2":str(norm2(cof(M)))}
    sharp=[]
    for n in range(1,9):
        eps=F(1,n)
        S=[[F(0),F(1),F(0)],[F(0),F(0),eps],[F(0),F(0),F(0)]]
        gram=[[sum(S[k][i]*S[k][j] for k in range(3)) for j in range(3)] for i in range(3)]
        check(gram == [[0,0,0],[0,1,0],[0,0,eps*eps]],"sharpness singular square values")
        check(norm2(cof(S)) == eps*eps,"sharpness cofactor")
        check(norm2(S) == 1+eps*eps,"sharpness energy")
        ratio=norm2(cof(S))/((1+eps*eps)*(eps*eps))
        check(ratio == 1/(1+eps*eps),"sharpness quotient")
        sharp.append([str(eps),str(ratio)])
    result["sharpness_epsilon_and_ratio"] = sharp
    result["gates"] = GATES
    text=json.dumps(result,sort_keys=True,indent=2)
    print("PERIODIC GRADIENT MINORS: EXACT FOURIER AND RATIONAL CUBATURE PASS")
    print(text)
    print("semantic_sha256="+hashlib.sha256(text.encode()).hexdigest())
    print("PASS: all explicit gates remain active under -O; no PDE realization claimed")


if __name__ == "__main__":
    main()
