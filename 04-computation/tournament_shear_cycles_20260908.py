"""Exact tournament shear-production and weighted-substitution audit.

Declared universes: all labeled unit tournaments through order six;
all order-four tournaments with edge magnitudes in {1,2}; all order-three
tournaments with magnitudes in {1,2,3}; 4160 weighted substitutions from
the four explicit blocks; and named boundary/strong-core hostiles.
Standard library only. Checks remain active with python -O.
"""
from collections import Counter
from fractions import Fraction as F
from itertools import combinations,product
import hashlib
import json

GATES=0


def check(ok,label):
    global GATES
    GATES+=1
    if not ok:
        raise RuntimeError(label)


def adjacency(n,mask,weights=None):
    pairs=list(combinations(range(n),2))
    if weights is None:
        weights=[1]*len(pairs)
    M=[[0]*n for _ in range(n)]
    for s,((i,j),w) in enumerate(zip(pairs,weights)):
        if (mask>>s)&1:
            M[j][i]=w
        else:
            M[i][j]=w
    return M


def transpose(M):
    return list(map(list,zip(*M)))


def matmul(A,B):
    return [[sum(x*y for x,y in zip(row,col)) for col in zip(*B)] for row in A]


def tr(A):
    return sum(A[i][i] for i in range(len(A)))


def cyclic(M,i,j,k):
    return bool(M[i][j] and M[j][k] and M[k][i]) or bool(M[i][k] and M[k][j] and M[j][i])


def triangle_value(M):
    total=0
    for i,j,k in combinations(range(len(M)),3):
        w=(M[i][j]+M[j][i])*(M[i][k]+M[k][i])*(M[j][k]+M[k][j])
        total+=(3 if cyclic(M,i,j,k) else -1)*w
    return total


def trace_value(M):
    n=len(M)
    B=[[M[i][j]+M[j][i] for j in range(n)] for i in range(n)]
    C=[[M[i][j]-M[j][i] for j in range(n)] for i in range(n)]
    trace_formula=F(tr(matmul(B,matmul(C,C))),2)
    M2=matmul(M,M)
    check(trace_formula == tr(matmul(M2,M))-tr(matmul(M2,transpose(M))),"expanded trace identity")
    BC=matmul(B,C)
    CB=matmul(C,B)
    derivative=-F(sum(C[i][j]*(BC[i][j]+CB[i][j]) for i in range(n) for j in range(n)),4)
    check(derivative == trace_formula,"Euler skew-norm derivative")
    return trace_formula


def strong(M):
    n=len(M)
    reach=[(1<<i)|sum((1<<j) for j in range(n) if M[i][j]) for i in range(n)]
    for k in range(n):
        for i in range(n):
            if (reach[i]>>k)&1:
                reach[i]|=reach[k]
    return all(row==(1<<n)-1 for row in reach)


def regular_after_extreme_deletion(M):
    n=len(M)
    for v in range(n):
        degree=sum(bool(t) for t in M[v])
        if degree not in (0,n-1):
            continue
        vertices=[i for i in range(n) if i!=v]
        if all(sum(bool(M[i][j]) for j in vertices) == F(n-2,2) for i in vertices):
            return True
    return False


def edge_mass(M):
    return sum(map(sum,M))


def substitute(Q,blocks):
    sizes=[len(B) for B in blocks]
    origins=[sum(sizes[:i]) for i in range(len(sizes))]
    n=sum(sizes)
    out=[[0]*n for _ in range(n)]
    for b,B in enumerate(blocks):
        for i in range(sizes[b]):
            for j in range(sizes[b]):
                out[origins[b]+i][origins[b]+j]=B[i][j]
    for b,c in combinations(range(len(blocks)),2):
        for i in range(sizes[b]):
            for j in range(sizes[c]):
                out[origins[b]+i][origins[c]+j]=Q[b][c]
                out[origins[c]+j][origins[b]+i]=Q[c][b]
    return out


def substitution_value(Q,blocks):
    sizes=[len(B) for B in blocks]
    masses=[edge_mass(B) for B in blocks]
    out=sum(triangle_value(B) for B in blocks)
    for i,j in combinations(range(len(blocks)),2):
        lam=Q[i][j]+Q[j][i]
        out-=lam*lam*(sizes[j]*masses[i]+sizes[i]*masses[j])
    for i,j,k in combinations(range(len(blocks)),3):
        weight=(Q[i][j]+Q[j][i])*(Q[i][k]+Q[k][i])*(Q[j][k]+Q[k][j])
        out+=(3 if cyclic(Q,i,j,k) else -1)*weight*sizes[i]*sizes[j]*sizes[k]
    return out


def attach_source(B,b):
    out=[[0]+list(b)]
    out.extend([[0]+list(row) for row in B])
    return out


def boundary_tax(B,b):
    return sum((B[i][j]+B[j][i])*b[i]*b[j] for i,j in combinations(range(len(B)),2))


def main():
    result={}
    unit_count=0
    census={}
    for n in range(1,7):
        hist=Counter()
        for mask in range(1<<(n*(n-1)//2)):
            M=adjacency(n,mask)
            f=triangle_value(M)
            c3=sum(cyclic(M,*tri) for tri in combinations(range(n),3))
            check(f == 4*c3-n*(n-1)*(n-2)//6,"unit cycle identity")
            degrees=list(map(sum,M))
            check(f == F(n*(n-1),2)-2*sum((F(d)-F(n-1,2))**2 for d in degrees),"degree-variance identity")
            is_strong=strong(M)
            hist[("strong" if is_strong else "reducible", "positive" if f>0 else "negative" if f<0 else "zero")]+=1
            if n>=2 and not is_strong:
                check(f<=0,"unit reducibility sign")
                check((f==0) == regular_after_extreme_deletion(M),"unit reducibility zero classification")
            if n<=5:
                check(f == trace_value(M),"unit matrix-triple identity")
            unit_count+=1
        census[n]={"/".join(k):v for k,v in sorted(hist.items())}
    result["unit_tournaments"]=unit_count
    result["unit_strong_sign_census"]=census
    weighted_count=0
    for n,weights in [(3,(1,2,3)),(4,(1,2))]:
        edges=n*(n-1)//2
        for mask in range(1<<edges):
            for w in product(weights,repeat=edges):
                M=adjacency(n,mask,w)
                f=triangle_value(M)
                check(f==trace_value(M),"weighted matrix-triple identity")
                if n==3:
                    omega=[M[2][1]-M[1][2],M[0][2]-M[2][0],M[1][0]-M[0][1]]
                    check(f == sum(omega[i]*M[i][j]*omega[j] for i in range(3) for j in range(3)),"actual 3D stretching")
                weighted_count+=1
    result["weighted_matrix_instances"]=weighted_count
    C3=[[0,1,0],[0,0,1],[1,0,0]]
    WC3=[[0,2,0],[0,0,1],[1,0,0]]
    bank=[[[0]],[[0,1],[0,0]],C3,WC3]
    substitutions=0
    for q in (2,3):
        edges=q*(q-1)//2
        for mask in range(1<<edges):
            for magnitudes in product((1,2),repeat=edges):
                Q=adjacency(q,mask,magnitudes)
                for indices in product(range(len(bank)),repeat=q):
                    blocks=[bank[i] for i in indices]
                    A=substitute(Q,blocks)
                    check(triangle_value(A) == substitution_value(Q,blocks),"weighted substitution compiler")
                    predicted_mass=sum(edge_mass(B) for B in blocks)+sum((Q[i][j]+Q[j][i])*len(blocks[i])*len(blocks[j]) for i,j in combinations(range(q),2))
                    check(edge_mass(A) == predicted_mass,"substitution edge-mass closure")
                    substitutions+=1
    result["weighted_substitutions"]=substitutions
    # Exact directed-cut upper bounds, including parity, on every pair of
    # unit tournaments of total order at most seven from orders at most four.
    joins=0
    for a in range(1,5):
        for b in range(1,5):
            if a+b>7:
                continue
            for ma in range(1<<(a*(a-1)//2)):
                A=adjacency(a,ma)
                for mb in range(1<<(b*(b-1)//2)):
                    B=adjacency(b,mb)
                    T=substitute([[0,1],[0,0]],[A,B])
                    f=triangle_value(T)
                    exact=triangle_value(A)+triangle_value(B)-F(a*b*(a+b-2),2)
                    upper=-F((a+b)*(a-1)*(b-1)+(a if a%2==0 else 0)+(b if b%2==0 else 0),2)
                    check(f==exact and f<=upper,"directed-cut exact defect and bound")
                    joins+=1
    result["directed_cut_pairs"]=joins
    # Strong connectivity alone does not imply production: almost ordered
    # order six with the long arc reversed has exactly four cyclic triples.
    near=adjacency(6,0)
    near[0][5]=0
    near[5][0]=1
    check(strong(near) and triangle_value(near)==-4,"strong negative hostile")
    # Uniform exterior magnitude dominance does not repair weighted failure.
    pair=[[0,F(1,2)],[0,0]]
    core=substitute(C3,[pair,pair,pair])
    joined=attach_source(core,[1]*6)
    check(triangle_value(core)==18 and edge_mass(core)==F(27,2),"weighted six-core data")
    check(not strong(joined) and triangle_value(joined)==F(9,2),"weighted reducible positive hostile")
    # Same (n,e,F), same exterior labels, opposite total signs.
    B1=WC3
    B2=[[0,1,0],[0,0,1],[2,0,0]]
    b=[F(2,3),F(4,3),F(2)]
    check((len(B1),edge_mass(B1),triangle_value(B1))==(len(B2),edge_mass(B2),triangle_value(B2))==(3,4,6),"same scalar block state")
    taxes=[boundary_tax(B,b) for B in (B1,B2)]
    full=[triangle_value(attach_source(B,b)) for B in (B1,B2)]
    check(taxes == [F(52,9),F(56,9)] and full == [F(2,9),F(-2,9)],"boundary quadratic sign flip")
    for B,t,f in zip((B1,B2),taxes,full):
        check(triangle_value(B)-t==f,"exact singleton boundary law")
    result["hostiles"]={"strong_negative_unit_F":"-4","reducible_positive_weights_half_to_one_F":"9/2",
                         "same_state_boundary_taxes":list(map(str,taxes)),"same_state_full_F":list(map(str,full))}
    result["gates"]=GATES
    body=json.dumps(result,sort_keys=True,indent=2)
    print("TOURNAMENT SHEAR CYCLES: EXACT TRACE, SUBSTITUTION, AND BOUNDARY AUDIT")
    print(body)
    print("semantic_sha256="+hashlib.sha256(body.encode()).hexdigest())
    print("PASS: explicit checks survive -O; instantaneous algebra, not an Euler blowup claim")


if __name__ == "__main__":
    main()
