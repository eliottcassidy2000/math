"""Exact Boolean orbit-graph controls; no weighted/Boolean spectral identity.

Universe: all Eulerian graphs on n=3,...,6, every cycle prefix, and all
isomorphism classes. Independent edge-subset recognizers audit the primal
anchored-triangle/cyclic-order constructors through n=5.
"""
from fractions import Fraction as Q
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path
import sympy as s

spec=spec_from_file_location("inherited_primal",Path(__file__).with_name(
    "even_graph_triangle_quotient_spectrum_thm4078.py"))
base=module_from_spec(spec);spec.loader.exec_module(base)
gates=0


def check(test,label):
    global gates
    if not test:
        raise RuntimeError(label)
    gates+=1


def literal_sets(n,edges):
    euler=set();cycles={k:set() for k in range(3,n+1)}
    for mask in range(1<<len(edges)):
        adj=[set() for _ in range(n)]
        for j,(a,b) in enumerate(edges):
            if mask>>j&1:
                adj[a].add(b);adj[b].add(a)
        if any(len(row)%2 for row in adj):
            continue
        euler.add(mask)
        used={v for v in range(n) if adj[v]}
        if not used or any(len(adj[v])!=2 for v in used):
            continue
        seen={min(used)};todo=list(seen)
        while todo:
            for v in adj[todo.pop()]-seen:
                seen.add(v);todo.append(v)
        if seen==used:
            cycles[len(used)].add(mask)
    return euler,cycles


records=[]
for n in range(3,7):
    reps,lookup,sizes,_=base.quotient_data(n)
    edges,idx=base.edge_system(n)
    check(len(reps)=={3:2,4:3,5:7,6:16}[n],"complete primal orbit count")
    cycles={k:base.simple_cycles(n,k,idx) for k in range(3,n+1)}
    if n<=5:
        exact_euler,exact_cycles=literal_sets(n,edges)
        check(exact_euler==set(lookup),"independent complete Eulerian universe")
        for k in cycles:
            check(set(cycles[k])==exact_cycles[k],"independent connected 2-regular cycles")
    matrices={k:s.Matrix(base.weighted_operator(reps,lookup,cycles[k])) for k in cycles}
    for k in matrices:
        for ell in matrices:
            check(matrices[k]*matrices[ell]==matrices[ell]*matrices[k],"weighted commuting lift")
    q=len(reps);one=s.ones(q,1);orb=s.Matrix(sizes)
    psi=s.Matrix([n*(n-1)//2-2*mask.bit_count() for mask in reps])
    M=s.zeros(q);separate=s.zeros(q);N=0;E=0
    for k in range(3,n+1):
        M+=matrices[k];N+=len(cycles[k])
        E+=s.factorial(n-2)//s.factorial(n-k)
        separate+=s.Matrix(base.boolean_support(matrices[k].tolist()))
        B=s.Matrix(base.boolean_support(M.tolist()))
        degrees=B*one;LB=s.diag(*degrees)-B;LW=N*s.eye(q)-M
        check(M*one==N*one,"weighted regular row total")
        check(B==B.T,"Boolean undirected support")
        check(s.diag(*sizes)*M==M.T*s.diag(*sizes),"native orbit detailed balance")
        check(s.diag(*[N-M[i,i] for i in range(q)])-(M-s.diag(*[M[i,i] for i in range(q)]))==LW,
              "degree-adjusted loop deletion leaves Laplacian unchanged")
        eigen=LW.eigenvals()
        check(min(value for value in eigen if value>0)==2*E,"proved weighted gap control")
        rank=s.Matrix.hstack(one,psi,LB*psi).rank()
        check(rank==(2 if (n,k) in ((3,3),(4,4)) else 3),"affine Fourier transfer hostile")
        overlap=sum(int(bool(separate[i,j]>1)) for i in range(q) for j in range(i))
        check(overlap=={(5,5):6,(6,5):23,(6,6):44}.get((n,k),0),"cross-length support collapse")
        conductance=[int(sizes[i]*M[i,j]) for i in range(q) for j in range(i) if B[i,j]]
        cmin,cmax=min(conductance),max(conductance)
        check(cmax<=s.factorial(n)*N,"universal native conductance ceiling")
        for i in range(q):
            for j in range(q):
                if i!=j and M[i,j]:
                    check(Q(int(M[i,j]),N)/int(M[i,j])==Q(1,N),"reciprocal multiplicity thinning")
        thin=s.eye(q)-LB/N
        check(thin==thin.T and thin*one==one and min(thin)>=0,"native uniform reversible thinned walk")
        variance=s.diag(*sizes)-orb*orb.T/sum(sizes)
        uniform=s.eye(q)-one*one.T/q
        for f in [psi]+[s.eye(q)[:,i] for i in range(q)]:
            qw=(f.T*s.diag(*sizes)*LW*f)[0]
            qb=(f.T*LB*f)[0]
            vs=(f.T*variance*f)[0];vu=(f.T*uniform*f)[0]
            check(cmin*qb<=qw<=cmax*qb,"native Dirichlet comparison controls")
            check(min(sizes)*vu<=vs<=max(sizes)*vu,"centered variance comparison controls")
        check(degrees[lookup[0]]==k-2,"marked empty degree")
        polynomial=s.factor(LB.charpoly().as_expr())
        if (n,k)==(5,3):
            x=s.Symbol("lambda")
            check(s.expand(polynomial-x*(x-3)*(x-1)*(x*x-8*x+14)*(x*x-4*x+2))==0,
                  "exact triangle n5 Boolean spectrum")
        if (n,k)==(5,5):
            x=s.Symbol("lambda")
            check(s.expand(polynomial-x*(x-7)*(x-6)**2*(x-4)*(x*x-9*x+16))==0,
                  "exact envelope n5 Boolean spectrum")
        row=(n,k-1,q,int(min(degrees)),int(max(degrees)),rank,overlap,
             int(sum(M[i,i] for i in range(q))),cmin,cmax,Q(2*int(E),cmax),
             Q((k-2)*q,q-1),sha256(str(polynomial).encode()).hexdigest())
        records.append(row)
        print("ROW n,D,q,degree_min,max,affine_rank,overlaps,raw_loop_sum,cmin,cmax,gap_lower,empty_upper,poly_sha",row)
    if n==6:
        tri=lookup[cycles[3][0]];empty=lookup[0]
        B3=s.Matrix(base.boolean_support(matrices[3].tolist()))
        L3=s.diag(*(B3*one))-B3
        check(sum(B3[tri,j] for j in range(q))==4,"four complete triangle neighbor types")
        check(((L3*psi)[empty],(L3*psi)[tri])==(6,8),"all-order negative-eigenvalue hostile")

print("EXACT_BOOLEAN_GAPS n5_D2=2-sqrt(2) n5_D4=(9-sqrt(17))/2")
print("semantic_sha256",sha256(repr(records).encode()).hexdigest())
print("explicit_failure_gates",gates)
print("PASS exact Boolean consumer controls; no exact all-order Boolean gap claim")
