#!/usr/bin/env python3
"""Exact small controls for ternary interfaces, golden defects and ROOT charges."""
import argparse
import itertools as it
import json
from fractions import Fraction as Q
from pathlib import Path
import networkx as nx
import sympy as sp


def add(x,y):return (x[0]+y[0],x[1]+y[1])
def mul(x,y):return (x[0]*y[0]-x[1]*y[1],x[0]*y[1]+x[1]*y[0]-x[1]*y[1])
def norm(x):return x[0]**2-x[0]*x[1]+x[1]**2
def scale(x,n):return (n*x[0],n*x[1])
def phase(x):return [2*x[0]-x[1],-x[0]-x[1],-x[0]+2*x[1]]
Z=[(1,0),(0,1),(-1,-1)]


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,required=True);args=ap.parse_args()
    # Explicit, labelled extension of the incoming Seidel dictionary.
    pairs=list(it.combinations(range(5),2)); P=nx.Graph();P.add_nodes_from(pairs)
    P.add_edges_from((a,b) for a,b in it.combinations(pairs,2) if not set(a)&set(b))
    assert nx.is_isomorphic(P,nx.petersen_graph())
    v=(3,4); shore=set(P[v]); R=P.copy()
    for x in shore:
        for y in set(P)-shore:
            if R.has_edge(x,y):R.remove_edge(x,y)
            else:R.add_edge(x,y)
    assert R.degree(v)==0
    coord={}
    for i in range(3):
        coord[tuple(j for j in range(3) if j!=i)]=(i,i)
        coord[(i,3)]=((i+1)%3,(i-1)%3)
        coord[(i,4)]=((i-1)%3,(i+1)%3)
    for x,y in it.combinations(coord,2):
        a,b=coord[x],coord[y];rook=a[0]==b[0] or a[1]==b[1]
        circle=((a[0]-b[0])**2+(a[1]-b[1])**2)%3==1
        assert R.has_edge(x,y)==rook==circle
        aa=((a[0]+a[1])%3,(a[0]-a[1])%3);bb=((b[0]+b[1])%3,(b[0]-b[1])%3)
        assert rook != (aa[0]==bb[0] or aa[1]==bb[1])
    # The norm-one/rook/Paley coincidence is exceptional to q=3.
    field_degrees=[]
    for p in (3,5,7,11,19):
        squares={a*a%p for a in range(1,p)}
        eta=next(a for a in range(1,p) if a not in squares)
        norms=[(a*a-eta*b*b)%p for a in range(p) for b in range(p)]
        circle=norms.count(1);paley=sum(n in squares for n in norms);rook=2*(p-1)
        assert circle==p+1 and paley==(p*p-1)//2
        assert (circle==rook)==(p==3)
        field_degrees.append(dict(q=p,nonsquare=eta,norm_one=circle,norm_square=paley,rook=rook))
    # Golden reduction is F9; Eisenstein reduction is ramified, not F9.
    def f9mul(x,y):return ((x[0]*y[0]-x[1]*y[1])%3,(x[0]*y[1]+x[1]*y[0])%3)
    phi9=(2,2);v=(1,0);orbit=[]
    for _ in range(8):v=f9mul(v,phi9);orbit.append(v)
    assert len(set(orbit))==8 and orbit[-1]==(1,0)
    eps=(2,1) # omega-1 in the Eisenstein basis modulo 3
    assert tuple(t%3 for t in mul(eps,eps))==(0,0) and eps!=(0,0)
    for a,b in it.product(range(3),repeat=2):assert norm((a,b))%3==(a+b)**2%3
    # Independent audit of the newly incoming conic dictionary. Group closure
    # uses projective matrices, not its permutations on projective points.
    generators=[((1,-2,2),(2,-1,2),(2,-2,3)),((1,2,2),(2,1,2),(2,2,3)),((-1,2,2),(-2,1,2),(-2,2,3))]
    conic_checks=[]
    for p in (5,7):
        def canon(v):
            v=tuple(x%p for x in v);u=pow(next(x for x in v if x),-1,p)
            return tuple(x*u%p for x in v)
        pts=sorted({canon(v) for v in it.product(range(p),repeat=3) if any(v)})
        def bf(x,y):return (x[0]*y[0]+x[1]*y[1]-x[2]*y[2])%p
        conic=[x for x in pts if bf(x,x)==0]
        interior=[x for x in pts if bf(x,x)!=0 and not any(bf(x,y)==0 for y in conic)]
        Gi=nx.Graph();Gi.add_nodes_from(interior)
        Gi.add_edges_from((x,y) for x,y in it.combinations(interior,2) if bf(x,y)==0)
        target=nx.petersen_graph() if p==5 else nx.line_graph(nx.heawood_graph())
        matcher=nx.algorithms.isomorphism.GraphMatcher(Gi,target)
        assert matcher.is_isomorphic()
        gens=[canon(tuple(x for r in g for x in r)) for g in generators]
        def mm(a,b):return canon(tuple(sum(a[3*i+k]*b[3*k+j] for k in range(3)) for i in range(3) for j in range(3)))
        ident=(1,0,0,0,1,0,0,0,1);group={ident};todo=[ident]
        for a in todo:
            for b in gens:
                v=mm(a,b)
                if v not in group:group.add(v);todo.append(v)
        assert len(group)==(60 if p==5 else 336)
        aut=None
        if p==5:
            aut=sum(1 for _ in nx.algorithms.isomorphism.GraphMatcher(Gi,Gi).isomorphisms_iter())
            assert aut==120 and len(group)<aut
        conic_checks.append(dict(p=p,interior_order=len(Gi),isomorphic_to='Petersen' if p==5 else 'L(Heawood)',
                                 berggren_projective_matrix_order=len(group),automorphism_order=aut))
    # K5 Delta-Y is K3,3 plus the extra edge on the opposite two vertices.
    D=nx.complete_graph(5);D.remove_edges_from(it.combinations(range(3),2));D.add_node(5)
    D.add_edges_from((5,i) for i in range(3));D.remove_edge(3,4)
    assert nx.is_isomorphic(D,nx.complete_bipartite_graph(3,3))
    # All three-port assignments, all q=4..9: no floating point comparisons.
    star_checks=0;table=[]
    for q in range(4,10):
        C=Q(q-3); defect=Q(q*q-3*q+1,(q-3)**2)
        for ports in it.product(range(q),repeat=3):
            actual=sum(all(c!=p for p in ports) for c in range(q))
            ds=[int(ports[i]==ports[j]) for i,j in ((0,1),(1,2),(0,2))]
            prediction=C
            for d in ds:prediction*=1+Q(d,C)
            prediction-=defect*int(len(set(ports))==1)
            assert actual==prediction==q-len(set(ports));star_checks+=1
        table.append(dict(q=q,numerator=q*q-3*q+1,defect=str(defect)))
    q=sp.symbols('q'); phi=(1+sp.sqrt(5))/2
    assert sp.expand((q-phi**2)*(q-phi**-2)).simplify()==q*q-3*q+1
    assert sp.simplify(phi**5-1-phi**3*(5-phi**2))==0
    assert sp.simplify(phi**7-1-phi**3*(7-phi**-2))==0
    assert sp.simplify(phi**9-1-2*phi**5*(6-phi**2))==0
    # The seven-colour star defect is exactly the inherited c=-29/16
    # parameter, and it is the unique rational 3-cycle in this integer-q family.
    x,c=sp.symbols('x c');f1=x*x+c;f2=f1*f1+c;f3=f2*f2+c
    dyn=sp.div(f3-x,f1-x,x)[0]
    trace=x+f1+f2
    assert sp.rem(trace*trace+trace+c+2,dyn,x)==0
    trace_discriminants={d:-3*d*d+12*d+4 for d in range(1,5)}
    assert trace_discriminants=={1:13,2:16,3:13,4:4}
    cyc=[Q(-1,4),Q(-7,4),Q(5,4)]
    assert [a*a-Q(29,16) for a in cyc]==cyc[1:]+cyc[:1]
    assert 8*cyc[0]*cyc[1]*cyc[2]==Q(35,8)
    # c=-11/4: every rational periodic point would have denominator 2
    # and |x| < 5/2; all four possible odd numerators escape that interval.
    for a in [Q(-3,2),Q(-1,2),Q(1,2),Q(3,2)]:
        for _ in range(3):
            a=a*a-Q(11,4)
            if abs(a)>Q(5,2):break
        assert abs(a)>Q(5,2)
    # Paths retain positive boundary weight, so bare subdivisions never force
    # a new boundary relation at q>=3.
    B=sp.ones(6)-sp.eye(6)
    assert B**2==4*sp.ones(6)+sp.eye(6)
    # Primitive positive zero relations of the six Eisenstein units.
    U=[(1,0),(1,1),(0,1),(-1,0),(-1,-1),(0,-1)]
    zeros=[]
    for c in it.product(range(4),repeat=6):
        if not any(c):continue
        if (sum(c[i]*U[i][0] for i in range(6)),sum(c[i]*U[i][1] for i in range(6)))==(0,0):zeros.append(c)
    primitive=[c for c in zeros if not any(d!=c and all(x<=y for x,y in zip(d,c)) for d in zeros)]
    expected=[tuple(int(j in (i,i+3)) for j in range(6)) for i in range(3)]
    expected += [(1,0,1,0,1,0),(0,1,0,1,0,1)]
    assert set(primitive)==set(expected)
    # Source shells and the explicit norm upper bound, independent integer
    # implementation (all N<=200, sparse large N).
    shell_checks=0;source_rows=[]
    for N in list(range(1,201))+[1024,4096,7910,16384]:
        total=0;mass=0
        for s in range(1,4*N,2):
            if s%8==5:continue
            A=(0,0);c=[0,0,0];m=s
            while m<4*N:
                A=add(A,mul(Z[m%3],(1,2) if m<2*N else (0,1)))
                for r in range(3):c[r]+=([0,-3,3] if m<2*N else [-1,-1,2])[(m+r)%3]
                m=4*m+1
            assert phase(A)==c and norm(A) in (0,1,3)
            assert sum(x*x for x in c)==6*norm(A)
            total+=norm(A);mass+=max(abs(A[0]),abs(A[1]),abs(A[0]-A[1]));shell_checks+=1
        J=1;v=1
        while v<12*N+1:v*=4;J+=1
        assert total<=Q(17,7)*N+12*J+1
        assert mass<=Q(41,21)*N+9*J+1
        if N>=1024:source_rows.append(dict(N=N,norm=total,atomic_mass=mass,bound=str(Q(17,7)*N+12*J+1)))
    # Zero parity value does not erase integer-clock lag: exact quotient.
    z=sp.symbols('z')
    for L in range(1,64,2):
        quotient=sum((-z)**j for j in range(L))
        assert sp.expand((1+z)*quotient)==1+z**L
        assert quotient.subs(z,-1)==L
    # Exact hostile to all-direction Gram contraction, even at a large N.
    p=Path(__file__).parent/'frontier_interfaces_20261018j/root_phase.json'
    dat=json.loads(p.read_text())['rows'][-1];S=sp.Matrix(dat['source_gram']);G=sp.Matrix(dat['arrival_gram'])
    hostile=None
    for v in it.product(range(-3,4),repeat=3):
        if sum(v) or not any(v):continue
        vv=sp.Matrix(v);debt=(vv.T*(G-S)*vv)[0]
        if debt>0:hostile=dict(N=dat['N'],direction=v,excess=int(debt));break
    assert hostile
    out=dict(status='FINITE-EXACT',petersen_coordinates=[dict(edge=k,grid=v) for k,v in coord.items()],
             seidel_shore=sorted(shore),finite_field_degrees=field_degrees,incoming_conic_audit=conic_checks,ternary_ring_controls=dict(golden_phi_order=8,eisenstein_nonzero_nilpotent=eps),star_assignment_checks=star_checks,star_table=table,
             source_shell_checks=shell_checks,source_rows=source_rows,unit_zero_circuits=primitive,
             odd_lag_quotient_checks=32,
             star_cycle=dict(parameter='-(q^2-3q+1)/(q-3)^2',integer_q_at_least_4_with_rational_3cycle=[7],
                             cycle=[str(a) for a in cyc],multiplier='35/8',trace_discriminants=trace_discriminants),
             gram_hostile=hostile,scope='Small graph maps and algebraic identities; no transfer of convergence or six-colour forcing.')
    args.out.parent.mkdir(parents=True,exist_ok=True);args.out.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:v for k,v in out.items() if k!='petersen_coordinates'},indent=2))


if __name__=='__main__':main()
