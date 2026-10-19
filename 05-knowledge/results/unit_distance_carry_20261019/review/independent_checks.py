#!/usr/bin/env python3
"""Independent integer checks; no imports from the reviewed implementation."""
from collections import deque
from itertools import product, combinations
from pathlib import Path
import json

HERE=Path(__file__).resolve().parent

def norm(v,m,T,N):
    x,y=v
    return (x*x+T*x*y+N*y*y)%m

def mul(u,v,m,T,N):
    x,y=u; a,b=v
    return ((x*a-N*y*b)%m,(x*b+y*a+T*y*b)%m)

def quotient_by_conjugate(w,m,T,N):
    # w / conjugate(w) = w^2 / Norm(w).
    den=pow(norm(w,m,T,N),-1,m)
    z=mul(w,w,m,T,N)
    return tuple(den*a%m for a in z)

def units(p,n,T,N):
    m=p**n
    return [v for v in product(range(m),repeat=2) if norm(v,m,T,N)==1]

def local_checks(p,n,T,N):
    m=p**n; U=set(units(p,n,T,N))
    quotients={quotient_by_conjugate(w,m,T,N) for w in product(range(m),repeat=2) if norm(w,p,T,N)}
    assert quotients==U
    assert len(U)==(p+1)*p**(n-1)
    return {'p':p,'depth':n,'unit_count':len(U),'hilbert90_quotients_equal_norm_one':True}

def mixed_depth_image(p,T,N):
    m1,m2=p,p*p
    steps=[a+b for a,b in product(units(p,1,T,N),units(p,2,T,N))]
    zero=(0,0,0,0); seen={zero}; pending=deque([zero])
    while pending:
        v=pending.popleft()
        for s in steps:
            w=tuple((v[i]+s[i])%(m1 if i<2 else m2) for i in range(4))
            if w not in seen:seen.add(w);pending.append(w)
    assert len(seen)==p**6
    return {'p':p,'depths':[1,2],'vertices_generated':len(seen),'full_vertex_image':True}

def defect_witness(p):
    data=json.loads((HERE/'section-witnesses.json').read_text())[str(p)]
    V=list(product(range(p),repeat=2)); H=data['h']
    D=next(a for a in range(1,p) if pow(a,(p-1)//2,p)==p-1)
    L=[(a+p*h,b+p*k) for (a,b),(h,k) in zip(V,H)]
    edge_count=bad=0
    failed=set()
    for i,j in combinations(range(p*p),2):
        dx,dy=V[j][0]-V[i][0],V[j][1]-V[i][1]
        if (dx*dx-D*dy*dy)%p!=1:continue
        edge_count+=1
        dx,dy=L[j][0]-L[i][0],L[j][1]-L[i][1]
        if (dx*dx-D*dy*dy)%(p*p)!=1:
            bad+=1;failed.add(frozenset((V[i],V[j])))
    directions=[u for u in product(range(p),repeat=2) if norm(u,p,0,-D)==1 and u<tuple(-x%p for x in u)]
    cycles=[]
    for u in directions:
        remaining=set(V)
        while remaining:
            v=min(remaining)
            line=[tuple((v[k]+t*u[k])%p for k in (0,1)) for t in range(p)]
            remaining.difference_update(line)
            cycle={frozenset((line[t],line[(t+1)%p])) for t in range(p)}
            cycles.append(cycle)
    all_edges=set().union(*cycles)
    assert sum(map(len,cycles))==len(all_edges)==edge_count
    assert all(len(c&failed)==1 for c in cycles)
    assert bad==len(cycles)==p*(p+1)//2
    return {'p':p,'edges':edge_count,'straight_cycles':len(cycles),'failed_edges':bad,'one_failure_per_cycle':True,'exact_minimum_defect':bad}

# F_9 = F_3[t]/(t^2+1), and its unramified lift R=(Z/9)[t]/(t^2+1).
def q9_carry():
    p=3;m=9
    def a(x,y):return tuple((x[i]+y[i])%m for i in (0,1))
    def neg(x):return tuple(-z%m for z in x)
    def mult(x,y):return mul(x,y,m,0,1)
    def power(x,e):
        out=(1,0)
        while e:
            if e&1:out=mult(out,x)
            x=mult(x,x);e//=2
        return out
    residues=list(product(range(3),repeat=2))
    D=next(z for z in residues if tuple(t%3 for t in power(z,4))==(2,0))
    # S=R[y]/(y^2-D) is the unramified quadratic extension of R.
    def nrm(x,y):return a(mult(x,x),neg(mult(D,mult(y,y))))
    xs=[(1+3*i,3*j) for i,j in residues]
    ys=[(3*i,3*j) for i,j in residues]
    steps=[(x,y) for x,y in product(xs,ys) if nrm(x,y)==(1,0)]
    reachable={((0,0),(0,0))}
    for _ in range(3):reachable={(a(x,u),a(y,v)) for (x,y),(u,v) in product(reachable,steps)}
    assert len(steps)==9 and len(reachable)==9
    assert {x for x,y in reachable}=={(3,0)}
    assert ((0,0),(0,0)) not in reachable
    return {'p':3,'q':9,'base_nonsquare':D,'steps_reducing_to_one':len(steps),'three_step_endpoints':len(reachable),'zero_endpoint_possible':False}

def higher_depth_cycles(p,r,D):
    low=p**r; high=p**(r+1)
    steps=[(1+low*a,low*b) for a,b in product(range(p),repeat=2)
           if ((1+low*a)**2-D*(low*b)**2)%high==1]
    endpoints={(0,0)}
    for _ in range(low):
        endpoints={((a+c)%high,(b+d)%high) for (a,b),(c,d) in product(endpoints,steps)}
    assert {a for a,b in endpoints}=={low}
    assert (0,0) not in endpoints
    if p==3:
        U=units(p,r,0,-D)
        assert all((a+b)%3 for a,b in U)
    return {'p':p,'lower_depth':r,'cycle_length':low,'upper_depth':r+1,
            'permitted_steps':len(steps),'zero_endpoint_possible':False,
            'depth_r_q3_three_coloring_checked': True if p==3 else None}

def rigidity_barrier_control():
    cert=json.loads((HERE/'rigidity-witness.json').read_text())
    p,D=17,3
    vertices=[tuple(v) for v in cert['vertices']]
    assert len(vertices)==len(set(vertices))==7
    edges=[(i,j) for i,j in combinations(range(7),2)
           if norm((vertices[j][0]-vertices[i][0],vertices[j][1]-vertices[i][1]),p,0,-D)==1]
    assert len(edges)==11
    rows=[]
    for i,j in edges:
        dx=vertices[j][0]-vertices[i][0];dy=vertices[j][1]-vertices[i][1]
        row=[0]*14
        row[2*i]=-2*dx;row[2*j]=2*dx
        row[2*i+1]=2*D*dy;row[2*j+1]=-2*D*dy
        rows.append([v%p for v in row])
    # Independent column-space elimination selects an invertible 11-column minor.
    basis={};columns=[]
    for j in range(14):
        v=[row[j] for row in rows]
        for pivot,b in sorted(basis.items()):
            a=v[pivot]
            v=[(x-a*y)%p for x,y in zip(v,b)]
        pivot=next((i for i,x in enumerate(v) if x),None)
        if pivot is not None:
            scale=pow(v[pivot],-1,p)
            basis[pivot]=[scale*x%p for x in v];columns.append(j)
    assert len(basis)==11
    minor=[[row[j] for j in columns] for row in rows]
    det=1
    for j in range(11):
        pivot=next(i for i in range(j,11) if minor[i][j]%p)
        if pivot!=j:minor[j],minor[pivot]=minor[pivot],minor[j];det=-det
        a=minor[j][j]%p;det=det*a%p
        for i in range(j+1,11):
            scale=minor[i][j]*pow(a,-1,p)%p
            minor[i]=[(x-scale*y)%p for x,y in zip(minor[i],minor[j])]
    assert det%p
    kernel=[[1,0]*7,[0,1]*7,[z for x,y in vertices for z in (D*y,x)]]
    assert all(sum(a*b for a,b in zip(row,v))%p==0 for row in rows for v in kernel)
    coloring=cert['four_coloring']
    assert all(coloring[i]!=coloring[j] for i,j in edges)
    assert not any(all(c[i]!=c[j] for i,j in edges) for c in product(range(3),repeat=7))
    max_min_degree=0
    for mask in range(1,1<<7):
        selected={i for i in range(7) if mask>>i&1}
        chosen=[(i,j) for i,j in edges if i in selected and j in selected]
        degrees=[sum(i in edge for edge in chosen) for i in selected]
        max_min_degree=max(max_min_degree,min(degrees))
        if len(selected)>=2:assert len(chosen)<=2*len(selected)-3
    assert max_min_degree==3
    for lift in cert['exact_lifts']:
        m=lift['modulus'];L=lift['vertices']
        assert [(x%p,y%p) for x,y in L]==vertices
        assert all(norm((L[j][0]-L[i][0],L[j][1]-L[i][1]),m,0,-D)==1 for i,j in edges)
    return {'p':p,'vertices':7,'edges':len(edges),'rigidity_rank':11,
            'independent_minor_columns':columns,'minor_determinant_mod_p':det%p,
            'chromatic_number':4,'degeneracy':max_min_degree,
            'exact_lift_moduli':[z['modulus'] for z in cert['exact_lifts']]}

out={
 'local_norm_one_checks':[local_checks(*args) for args in [(2,1,-1,1),(2,2,-1,1),(2,3,-1,1),(3,1,0,-2),(3,2,0,-2),(5,2,0,-2),(7,2,0,-3)]],
 'mixed_depth_additive_generation':[mixed_depth_image(2,-1,1),mixed_depth_image(3,0,-2)],
 'prime_power_residue_carry':q9_carry(),
 'independently_checked_defect_witnesses':[defect_witness(3),defect_witness(5)],
 'rigidity_barrier_sharpness':rigidity_barrier_control(),
 'higher_depth_cycle_checks':[higher_depth_cycles(3,2,2),higher_depth_cycles(3,3,2),higher_depth_cycles(5,2,2)],
 'field_application':{'p':29,'cubic_roots':[a for a in range(29) if (a**3+2*a*a-4*a-4)%29==0],'remaining_quadratic_discriminant':10,'discriminant_is_nonsquare':pow(10,14,29)==28}
}
assert out['field_application']['cubic_roots']==[3]
assert out['field_application']['discriminant_is_nonsquare']
(HERE/'independent-results.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
