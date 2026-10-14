"""Recovered independent exact checks, 2026-10-14.
Run: python3 liouville_planar_obstructions_20261014.py > result.json
No network or repository dependency. Graph literals recovered from retained tool outputs.
These are new recovery bytes, not a claim of identity with the lost original files.
"""
import json
from functools import lru_cache
import sympy as s

SOURCE_COMMIT='bd2d0d6660b3889179e226172a9572351bd122bc'
GRAPHS={
'EI17': [[11,12,15],[10,12,16],[11,13,15],[10,14,16],[9,13,15],[8,14,16],[7,9,14,15],[6,8,13,16],[5,7,10,12],[4,6,11,12],[1,3,8,11],[0,2,9,10],[0,1,8,9],[2,4,7,14],[3,5,6,13],[0,2,4,6],[1,3,5,7]],
'EI19': [[4,8,15],[5,8,15],[3,6,18],[2,16,17],[0,5,14],[1,4,10],[2,7,13],[6,11,17],[0,1,10],[11,12,16],[5,8,16,18],[7,9,14,18],[9,14,15,17],[6,14,16,17],[4,11,12,13],[0,1,12,16,18],[3,9,10,13,15],[3,7,12,13,18],[2,10,11,15,17]],
'FISH': [[10,20,22],[11,21,22],[8,19,22],[9,18,22],[13,17,22],[12,16,22],[15,17,22],[14,16,22],[2,10,20],[3,11,21],[0,8,19],[1,9,18],[5,14,20],[4,15,21],[7,12,19],[6,13,18],[5,7,20],[4,6,21],[3,11,15,19,20],[2,10,14,18,21],[0,8,12,16,18,21],[1,9,13,17,19,20],[0,1,2,3,4,5,6,7]],
'EI21': [[4,8,19],[5,9,20],[7,9,20],[6,8,19],[0,6,14],[1,7,13],[3,4,11],[2,5,12],[0,3,11],[1,2,12],[13,14,18],[6,8,17,18],[7,9,16,18],[5,10,15,17],[4,10,15,16],[13,14,19,20],[12,14,17,20],[11,13,16,19],[10,11,12,19,20],[0,3,15,17,18],[1,2,15,16,18]]}

def check(ok,label):
    if not ok: raise RuntimeError(label)

def peel(adj):
    remaining=set(range(len(adj))); order=[]; width=0
    while remaining:
        v=min(remaining,key=lambda v:(len(set(adj[v]) & remaining),v))
        degree=len(set(adj[v]) & remaining);width=max(width,degree)
        order.append([v,degree]);remaining.remove(v)
    return width,order

def colour(adj,k):
    colours={}; calls=0
    def rec():
        nonlocal calls
        calls+=1
        if len(colours)==len(adj): return dict(colours)
        remaining=[v for v in range(len(adj)) if v not in colours]
        v=max(remaining,key=lambda v:(len({colours[w] for w in adj[v] if w in colours}),len(adj[v]),-v))
        used={colours[w] for w in adj[v] if w in colours}
        # Any colouring can be relabelled by first-use order, so this loses no solutions.
        cap=min(k,max(colours.values(),default=-1)+2)
        for c in range(cap):
            if c in used: continue
            colours[v]=c; answer=rec()
            if answer is not None:return answer
            del colours[v]
        return None
    answer=rec()
    if answer is not None:check(all(answer[v]!=answer[w] for v in range(len(adj)) for w in adj[v]),'proper colouring')
    return answer,calls

result={'provenance':{'type':'reconstructed_from_retained_history','original_source_commit':SOURCE_COMMIT,'original_byte_identity_claimed':False,'sympy':s.__version__},'graphs':{}}
expected_edges={'EI17':31,'EI19':35,'EI21':38,'FISH':42}
for name,adj in GRAPHS.items():
    check(all(v not in nb and len(nb)==len(set(nb)) for v,nb in enumerate(adj)),'simple graph')
    check(all(v in adj[w] for v in range(len(adj)) for w in adj[v]),'symmetric graph')
    check(sum(map(len,adj))//2==expected_edges[name],'edge count')
    check(not any(set(adj[v]) & set(adj[w]) for v in range(len(adj)) for w in adj[v]),'triangle free')
    width,order=peel(adj);c3,n3=colour(adj,3);c4,n4=colour(adj,4)
    check(c3 is None and c4 is not None,'four chromatic '+name)
    check(width==3,'three degenerate '+name)
    remaining=set(range(len(adj)))
    for v,degree in order:
        check(degree==sum(w in remaining for w in adj[v]) and degree<=3,'peeling certificate')
        remaining.remove(v)
    result['graphs'][name]={'vertices':len(adj),'edges':expected_edges[name],'degeneracy':width,'peel_certificate':order,'3colour_search_nodes':n3,'4colour_search_nodes':n4,'4colour_witness':[c4[v] for v in range(len(adj))]}

x,y=s.symbols('x y');Q=(x*x+y*y)/2
D=lambda f,a:s.expand(s.diff(f,x)*s.diff(a,y)-s.diff(f,y)*s.diff(a,x))
controls=[('Keller_linear',x),('Keller_triangular',x+y*y),('Keller_composed',x+(y+x*x)**2),('smooth_nonKeller',x+x*x*y),('radial_quadratic',Q),('radial_quartic',Q*Q+Q)]
result['gaussian_controls']=[]
for name,f in controls:
    for d in (0,2,4,6,8):
        basis=[x**a*y**b for a in range(d+1) for b in range(d+1-a)]
        outputs=[s.Poly(D(f,m)-m*D(f,Q),x,y) for m in basis]
        mons=sorted(set(m for p in outputs for m in p.monoms()))
        mat=s.polys.matrices.DomainMatrix.from_Matrix(s.Matrix([[p.coeff_monomial(m) for p in outputs] for m in mons]))
        nullity=len(basis)-mat.rank();expected=d//2+1 if name.startswith('radial') else 0
        check(nullity==expected,'Gaussian nullity '+name+str(d))
        result['gaussian_controls'].append({'f':name,'degree_bound':d,'domain_dimension':len(basis),'nullity':nullity})

@lru_cache(None)
def gaussian_moment(a,b):
    if a%2 or b%2:return s.Integer(0)
    return s.factorial2(a-1)*s.factorial2(b-1)
result['direct_gaussian_pairing_controls']=[]
for name,f in controls:
    d=6; densities=[(a,b) for a in range(d+1) for b in range(d+1-a)]
    max_test=s.Poly(f,x,y).total_degree()+d
    tests=[(a,b) for a in range(max_test+1) for b in range(max_test+1-a)]
    derivs=[s.Poly(D(f,x**a*y**b),x,y).terms() for a,b in tests]
    residuals=[s.Poly(x**a*y**b*D(f,Q)-D(f,x**a*y**b),x,y).terms() for a,b in densities]
    direct=s.Matrix([[sum(coef*gaussian_moment(i+ha,j+hb) for (i,j),coef in terms) for ha,hb in densities] for terms in derivs])
    byparts=s.Matrix([[sum(coef*gaussian_moment(i+ta,j+tb) for (i,j),coef in terms) for terms in residuals] for ta,tb in tests])
    check(direct==byparts,'direct IBP '+name)
    nullity=len(densities)-s.polys.matrices.DomainMatrix.from_Matrix(direct).rank()
    check(nullity==(4 if name.startswith('radial') else 0),'direct Gaussian nullity '+name)
    result['direct_gaussian_pairing_controls'].append({'f':name,'density_degree':d,'test_degree':int(max_test),'number_of_tests':len(tests),'nullity':nullity})

phi=(1+s.sqrt(5))/2
radius2=s.simplify(1/(4*s.sin(s.pi/5)**2))
shadow2=s.simplify(phi**(-4)*radius2)
check(s.simplify(shadow2-(5-2*s.sqrt(5))/5)==0,'regular pentagon exact radius')
check(shadow2!=1,'regular pentagon no common unit apex')
result['regular_pentagon']={'original_radius_squared':str(radius2),'shadow_radius_squared':str(shadow2),'no_common_unit_apex':True}
core={}
for n in range(1,25):
    core[n]=[marked for marked in range(min(n,9)+1) if 30*n-10*marked<=2*n*(n-1) and 6*marked<=(n-marked)*(n-marked-1)]
check(all(not core[n] for n in range(1,14)),'no core below fourteen')
check(core[14]==[6,7],'fourteen vertex necessary marked counts')
result['core_inequality_candidates']=core

c=s.Integer(2);t=s.symbols('t')
def ell_formula(a,b):
    k=2*b+1-a
    if 0<=k<=b:return (-1)**(k+1)*s.binomial(b,k)*c**(a-b-1)
    return s.Integer(0)
def ell(poly):
    return s.expand(sum(coef*ell_formula(a,b) for (a,b),coef in s.Poly(poly,x,y).terms()))
f=x+x*x*y;count=0
for a in range(9):
    for b in range(9-a):
        direct=s.residue(-t**a*((c-t)/t**2)**b/t**2,t,0)
        check(direct==ell_formula(a,b),'residue formula')
        check(ell(D(f,x**a*y**b))==0,'residue annihilator')
        count+=1
check(ell(1)==0 and ell(x)==-1,'nonzero period')
check(s.simplify(D(f,1/x))==1,'rational slice')
result['residue_control']={'f':str(f),'c':str(c),'monomials_checked':count,'ell_1':str(ell(1)),'ell_x':str(ell(x)),'rational_slice':'1/x'}
result['status']='ALL EXACT CHECKS PASSED'
print(json.dumps(result,indent=2))
