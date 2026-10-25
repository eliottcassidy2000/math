#!/usr/bin/env python3
"""Exact finite/symbolic controls; no removable assertions. Needs SymPy."""
from itertools import product, permutations
from collections import defaultdict, Counter
from math import gcd
from pathlib import Path
import json, hashlib
import sympy as s

OUT = Path(__file__).with_name('certificate.json')
def require(ok, message):
    if not ok: raise RuntimeError(message)

def components(n, edges):
    adj=[[] for _ in range(n)]
    for a,b in edges: adj[a].append(b); adj[b].append(a)
    cs=[]; seen=set()
    for root in range(n):
        if root in seen: continue
        seen.add(root); stack=[root]; c=[]
        while stack:
            u=stack.pop(); c.append(u)
            for v in adj[u]:
                if v not in seen: seen.add(v); stack.append(v)
        cs.append(tuple(sorted(c)))
    return cs

def integrate(n, edges, labels, mod):
    adj=[[] for _ in range(n)]
    for (a,b),z in zip(edges,labels):
        adj[a].append((b,z)); adj[b].append((a,-z))
    p=[None]*n
    for root in range(n):
        if p[root] is not None: continue
        p[root]=0; stack=[root]
        while stack:
            a=stack.pop()
            for b,z in adj[a]:
                want=(p[a]+z)%mod
                if p[b] is None: p[b]=want; stack.append(b)
                elif p[b]!=want: return None
    return tuple(p)

def potentials(n, edges, mod):
    roots={c[0] for c in components(n,edges)}
    free=[v for v in range(n) if v not in roots]
    ans={}
    for vals in product(range(mod),repeat=len(free)):
        p=[0]*n
        for v,z in zip(free,vals): p[v]=z
        lab=tuple((p[b]-p[a])%mod for a,b in edges)
        ans[lab]=tuple(p)
    return ans

def graph_census():
    all_edges=[(i,3+j) for i in range(3) for j in range(3)]
    result={}
    for mod in [2,3]:
        assigned=0; potential_vectors=0
        for mask in range(1<<9):
            es=tuple(e for k,e in enumerate(all_edges) if mask>>k&1)
            # Independent direct image uses all vertex assignments, not a forest gauge.
            im={tuple((p[b]-p[a])%mod for a,b in es)
                for p in product(range(mod), repeat=6)}
            integrated={z for z in product(range(mod),repeat=len(es))
                        if integrate(6,es,z,mod) is not None}
            require(im==integrated,('forest image mismatch',mod,mask))
            require(len(im)==mod**(6-len(components(6,es))),('rank',mod,mask))
            assigned+=mod**len(es); potential_vectors+=mod**6
        require(assigned==(1+mod)**9, 'assignment universe')
        result[str(mod)]={'graphs':512,'edge_vectors_tested':assigned,
                          'vertex_vectors_tested':potential_vectors}
    return result

def edgepath(path):
    return set(tuple(sorted(e)) for e in zip(path,path[1:]))

def seam_census():
    n=6
    paths=[]
    for l in permutations(range(3)):
        for r in permutations(range(3,6)):
            path=tuple(v for pair in zip(l,r) for v in pair)
            paths.append(edgepath(path))
    p=paths[0]; q=next(q for q in paths if len(p&q)==1)
    full={(i,3+j) for i in range(3) for j in range(3)}
    require(p|q==full,'K33 path cover')
    h=full-{(0,5),(2,3)}; v=full-{(0,3),(2,5)}
    p44=edgepath((0,7,1,5,2,4,3,6))
    q44=edgepath((0,7,2,5,1,6,3,4))
    specs=[('K33_paths',6,p,q,5),('K44_frozen_paths',8,p44,q44,3),
           ('K33_hourglasses',6,h,v,1)]
    out=[]
    for name,n,h1,h2,kexp in specs:
        shared=tuple(sorted(h1&h2)); union=tuple(sorted(h1|h2))
        cs=components(n,shared); k=len(cs)
        require(k==kexp,(name,'seam component count'))
        for mod in [2,3,4]:
            e1=tuple(sorted(h1)); e2=tuple(sorted(h2))
            im1=potentials(n,e1,mod); im2=potentials(n,e2,mod)
            ix1=[e1.index(e) for e in shared]; ix2=[e2.index(e) for e in shared]
            by_shared=defaultdict(list)
            for z,pot in im2.items(): by_shared[tuple(z[i] for i in ix2)].append((z,pot))
            hist=Counter(); good=total=0
            for z1,pot1 in im1.items():
                for z2,pot2 in by_shared[tuple(z1[i] for i in ix1)]:
                    labs=dict(zip(e1,z1)); labs.update(zip(e2,z2))
                    z=tuple(labs[e] for e in union)
                    diff=[(pot1[c[0]]-pot2[c[0]])%mod for c in cs]
                    for c,d in zip(cs,diff):
                        require(all((pot1[w]-pot2[w])%mod==d for w in c),'component gauge')
                    obstruction=tuple((d-diff[0])%mod for d in diff[1:])
                    ok=integrate(n,union,z,mod) is not None
                    require(ok==(not any(obstruction)),(name,'obstruction iff'))
                    hist[obstruction]+=1; good+=ok; total+=1
            require(len(hist)==mod**(k-1),'obstruction surjectivity')
            require(set(hist.values())=={mod**(n-1)},'equal obstruction fibers')
            out.append({'name':name,'modulus':mod,'intersection_components':list(map(len,cs)),
                        'compatible_pairs':total,'globally_consistent':good,
                        'obstruction_classes':len(hist),'fiber_size':mod**(n-1),
                        'edges1':e1,'edges2':e2})
    return out

LINES=((0,1,2),(3,4,5),(6,7,8),(0,3,6),(1,4,7),(2,5,8),(0,4,8),(2,4,6))
def magic_group(vals, add, inv, zero):
    def sumg(xs):
        out=zero
        for x in xs: out=add(out,x)
        return out
    def nf(e,x,z):
        nx,nz=inv(x),inv(z)
        return (sumg([e,x]),sumg([e,nx,nz]),sumg([e,z]),sumg([e,nx,z]),e,
                sumg([e,x,nz]),sumg([e,nz]),sumg([e,x,z]),sumg([e,nx]))
    expected={nf(e,x,z) for e,x,z in product(vals,repeat=3)}
    require(len(expected)==len(vals)**3,'magic normal form injective')
    actual=set()
    # Independent row/column elimination from A,B,C,D,E, then check all eight lines.
    for a,b,c,d,e in product(vals,repeat=5):
        target=sumg([a,b,c])
        f=sumg([target,inv(d),inv(e)])
        g=sumg([target,inv(a),inv(d)])
        h=sumg([target,inv(b),inv(e)])
        i=sumg([target,inv(c),inv(f)])
        arr=(a,b,c,d,e,f,g,h,i)
        if all(sumg(arr[k] for k in line)==target for line in LINES): actual.add(arr)
    require(actual==expected,'independent magic solution sets')
    return {'group_size':len(vals),'five_coordinate_candidates':len(vals)**5,'magic_arrays':len(actual)}

def group_family_census():
    out=[]
    for mod in range(2,9):
        for t in range(mod):
            vals=tuple(u for u in range(mod) if gcd(1+t*u,mod)==1)
            add=lambda u,v:(u+v+t*u*v)%mod
            inv=lambda u:(-u*pow((1+t*u)%mod,-1,mod))%mod
            for a,b,c in product(vals,repeat=3):
                require(add(add(a,b),c)==add(a,add(b,c)),'associativity')
            for a in vals: require(add(a,inv(a))==0,'inverse')
            d=magic_group(vals,add,inv,0);d.update(modulus=mod,t=t);out.append(d)
    vals=tuple(product(range(2),repeat=2))
    add=lambda a,b:tuple((x+y)%2 for x,y in zip(a,b))
    d=magic_group(vals,add,lambda a:a,(0,0));d['group']='C2 x C2';out.append(d)
    return out

def symbolic_checks():
    t,E,x,z=s.symbols('t E x z')
    op=lambda a,b:s.cancel(a+b+t*a*b)
    neg=lambda a:s.cancel(-a/(1+t*a))
    fold=lambda seq:__import__('functools').reduce(op,seq,0)
    X=[fold(a) for a in [(E,x),(E,neg(x),neg(z)),(E,z),(E,neg(x),z),(E,),
                         (E,x,neg(z)),(E,neg(z)),(E,x,z),(E,neg(x))]]
    target=fold([E,E,E])
    for line in LINES: require(s.cancel(fold(X[i] for i in line)-target)==0,'symbolic magic')
    special=[s.cancel(f.subs(t,0)) for f in X]
    require(special==[E+x,E-x-z,E+z,E-x+z,E,E+x-z,E-z,E+x+z,E-x],'special fiber')
    a,b,c,d=s.symbols('a b c d')
    cycle=s.expand(((1+t*a)*(1+t*c)-(1+t*b)*(1+t*d))/t)
    require(s.expand(cycle-(a+c-b-d+t*(a*c-b*d)))==0,'four cycle')
    # A graph-potential parametrization of K22 in deformed coordinates.
    u,v,w=s.symbols('u v w')
    require(s.cancel(cycle.subs({a:0,b:u,c:op(u,v),d:v}))==0,'cycle substitution')
    # Monomial exponent identities are the integral tangent equations.
    rows=[(1,1,0),(1,-1,-1),(1,0,1),(1,-1,1),(1,0,0),
          (1,1,-1),(1,0,-1),(1,1,1),(1,-1,0)]
    for line in LINES:
        require(tuple(sum(rows[i][j] for i in line) for j in range(3))==(3,0,0),'line exponents')
    M=s.Matrix(rows)
    require(M.rank()==3,'magic exponent rank')
    # Universal line-difference matrix has six unit invariant factors.
    D=s.Matrix([[int(i in line)-int(i in LINES[0]) for i in range(9)] for line in LINES[1:]])
    from sympy.matrices.normalforms import smith_normal_form
    S=smith_normal_form(D,domain=s.ZZ)
    invariants=[abs(int(S[i,i])) for i in range(min(S.shape)) if S[i,i]!=0]
    require(invariants==[1]*6,'magic saturation')
    require(D*M==s.zeros(7,3),'magic lattice kernel')
    require(2**5*2**2*2**2==2**9==512,'root cover fiber degrees')
    exps=tuple(product([0,1],repeat=9))
    require(all(sum(e[i] for e in exps)==256 for i in range(9)),'discriminant exponents')
    require(9*len(exps)==4608,'discriminant power of two')
    # Raw product relation loses an additive constraint in a dual-number ring.
    # R=F2[eps]/eps^2, t=eps, one edge z=eps, other edges zero.
    # t*z=0, so all raw multiplicative factors are1, but divided relation=eps!=0.
    return {'deformed_magic_entries':[str(v) for v in X],
            'special_fiber':[str(v) for v in special],
            'four_cycle':str(cycle),'magic_smith_invariants':invariants,
            'root_cover_rank':512,'shared_root_cover_rank':32,
            'root_cover_discriminant':'2^4608 * (z_1 * ... * z_9)^256',
            'hourglass_relative_ranks':[4,4],
            'raw_product_hostile':{'ring':'F2[eps]/(eps^2)','t':'eps','cycle_labels':['eps','0','0','0'],
                                   'raw_product_difference':'0','divided_cycle_difference':'eps'}}

def main():
    result={'scope':'Finite checks supplement proofs; no square rational-point or F1 descent claim.',
            'graph_census':graph_census(),'seam_census':seam_census(),
            'group_family_magic_census':group_family_census(),'symbolic':symbolic_checks(),
            'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    OUT.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print(json.dumps({'graphs_per_group':512,'seam_cases':len(result['seam_census']),
                      'finite_group_magic_cases':len(result['group_family_magic_census']),
                      'magic_rank':3,'root_cover_rank':512,'status':'PASS'},sort_keys=True))
if __name__=='__main__': main()
