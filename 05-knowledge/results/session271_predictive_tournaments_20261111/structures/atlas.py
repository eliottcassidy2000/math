#!/usr/bin/env python3
"""Small reference classifier + polynomial spectral substitution, not a census."""
import itertools, json, random
from pathlib import Path
import sympy as s
z = s.Symbol('z')

def tournament(n, mask):
    A = s.zeros(n)
    for b, (i,j) in enumerate(itertools.combinations(range(n), 2)):
        A[i,j], A[j,i] = ((1,0) if (mask>>b)&1 else (0,1))
    return A

def sub(Q, children):
    sizes = [len(A.tolist()) for A in children]
    offsets = [sum(sizes[:i]) for i in range(len(sizes))]
    A = s.diag(*children)
    for i,j in itertools.permutations(range(len(children)),2):
        if Q[i,j]:
            for a in range(sizes[i]):
                for b in range(sizes[j]):
                    A[offsets[i]+a,offsets[j]+b]=1
    return A

def chi(A): return A.charpoly(z).as_expr().expand()
def spectral_compose(Q, fs):
    gs=[((-1)**s.degree(f,z)*f.subs(z,-z-1)-f).expand() for f in fs]
    return (s.diag(*fs)-Q*s.diag(*gs)).det().expand()

def is_module(A, M):
    return all(len({int(A[v,u]) for u in M})==1 for v in range(A.rows) if v not in M)

def strong_children(A):
    n=A.rows
    mods=[frozenset(i for i in range(n) if mask>>i&1) for mask in range(1,1<<n)]
    mods=[M for M in mods if is_module(A,M)]
    strong=[M for M in mods if not any(M&N and M-N and N-M for N in mods)]
    proper=[M for M in strong if len(M)<n]
    return [M for M in proper if not any(M<N for N in proper)]

def canonical_tree(A):
    """Exponential subsets/permutations: transparent reference, small n only."""
    if A.rows==1: return {'kind':'leaf'}
    blocks=strong_children(A)
    if set().union(*blocks)!=set(range(A.rows)) or sum(map(len,blocks))!=A.rows:
        raise ValueError('strong module partition failed')
    blocks=[sorted(M) for M in blocks]
    Q=s.Matrix(len(blocks),len(blocks),lambda i,j:0 if i==j else A[blocks[i][0],blocks[j][0]])
    kids=[canonical_tree(A.extract(B,B)) for B in blocks]
    scores=[sum(Q.row(i)) for i in range(Q.rows)]
    linear=sorted(scores)==list(range(Q.rows))
    orders=[tuple(sorted(range(Q.rows),key=lambda i:-scores[i]))] if linear else itertools.permutations(range(Q.rows))
    forms=[]
    for p in orders:
        bits=''.join(str(Q[p[i],p[j]]) for i,j in itertools.combinations(range(Q.rows),2))
        item={'kind':'linear' if linear else 'prime','quotient_bits':bits,'children':[kids[i] for i in p]}
        forms.append(json.dumps(item,sort_keys=True,separators=(',',':')))
    return json.loads(min(forms))

def realize(node):
    if node['kind']=='leaf': return s.zeros(1)
    Q=tournament(len(node['children']),sum(int(b)<<i for i,b in enumerate(node['quotient_bits'])))
    return sub(Q,[realize(c) for c in node['children']])

def aut_order(node):
    if node['kind']=='leaf': return 1
    kids=node['children']; k=len(kids)
    Q=tournament(k,sum(int(b)<<i for i,b in enumerate(node['quotient_bits'])))
    colors=[json.dumps(c,sort_keys=True) for c in kids]
    stab=sum(all(colors[i]==colors[p[i]] for i in range(k)) and all(Q[i,j]==Q[p[i],p[j]] for i in range(k) for j in range(k)) for p in itertools.permutations(range(k)))
    return stab*s.prod(aut_order(c) for c in kids)

def direct_aut(A):
    return sum(all(A[i,j]==A[p[i],p[j]] for i in range(A.rows) for j in range(A.rows)) for p in itertools.permutations(range(A.rows)))

def boundary(A, vertices):
    U=s.Matrix.hstack(s.ones(A.rows,1),*[s.eye(A.rows)[:,v] for v in vertices])
    return U.T*(z*s.eye(A.rows)-A).inv()*U

def main():
    rng=random.Random(271)
    graphs=[tournament(n,rng.randrange(1<<(n*(n-1)//2))) for n in range(1,7) for _ in range(4)]
    checks={'canonical_roundtrips':0,'automorphism_orders':0,'coronal_identities':0,'substitution_polynomials':0}
    for A in graphs:
        node=canonical_tree(A)
        if canonical_tree(realize(node))!=node: raise ValueError('tree roundtrip')
        checks['canonical_roundtrips']+=1
        if aut_order(node)!=direct_aut(A): raise ValueError('automorphism recursion')
        checks['automorphism_orders']+=1
        f=chi(A)
        r=(s.ones(1,A.rows)*(z*s.eye(A.rows)-A).inv()*s.ones(A.rows,1))[0]
        if s.cancel(1+r-(-1)**A.rows*f.subs(z,-z-1)/f)!=0: raise ValueError('coronal')
        checks['coronal_identities']+=1
    C3=tournament(3,5)
    for Q in [tournament(2,1),C3,tournament(4,17)]:
        for trial in range(3):
            kids=[tournament(n,rng.randrange(1<<(n*(n-1)//2))) for n in [rng.randrange(1,4) for _ in range(Q.rows)]]
            if spectral_compose(Q,[chi(A) for A in kids])!=chi(sub(Q,kids)): raise ValueError('substitution')
            checks['substitution_polynomials']+=1
    # The same TT2 with two different chosen roots: reverse root -> sink.
    A=tournament(2,1); B=s.zeros(1); J=sub(tournament(2,1),[A,B])
    seams=[]
    for a in range(2):
        T=J.copy(); T[a,2]=0; T[2,a]=1
        E=s.Matrix.hstack(s.eye(3)[:,a],s.eye(3)[:,2]); K=s.Matrix([[0,-1],[1,0]])
        R=E.T*(z*s.eye(3)-J).inv()*E
        if s.cancel(chi(J)*(s.eye(2)-K*R).det()-chi(T))!=0: raise ValueError('seam')
        seams.append({'root':a,'polynomial':str(chi(T))})
    examples={}
    for name,A in [('linear3',tournament(3,7)),('cyclic3',C3),('cyclic_tower9',sub(C3,[C3]*3)),('unequal_cyclic4',sub(C3,[tournament(2,1),s.zeros(1),s.zeros(1)]))]:
        node=canonical_tree(A); examples[name]={'tree':node,'chi':str(chi(A)),'automorphism_order':int(aut_order(node))}
    result={'status':'FINITE-EXACT controls for new reference implementation, not a full census','checks':checks,'seam_hostile':seams,'examples':examples,'prime_union_bound':{str(n):str(s.Rational(4*s.binomial(n,2)+4*n+64,2**n)) for n in [16,32,64,128]}}
    print(json.dumps(result,indent=2))
if __name__=='__main__': main()
