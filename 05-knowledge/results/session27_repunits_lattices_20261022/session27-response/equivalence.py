"""Exact rational reachability test for full labelled attachment response.
Standalone implementation; small k only. No graph classification claim.
"""
from itertools import permutations
from fractions import Fraction
from collections import deque
from math import factorial

def core(k,mask):
    Q=[[0]*k for _ in range(k)];t=0
    for u in range(k):
        for v in range(u+1,k):
            a,b=(u,v) if mask>>t&1 else (v,u)
            Q[a][b]=1;t+=1
    return Q

def validate(Q):
    k=len(Q)
    if any(len(r)!=k for r in Q) or any(x not in (0,1) for r in Q for x in r) or any(Q[i][i] for i in range(k)) or any(Q[i][j]+Q[j][i]!=1 for i in range(k) for j in range(i)):
        raise ValueError('not a tournament')

def realization(k):
    if k==0:return [()],1,[[(0,)]],{0:1}
    ps=list(permutations(range(k)));t=3**k;size=2*t
    matrices=[]
    for a in range(1<<k):
        rows=[]
        for block,p in enumerate(ps):
            b=sum(((a>>u)&1)<<j for j,u in enumerate(p))
            for state in range(size):
                dest=[]
                for j in range(k):
                    d=state//(3**j)%3
                    if j and d==0 and not (b>>(j-1)&1):continue
                    nd=1+((b>>j)&1)
                    dest.append(block*size+state+(nd-d)*3**j)
                if state>=t or b>>(k-1)&1:dest.append(block*size+(state if state>=t else state+t))
                rows.append(tuple(dest))
        matrices.append(rows)
    return ps,size,matrices,{i*size:1 for i in range(len(ps))}

def final(Q,ps,size):
    k=len(Q)
    if k==0:return [1]
    return [int(all(state//3**j%3!=2 and (not j or state//3**j%3!=0 or Q[p[j-1]][p[j]]) for j in range(k))) for p in ps for state in range(size)]

def advance(v,rows):
    r={}
    for s,c in v.items():
        for t in rows[s]:r[t]=r.get(t,0)+c
    return r

def independent(v,basis):
    v={i:Fraction(x) for i,x in v.items() if x}
    for p,b in basis:
        c=v.get(p,0)
        if c:
            for j,x in b.items():
                z=v.get(j,0)-c*x
                if z:v[j]=z
                elif j in v:del v[j]
    if not v:return False
    p=min(v);c=v[p];v={j:x/c for j,x in v.items()}
    basis.append((p,v));basis.sort(key=lambda b:b[0]);return True

def equivalent(Q,R):
    """Return exact certificate summary; word is a distinguishing mask tuple.
    Equality exhausts a transition-invariant rational reachable row basis.
    No shortcut on Q==R: tests exercise the actual equivalence certification.
    """
    validate(Q);validate(R)
    if len(Q)!=len(R):raise ValueError('label sets must agree')
    k=len(Q);ps,size,matrices,alpha=realization(k)
    f=final(Q,ps,size);g=final(R,ps,size);delta=[x-y for x,y in zip(f,g)]
    basis=[];independent(alpha,basis);queue=deque([(alpha,())]);maxdepth=0
    while queue:
        v,w=queue.popleft();maxdepth=max(maxdepth,len(w))
        diff=sum(c*delta[s] for s,c in v.items())
        if diff:return {'equivalent':False,'word':w,'difference':diff,'basis_size':len(basis),'bound':len(f)-1}
        for a,rows in enumerate(matrices):
            new=advance(v,rows)
            if independent(new,basis):queue.append((new,w+(a,)))
    return {'equivalent':True,'reachable_dimension':len(basis),'max_basis_word_length':maxdepth,'bound':len(f)-1}

if __name__=='__main__':
    import argparse,json
    p=argparse.ArgumentParser();p.add_argument('k',type=int);p.add_argument('left',type=int);p.add_argument('right',type=int);a=p.parse_args()
    print(json.dumps(equivalent(core(a.k,a.left),core(a.k,a.right))))
