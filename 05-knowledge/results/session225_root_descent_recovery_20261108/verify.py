"""Fresh reconstruction; Python standard library; no external data or removable checks."""
from itertools import product
from math import gcd
from fractions import Fraction
import json

def require(ok, msg):
    if not ok: raise RuntimeError(msg)

def move(s, c):
    D,E,U,d=s
    if c=='O': return D,0,2*U,d
    if c=='Z': return D-E+U,d,2*U-d,d
    if c=='J': return D-E+U,D-E,D+U,D
    raise ValueError(c)

def run(s,w):
    for c in w: s=move(s,c)
    return s

def mat(w):
    p=run((1,1,1,1),w);q=run((0,-1,0,-1),w)
    require(p[0]==p[2] and p[1]==p[3] and q[0]==q[2] and q[1]==q[3], 'chamber return')
    return p[0],q[0],p[0]-p[1],q[0]-q[1]

def digital(B):
    bits=[int(c=='Z') for c in B];n=len(bits);L=2**n
    q=sum(b*2**(n-i-1) for i,b in enumerate(bits))
    a=sum(b*2**i for i,b in enumerate(bits))
    c=sum(bits[i]*bits[j]*2**(i-1-j) for i in range(n) for j in range(i))
    c+=sum(bits[i]*bits[i-1] for i in range(1,n))
    i=bits[0]
    return L+a+1-q-i-c,q+i+c,L-q,q

def literal(w):
    rows=[2,0];root=0
    for c in w:
        v=len(rows);rows.append(0)
        for u in range(v):
            if c=='O' or u!=root: rows[u]|=1<<v
            else: rows[v]|=1<<u
        if c=='J':root=v
    return rows,root

def profile(rows,root):
    n=len(rows);full=(1<<n)-1
    pred=[sum(1<<u for u in range(n) if rows[u]>>v&1) for v in range(n)]
    ends=[[0]*n for _ in range(full+1)];starts=[[0]*n for _ in range(full+1)]
    H=[0]*(full+1);S=[0]*(full+1);H[0]=1
    for mask in range(1,full+1):
        bs=mask
        while bs:
            bit=bs&-bs;v=bit.bit_length()-1;bs-=bit;old=mask^bit
            if not old:ends[mask][v]=1;starts[mask][v]=int(v==root)
            else:
                ps=old&pred[v]
                while ps:
                    ubit=ps&-ps;u=ubit.bit_length()-1;ps-=ubit
                    ends[mask][v]+=ends[old][u];starts[mask][v]+=starts[old][u]
        H[mask]=sum(ends[mask]);S[mask]=sum(starts[mask])
    U=sum(S[m]*H[full^m] for m in range(1,full+1) if m>>root&1)
    return H[full],ends[full][root],U,H[full^(1<<root)]

def reachable(start,adj):
    seen={start};todo=[start]
    while todo:
        v=todo.pop()
        for u in adj(v):
            if u not in seen:seen.add(u);todo.append(u)
    return seen

def main():
    matrices={};blocks=0
    for n in range(1,12):
        for p in product('OZ',repeat=n-1):
            B=''.join(p)+'O';T=mat(B+'J');blocks+=1
            require(T==digital(B),('digital',B))
            A,Bb,C,D=T;L=2**n
            require(min(T)>=0 and C+D==L and L<A+Bb<2*L,('rows',B))
            if D:require(Bb>0 and C>0,('irreducible',B))
            if A*D-Bb*C>1:matrices.setdefault(T,B+'J')
    vertices={(a,-b) for b in range(2,81) for a in range(1,b) if gcd(a,b)==1}
    edges={v:[] for v in vertices};rev={v:[] for v in vertices}
    for T,w in matrices.items():
        A,B,C,D=T
        for v in vertices:
            x,y=v;u=(x*A+y*C,x*B+y*D)
            if u in vertices:edges[v].append((u,w,T));rev[u].append(v)
    root=(3,-4)
    comp=reachable(root,lambda v:(e[0] for e in edges[v]))&reachable(root,lambda v:iter(rev[v]))
    kept=sorted((v,u,w,T) for v in comp for u,w,T in edges[v] if u in comp)
    branch=max(sum(u in comp for u,w,T in edges[v]) for v in comp)
    mindet=min(T[0]*T[3]-T[1]*T[2] for v,u,w,T in kept)
    require((blocks,len(comp),len(kept),branch,mindet)==(2047,21,67,5,6),'graph summary')
    expected={(3,-4),(5,-6),(5,-7),(7,-8),(7,-9),(7,-10),(11,-12),(11,-13),(11,-14),(11,-15),(13,-16),(13,-19),(15,-16),(15,-17),(15,-19),(19,-22),(19,-25),(23,-32),(25,-32),(31,-38),(31,-45)}
    require(comp==expected,'component vertices')
    families=[('OOZOZOJ',(73,12,54,10),(3,-4),82,15,-1),('OZOOZZOOJ',(217,90,180,76),(5,-6),292,21,1),('OOOZOOZOOJ',(545,40,476,36),(7,-8),580,63,3)]
    family_checks=[]
    for w,T,l,radix,kappa,I in families:
        require(mat(w)==T,('family matrix',w));A,B,C,D=T;x,y=l
        require((x*A+y*C,x*B+y*D)==l,'invariant')
        require((A+C-radix,B+D-radix)==(kappa*x,kappa*y),'readout')
        s=run((1,0,2,1),'OJ');N=9;terms=[]
        for m in range(11):
            h,e,U,d=s;g=h-e
            require(h==U and e==d and x*h+y*g==I and h+g==N,'iterate')
            require(move(s,'J')[0]==N,'closing J')
            terms.append(N);N=radix*N+kappa*I;s=run(s,w)
        family_checks.append({'word':w,'matrix':T,'terms':terms})
    literal_checks=0
    for n in range(5):
        for w in map(''.join,product('OZJ',repeat=n)):
            require(profile(*literal(w))==run((1,0,2,1),w),('literal profile',w));literal_checks+=1
    closed=[]
    for w,expected_H in [('OJJ',9),('OJOOZOZOJJ',723),('OJOZOOZZOOJJ',2649),('OJOOOZOOZOOJJ',5409)]:
        p=profile(*literal(w));require(p==run((1,0,2,1),w) and p[0]==expected_H,('closed witness',w))
        closed.append({'word':w,'order':len(w)+2,'H':p[0]})
    # Hostile controls: omission of closing J and noninvariance of scalar alone.
    require(run((1,0,2,1),'OJ')[0]==5,'unclosed normalization')
    require(run((1,0,2,1),'OJJ')[0]==9,'closed normalization')
    # Continuation controls, not substitutes for its all-order proof.
    cone_checks=0;boundary=[]
    for n in range(8):
        for w in map(''.join,product('OZJ',repeat=n)):
            s=run((1,0,2,1),w);D,E,U,d=s;x=D-E;W=U-d
            require(d>0 and W>0 and 0<=E<=d and 0<x<=W,('cone',w));cone_checks+=1
            h,e,hh,ee=run(s,'OJ')
            require(h==hh and e==ee and 3*(h-e)>=2*h,('ratio',w))
            if e==1:
                require(set(w)<=set('O') and 2*h-e==2**(n+3)+1,('boundary',w))
                boundary.append(2*h-e)
    for a,negative_b in comp:
        b=-negative_b
        for c in (-1,1):require((b+c)%(a+b)!=0,'unit-line hostile residue')
    return {'status':'PASS','blocks':blocks,'graph':{'vertices':sorted(comp),'edges':kept,'maximum_branching':branch,'minimum_determinant':mindet},'families':family_checks,'literal_profile_checks':literal_checks,'closed_literal_witnesses':closed,'cone_controls':cone_checks,'boundary_counts':boundary,'scope':'Fresh reconstruction, finite controls; all-order claims require the accompanying proof. No cutoff or coverage claim.'}

if __name__=='__main__': print(json.dumps(main(),sort_keys=True,indent=2))
