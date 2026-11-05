from fractions import Fraction
from collections import deque
from math import gcd
from pathlib import Path
import hashlib,json

P=Path(__file__).resolve().parent
def need(x,m):
    if not x: raise RuntimeError(m)
I=((1,0,0),(0,1,0),(0,0,1))
def mm(a,b,m=None):
    c=tuple(tuple(sum(a[i][k]*b[k][j] for k in range(3)) for j in range(3)) for i in range(3))
    return tuple(tuple(x%m for x in row) for row in c) if m else c
def det(a):
    return sum((1 if j%2==0 else -1)*a[0][j]*(a[1][(j+1)%3]*a[2][(j+2)%3]-a[1][(j+2)%3]*a[2][(j+1)%3])*(1 if j%2==0 else -1) for j in range(3))
def inverse(a,m=None):
    # Adjugate, not the submitted Gaussian elimination.
    d=det(a); z=[]
    for i in range(3):
        row=[]
        for j in range(3):
            rs=[x for x in range(3) if x!=j]; cs=[x for x in range(3) if x!=i]
            co=(-1)**(i+j)*(a[rs[0]][cs[0]]*a[rs[1]][cs[1]]-a[rs[0]][cs[1]]*a[rs[1]][cs[0]])
            row.append(co*pow(int(d),-1,m)%m if m else Fraction(co,d))
        z.append(tuple(row))
    return tuple(z)
def powm(a,n,m=None):
    out=I
    while n:
        if n&1: out=mm(out,a,m)
        a=mm(a,a,m); n//=2
    return out
def comm(a,b,m=None): return mm(mm(mm(a,b,m),inverse(a,m),m),inverse(b,m),m)
def E(i,j,t=1):return tuple(tuple(I[x][y]+t*(x==i and y==j) for y in range(3)) for x in range(3))
def mv(a,v,m=None):
    z=tuple(sum(a[i][j]*v[j] for j in range(3)) for i in range(3))
    return tuple(x%m for x in z) if m else z
def matrix(L):
    n=L-1;s=n.bit_length()
    g=1+sum((n>>(i+1))+((n>>(i+1))&1) for i in range(s) if not(n>>i&1))
    q=1+sum(1<<(s-i-1) for i in range(s) if not(n>>i&1))
    return ((g,-(L&1),L),(1,-int(L==1),0),(q,0,1<<s))
A=matrix(1)
def directions(t12,m=None):
    t32base=mm(comm(t12,mm(mm(powm(A,2,m),t12,m),inverse(powm(A,2,m),m),m),m),powm(t12,2,m),m)
    if m is None: return t32base
    # The caller supplies q separately for the inverse of -4.
    return t32base
def all_directions(t32,m):
    s=mm(mm(inverse(A,m),t32,m),A,m)
    t31=mm(comm(t32,s,m),t32,m)
    rr=mm(t31,inverse(t32,m),m)
    t12=mm(mm(mm(A,rr,m),inverse(A,m),m),inverse(t32,m),m)
    t21=mm(mm(inverse(t12,m),s,m),t12,m)
    t13=mm(mm(inverse(t31,m),mm(mm(A,t32,m),inverse(A,m),m),m),t31,m)
    return {(2,1):t32,(2,0):t31,(0,1):t12,(1,0):t21,(0,2):t13,(1,2):comm(t21,t13,m)}

def gadget(L):
    n=L-1;s=n.bit_length(); rows=[0]*(s+1)
    for i in range(s):
        for j in range(i+1,s):rows[i+1]|=1<<(j+1)
        if n>>i&1:rows[i+1]|=1
        else:rows[0]|=1<<(i+1)
    return rows
def seam(T,r,L):
    G=gadget(L);n=len(T);rows=[x|(((1<<len(G))-1)<<n) for x in T]+[x<<n for x in G]
    rows[r]^=1<<n;rows[n]|=1<<r
    return rows,n
def flip(T):
    f=(1<<len(T))-1
    return [f^x^(1<<i) for i,x in enumerate(T)]
def pair_graft(G,r):
    # Left triangle is Tr2 with a sole-neighbor birth, entries 0 and 2.
    T=[2,4,1]; n=3
    rows=[x|(((1<<len(G))-1)<<n) for x in T]+[x<<n for x in G]
    for a in (0,2):rows[a]^=1<<(n+r);rows[n+r]|=1<<a
    return rows

def main():
    counts=dict(relative=0,elementary=0,primitive_vectors=0,positive_inverse=0,seeds=0,seed_maps=0,targets=0,physical=0)
    P9=mm(matrix(9),inverse(matrix(11)));Q11=mm(matrix(11),inverse(matrix(13)))
    H=comm(P9,Q11)
    need(H==E(0,1,Fraction(3,4)),'exact rational shear')
    need(directions(E(0,1))==E(2,1,-4),'exact E32 commutator')
    need(all(t==E(i,j) for (i,j),t in all_directions(E(2,1),None).items()),'six exact rational elementary directions')
    for z in range(-50,52):
        v=(z,1,z)
        for t in [E(0,1,-z),E(2,1,-z),E(0,1),E(1,0,-1)]:v=mv(t,v)
        need(v==(1,0,0),'first exact seed reduction')
        v=(10*z-1,z,15*z)
        for t in [E(0,1,-10),E(1,0,z),E(2,0,15*z),((-1,0,0),(0,-1,0),(0,0,1))]:v=mv(t,v)
        need(v==(1,0,0),'second exact seed reduction')
        counts['seed_maps']+=2
    ratio=mm(matrix(5),inverse(matrix(7)))
    need(ratio==((1,Fraction(5,4),Fraction(-1,4)),(0,1,0),(0,2,1)) and ratio!=E(2,1,2),'false shortcut hostile')
    for m in range(3,256,2):
        q=m;a=0
        while q%3==0:q//=3;a+=1
        if q==1:continue
        d=3**max(a-1,0);k=d*((4*pow(3*d,-1,q))%q)
        hmod=comm(mm(matrix(9),inverse(matrix(11),m),m),mm(matrix(11),inverse(matrix(13),m),m),m)
        t12=powm(hmod,k,m)
        t32=powm(directions(t12,m),pow(-4,-1,q),m)
        dirs=all_directions(t32,m)
        for (i,j),g in dirs.items():
            need(tuple(tuple(x%q for x in row) for row in g)==tuple(tuple(x%q for x in row) for row in E(i,j)),'relative elementary q projection')
            mod3=3**a
            need(tuple(tuple(x%mod3 for x in row) for row in g)==tuple(tuple(x%mod3 for x in row) for row in I),'relative elementary identity on 3 power')
            counts['elementary']+=1
        counts['relative']+=1
    # Actual positive generator orbits, compared to every primitive vector.
    orbit_sizes=[]
    for q in (5,7,11,13,25,35):
        start=(3%q,1,3%q); seen={start}; todo=deque([start])
        mats=[matrix(L) for L in (1,9,11,13)]
        while todo:
            v=todo.popleft()
            for g in mats:
                w=mv(g,v,q)
                if w not in seen:seen.add(w);todo.append(w)
        exact={(x,y,z) for x in range(q) for y in range(q) for z in range(q) if gcd(gcd(gcd(x,y),z),q)==1}
        need(seen==exact,'entire primitive-vector positive orbit')
        counts['primitive_vectors']+=len(seen);orbit_sizes.append([q,len(seen)])
    for m in range(3,16,2):
        for L in (1,9,11,13):
            g=matrix(L);cur=I;previous=I;order=0
            while True:
                previous=cur;cur=mm(cur,g,m);order+=1
                if cur==I:break
                need(order<=m**9,'finite order bound')
            need(previous==inverse(g,m),'positive-power inverse')
            counts['positive_inverse']+=1
    for a in range(1,8):
        mod=3**a;sets=[set() for _ in range(4)]
        for j in range(1,2*3**(a-1)+1):
            z=pow(2,j,mod)+1;v=(z,1,z);w=mv(matrix(5),v)
            sets[0].add((53*v[0]-3*v[1]+63*v[2])%mod)
            sets[1].add((53*w[0]-3*w[1]+63*w[2])%mod)
            sets[2].add(v[0]%mod);sets[3].add(w[0]%mod);counts['seeds']+=1
        for ss,excluded in zip(sets,(2,0,1,0)):
            need(ss=={r for r in range(mod) if r%3!=excluded},'full exact 3-adic ray set')
    for r in range(-200,201):
        g=I
        for t in [E(1,0,-r),E(0,1,-1),E(1,0,2),E(0,1),E(1,0,17)]:g=mm(t,g)
        w=mv(g,(1,0,0));need(w==(2*r+3,35*r+53,0) and 53*w[0]-3*w[1]+63*w[2]==r,'signed target and order')
        counts['targets']+=1
    lines=[]
    for word in [(),(1,),(5,),(7,),(9,),(11,),(13,)]:
        T=[2,4,1];r=0;v=(3,1,3)
        for L in word:T,r=seam(T,r,L);v=mv(matrix(L),v)
        lines.append(' '.join(map(str,[len(T),v[0],1]+T)))
        T,r=seam(T,r,5);T=pair_graft(flip(T),r)
        want=53*v[0]-3*v[1]+63*v[2]
        lines.append(' '.join(map(str,[len(T),want,2]+T)));counts['physical']+=2
    (P/'graphs.txt').write_text('\n'.join(lines)+'\n')
    print(json.dumps(dict(status='PASS',counts=counts,positive_orbit_sizes=orbit_sizes,graphs_sha256=hashlib.sha256(('\n'.join(lines)+'\n').encode()).hexdigest()),sort_keys=True,indent=2))
if __name__=='__main__':main()
