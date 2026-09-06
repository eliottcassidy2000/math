"""Independent Boolean-sector audit: literal toggles, integer ranks, wedge vectors.

Python + SymPy. No producer/repository imports, all checks remain live under -O.
"""
from collections import Counter, deque
from fractions import Fraction as Q
from itertools import combinations, combinations_with_replacement, permutations
from math import factorial
import hashlib,json,sys
import sympy as s
sys.stdout.reconfigure(newline='\n')
GATES=0
def need(ok,label):
    global GATES
    GATES+=1
    if not ok:raise RuntimeError(label)
def states(n,c,M):
    return [t for t in combinations_with_replacement(range(3,M+1),c) if sum(t)<=n]
def edges(t):
    E=set();start=0
    for length in t:
        for j in range(length):
            E.add(tuple(sorted((start+j,start+(j+1)%length))))
        start+=length
    return E
def classify(E,n):
    adj=[set() for _ in range(n)]
    for u,v in E:adj[u].add(v);adj[v].add(u)
    if any(len(x) not in (0,2) for x in adj):return None
    seen=set();lengths=[]
    for root in range(n):
        if root in seen or not adj[root]:continue
        component={root};stack=[root]
        while stack:
            u=stack.pop()
            for v in adj[u]-component:component.add(v);stack.append(v)
        seen.update(component);lengths.append(len(component))
    return tuple(sorted(lengths))
def token(t):return tuple(l-3+i for i,l in enumerate(t,1))
def matrix(n,c,M):
    st=states(n,c,M)
    return st,s.Matrix([[int(sum(abs(x-y) for x,y in zip(a,b))==1)
                         for b in st] for a in st])
def volume(t,n):
    denominator=factorial(n-sum(t))
    for length,m in Counter(t).items():denominator*=(2*length)**m*factorial(m)
    return factorial(n)//denominator
def literal(n,c,M):
    st,A=matrix(n,c,M);index={t:i for i,t in enumerate(st)}
    multi=s.zeros(len(st));checks=0
    for i,t in enumerate(st):
        E=edges(t)
        for tri in combinations(range(n),3):
            E2=E ^ set(combinations(tri,2));target=classify(E2,n);checks+=1
            if target in index:multi[i,index[target]]+=1
    for i,t in enumerate(st):
        for j,u in enumerate(st):
            need(bool(multi[i,j])==bool(A[i,j]),'literal XOR versus length move')
            S,T=token(t),token(u)
            diff=set(S)^set(T)
            adjacency=len(diff)==2 and abs(min(diff)-max(diff))==1
            need(bool(A[i,j])==adjacency,'unweighted token move including repeats')
            need(multi[i,j]*volume(t,n)==multi[j,i]*volume(u,n),'actual orbit detailed balance')
            expected=0
            for length,m in Counter(t).items():
                values=list(t);values.remove(length)
                if length<M and tuple(sorted(values+[length+1]))==u:
                    expected+=m*length*(n-sum(t))
                if length>=4 and tuple(sorted(values+[length-1]))==u:expected+=m*length
            need(multi[i,j]==expected,'literal directed multiplicity')
    for start in range(len(st)):
        dist={start:0};todo=deque([start])
        while todo:
            i=todo.popleft()
            for j in range(len(st)):
                if A[i,j] and j not in dist:dist[j]=dist[i]+1;todo.append(j)
        for j,u in enumerate(st):
            need(dist[j]==sum(abs(x-y) for x,y in zip(st[start],u)),'cap-preserving distance')
    return {'n':n,'c':c,'M':M,'vertices':len(st),'toggles':checks,'edges':sum(A)//2}

# A small exact cyclotomic vector engine, using SymPy only for the modulus.
def ring(N):
    z=s.Symbol('z');p=s.Poly(s.cyclotomic_poly(2*(N+1),z),z)
    mod=list(reversed([int(x) for x in p.all_coeffs()]));degree=len(mod)-1
    zero=(0,)*degree
    def reduce(raw):
        raw=list(raw)+[0]*max(0,degree-len(raw))
        for j in range(len(raw)-1,degree-1,-1):
            c=raw[j]
            if c:
                for k in range(degree+1):raw[j-degree+k]-=c*mod[k]
        return tuple(raw[:degree])
    def add(a,b):return tuple(x+y for x,y in zip(a,b))
    def neg(a):return tuple(-x for x in a)
    def mul(a,b):
        raw=[0]*(2*degree-1)
        for i,x in enumerate(a):
            for j,y in enumerate(b):raw[i+j]+=x*y
        return reduce(raw)
    def power(k):return reduce([0]*(k%(2*(N+1)))+[1])
    def phi(j,b):return add(power(j*b),neg(power(-j*b)))
    def eig(J):
        out=zero
        for j in J:out=add(out,add(power(j),power(-j)))
        return out
    def det(rows,cols):
        out=zero
        for p in permutations(range(len(rows))):
            term=reduce([1]);inversions=sum(p[i]>p[j] for i in range(len(p)) for j in range(i+1,len(p)))
            for i in range(len(rows)):term=mul(term,phi(cols[p[i]],rows[i]))
            out=add(out,neg(term) if inversions%2 else term)
        return out
    return zero,add,mul,eig,det
def modes(N,c,only_zero=False):
    st=list(combinations(range(1,N+1),c));zero,add,mul,eig,det=ring(N)
    adjacency=[[j for j,T in enumerate(st) if len(set(S)^set(T))==2
                 and max(set(S)^set(T))-min(set(S)^set(T))==1] for S in st]
    tested=0;zeros=[]
    for J in st:
        lam=eig(J)
        if lam==zero:zeros.append(J)
        if only_zero and lam!=zero:continue
        vector=[det(S,J) for S in st]
        need(any(v!=zero for v in vector),'Slater vector nonzero')
        for i in range(len(st)):
            result=zero
            for j in adjacency[i]:result=add(result,vector[j])
            need(result==mul(lam,vector[i]),'literal exact cyclotomic eigenvector')
        tested+=1
    return {'N':N,'c':c,'tested_vectors':tested,'zero_modes':zeros}
def main():
    bank=[literal(*p) for p in ((9,2,6),(11,3,5),(15,3,5),(15,4,6),(24,3,8))]
    spectra=[modes(N,c) for N,c in ((3,2),(4,2),(5,2),(6,3))]
    spectra.append(modes(8,3,True))
    t,B=matrix(24,3,8)
    need(len(t)==56 and len(t)-B.rank()==2,'independent integer rank of resonant sector')
    need(sum((-1)**sum(x) for x in t)==0,'native parity is balanced')
    need(spectra[-1]['zero_modes']==[(1,5,7),(2,4,8)],'complete critical zero-mode list')
    # Fresh exact capped-sector ranks are finite controls, not all-cap extrapolation.
    caps=[]
    for n in range(19,25):
        t,A=matrix(n,3,8);null=len(t)-A.rank();index=abs(sum((-1)**sum(x) for x in t))
        need(null>=index,'parity lower bound on every capped matrix')
        caps.append((n,len(t),index,null))
    # Build the additive Laplacian directly in the token basis, with occupied-edge correction.
    N,c=5,2;st=list(combinations(range(1,N+1),c));D=s.zeros(len(st));L=s.zeros(len(st))
    for i,S in enumerate(st):
        D[i,i]=sum(1 if j in (1,N) else 2 for j in S)
        moves=0
        for j,T in enumerate(st):
            diff=set(S)^set(T)
            if len(diff)==2 and max(diff)-min(diff)==1:D[i,j]=-1;L[i,j]=-1;moves+=1
        L[i,i]=moves
        occupied=sum(j+1 in S for j in S)
        need(D[i,i]-L[i,i]==2*occupied,'exact missing diagonal potential')
    need(D.det()!=0 and L.det()==0,'free-sum Laplacian hostile')
    pendant=[]
    for c in range(1,5):
        for d in range(1,5):
            if c*d<2:continue
            M=d+3;t,A=matrix(c*M,c,M)
            u,v=len(t)-1,len(t)-2
            need(sum(A[u,j] for j in range(len(t)))==1 and A[u,v]==1,'unique top leaf and neighbor')
            keep=list(range(len(t)-2));target,_=matrix(c*M-2,c,M)
            need([t[i] for i in keep]==target,'exact two-state resource-cap deletion')
            order=[u,v]+keep;Ar=A.extract(order,order);Qm=s.eye(len(t))
            for j in range(2,len(t)):Qm[0,j]=-Ar[1,j]
            want=s.diag(s.Matrix([[0,1],[1,0]]),A.extract(keep,keep))
            need(Qm.T*Ar*Qm==want,'integral unitriangular pendant congruence')
            need(Qm.det()==1,'unimodular change preserves integral cokernel')
            pendant.append((c,M,len(t)))
    payload={'literal':bank,'eigenvectors':spectra,'caps':caps,'pendant':pendant}
    print('STATUS: PASS; independent analytic review plus exact native controls')
    for row in bank:print('LITERAL:',row)
    for row in spectra:print('MODES:',row)
    print('CAPPED FINITE BANK (n,vertices,parity index,nullity):',caps)
    print('PENDANT INTEGRAL CONGRUENCES (c,M,vertices):',pendant)
    print('SCOPE: induced Boolean triangle adjacency; no ambient or multiplicity/Laplacian transfer')
    print('SEMANTIC SHA256:',hashlib.sha256(json.dumps(payload,sort_keys=True,default=str).encode()).hexdigest())
    print('ACTIVE GATES:',GATES)
if __name__=='__main__':main()
