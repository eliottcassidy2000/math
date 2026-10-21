"""Exact polynomial-basis finite-field tests; Python standard library only."""
import itertools, json, math
from collections import Counter

def need(ok, msg):
    if not ok: raise RuntimeError(msg)

def remainder(f,g,p):
    f=list(f)
    while len(f)>=len(g):
        a=f[-1]
        for j,c in enumerate(g): f[len(f)-len(g)+j]=(f[len(f)-len(g)+j]-a*c)%p
        while f and not f[-1]: f.pop()
    return f

def irreducible(f,p):
    return all(remainder(f,tuple(c)+(1,),p) for d in range(1,(len(f)-1)//2+1) for c in itertools.product(range(p),repeat=d))

def mobius(n):
    factors=0;d=2
    while d*d<=n:
        if n%d==0:
            n//=d;factors+=1
            if n%d==0:return 0
            while n%d==0:n//=d
        d+=1
    if n>1:factors+=1
    return (-1)**factors

def predicted_counts(p,e):
    divs=[d for d in range(1,e+1) if e%d==0]
    G=sum(mobius(e//d)*p**d for d in divs)
    if e%3: return G,0,G//(6*e)
    u=e//3;D=1 if p==3 else (2 if p%3==1 else 0)
    L=sum(mobius(e//d)*(D if d//math.gcd(u,d)==1 else p**math.gcd(u,d)+1) for d in divs)
    need((G-2*L)%(6*e)==0,'orbit integrality')
    return G,2*L,(G-2*L)//(6*e)

class Field:
    def __init__(self,p,e):
        self.p,self.e,self.q=p,e,p**e
        self.mod=next(tuple(c)+(1,) for c in itertools.product(range(p),repeat=e) if c[0] and irreducible(tuple(c)+(1,),p))
        self.digits=[tuple((x//p**i)%p for i in range(e)) for x in range(self.q)]
        self.add=[[self.encode([(a+b)%p for a,b in zip(u,v)]) for v in self.digits] for u in self.digits]
        self.neg=[self.encode([(-a)%p for a in u]) for u in self.digits]
        self.mul=[]
        for u in self.digits:
            row=[]
            for v in self.digits:
                c=[0]*(2*e-1)
                for i,a in enumerate(u):
                    for j,b in enumerate(v): c[i+j]=(c[i+j]+a*b)%p
                c=remainder(c,self.mod,p)
                row.append(self.encode(c))
            self.mul.append(row)
    def encode(self,c): return sum(a*self.p**i for i,a in enumerate(c))
    def sub(self,a,b): return self.add[a][self.neg[b]]
    def power(self,a,k):
        r=1
        while k:
            if k&1:r=self.mul[r][a]
            a=self.mul[a][a];k//=2
        return r
    def inv(self,a):
        need(a!=0,'zero inversion')
        return self.power(a,self.q-2)
    def degree(self,a):return next(i for i in range(1,self.e+1) if self.power(a,self.p**i)==a)
    def normalized(self,t):
        s=(0,1,t);out=[]
        for j in range(self.e):
            tj=self.power(t,self.p**j)
            for b,c in itertools.permutations(s,2):
                a=self.sub(c,b)
                if self.add[self.mul[a][tj]][b] in set(s)-{b,c}:out.append((a,b,j))
        return out
    def direct(self,t):
        s={0,1,t};out=[]
        for j in range(self.e):
            sj=[self.power(x,self.p**j) for x in s]
            for a in range(1,self.q):
                for b in range(self.q):
                    if {self.add[self.mul[a][x]][b] for x in sj}==s:out.append((a,b,j))
        return out

def run(p,e):
    F=Field(p,e);q=F.q
    for x in range(q):
        need(F.power(x,q)==x,'Frobenius q')
        if x:need(F.mul[x][F.inv(x)]==1,'inverse')
    # Frobenius multiplicativity/additivity and quadratic-character transport.
    fp=[F.power(x,p) for x in range(q)]
    for x in range(q):
        for y in range(q):
            need(fp[F.add[x][y]]==F.add[fp[x]][fp[y]],'Frobenius addition')
            need(fp[F.mul[x][y]]==F.mul[fp[x]][fp[y]],'Frobenius product')
    full=[t for t in range(q) if F.degree(t)==e];counts=Counter();bad=[];good=[]
    for t in full:
        stab=F.normalized(t);counts[len(stab)]+=1
        need(sum(j==0 for a,b,j in stab)==1,'affine kernel nontrivial')
        need(len(stab) in (1,3),'unexpected order')
        predicted=False
        if e%3==0:
            tr=F.power(t,p**(e//3))
            predicted=tr in (F.inv(F.sub(1,t)),F.mul[F.sub(t,1)][F.inv(t)])
        need((len(stab)>1)==predicted,'obstruction equivalence')
        (bad if predicted else good).append(t)
    need(good,'no regular triple')
    G,B,O=predicted_counts(p,e)
    need((G,B)==(len(full),len(bad)),'Mobius count')
    witness=good[0]
    need(set(F.direct(witness))==set(F.normalized(witness))=={(1,0,0)},'independent exhaustive witness')
    hostile=bad[0] if bad else 2
    need(set(F.direct(hostile))==set(F.normalized(hostile)),'independent exhaustive hostile')
    # Prime-subfield three-set is pointwise fixed by nonidentity Frobenius.
    need(all(F.power(x,p)==x for x in (0,1,2)),'prime-subfield hostile')
    need(fp!=list(range(q)),'Frobenius must be nonidentity')
    return dict(p=p,e=e,q=q,free_signed_triple_orbits=O,modulus_low_to_high=F.mod,full_degree_count=len(full),full_degree_stabilizer_counts=dict(sorted(counts.items())),regular_t=witness,regular_coefficients=F.digits[witness],hostile_t=hostile,hostile_coefficients=F.digits[hostile],hostile_stabilizer=F.normalized(hostile),bad_full_degree_parameters=bad,prime_subfield_hostile=[0,1,2],signed_group_order=q*(q-1)*e)

if __name__=='__main__': print(json.dumps([run(3,3),run(3,5),run(7,3)],indent=2,sort_keys=True))
