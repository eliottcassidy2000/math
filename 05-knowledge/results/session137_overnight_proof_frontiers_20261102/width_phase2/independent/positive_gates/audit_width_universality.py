"""Independent positive width-three controller and literal profile audit.

No producer imports. Small complete recipes are replayed as full polynomial
functions; larger CRT cases use independently derived controlled matrices.
"""
from itertools import product
from math import gcd,lcm
from collections import deque
from hashlib import sha256
import json
from tournament_kernel import matrix,require
from audit_width_states import profile

def factors(n):
    out=[];p=2
    while p*p<=n:
        if n%p==0:
            q=1
            while n%p==0:n//=p;q*=p
            out.append((p,q))
        p+=1
    if n>1:out.append((n,n))
    return out

def phi(n):
    out=n
    for p,_ in factors(n):out=out//p*(p-1)
    return out

def order(x,m):
    x%=m;y=1%m;k=0
    while True:
        y=y*x%m;k+=1
        if y==1%m:return k

def mm(a,b,m):
    return ((a[0]*b[0]+a[1]*b[2])%m,(a[0]*b[1]+a[1]*b[3])%m,
            (a[2]*b[0]+a[3]*b[2])%m,(a[2]*b[1]+a[3]*b[3])%m)

def mp(a,k,m):
    b=(1%m,0,0,1%m)
    while k:
        if k&1:b=mm(b,a,m)
        a=mm(a,a,m);k//=2
    return b

def matorder(a,m):
    eye=(1%m,0,0,1%m);b=eye;k=0
    while True:
        b=mm(b,a,m);k+=1
        if b==eye:return k

def transitive(n):
    return tuple(tuple(int(i<j) for j in range(n)) for i in range(n))

def combine(a,b,cyclic):
    n=len(a);m=len(b);N=n+m+int(cyclic);c=[[0]*N for _ in range(N)]
    for i in range(n):
        for j in range(n):c[i][j]=a[i][j]
    for i in range(m):
        for j in range(m):c[n+i][n+j]=b[i][j]
    for i in range(n):
        for j in range(m):c[i][n+j]=1
    if cyclic:
        for i in range(m):c[n+i][-1]=1
        for i in range(n):c[-1][i]=1
    return tuple(map(tuple,c))

def evaluate(p,x):
    ans=0;fall=1
    for k,value in enumerate(p,1):
        fall*=x-k+1;ans+=value*fall
    return ans

def unit_words(m):
    bank={1%m:(0,0,0)};queue=deque(bank)
    while queue:
        x=queue.popleft()
        for i,g in enumerate((5,17,23)):
            y=x*g%m
            if y not in bank:
                v=list(bank[x]);v[i]+=1;bank[y]=tuple(v);queue.append(y)
    return bank

class Controller:
    def __init__(self,M,full):
        self.M=M;self.full=full
        n=M;a=b=0
        while n%2==0:n//=2;a+=1
        while n%3==0:n//=3;b+=1
        self.two=2**a;self.three=3**b;self.n=n;self.a=a;self.b=b
        N=lcm(phi(self.three*n),phi(n))
        while N<max(a,b)+4:N+=lcm(phi(self.three*n),phi(n))
        self.N=N;self.r=self.two*self.three;self.tau=15*self.r*pow(4,-1,n)%n
        self.n0=n
        while self.n0%5==0:self.n0//=5
        i2=pow(2,-1,self.three);i4=pow(4,-1,n);i8=pow(8,-1,n)
        self.maps={'D1':(3,2,2,(2,0,0,1)),
                   'A':(9,4,4,(5,-1,1,0)),
                   'A4':(27,8,8,(11,-3,3,-2)),
                   'V':(-3,4,-i2,(1,0,0,i4))}
        self.orders={g:lcm(order(x,self.two),order(y,self.three),order(z,self.three),matorder(mat,n))
                     for g,(x,y,z,mat) in self.maps.items() if g!='A4'}
        self.primitive_calls=0;self.last=None;self.S_calls=0;self.inverse_calls=0
        if full:
            self.U0=[pow(x,N,M) for x in range(M)]
            self.U2=self.gate_fun('D1',self.U0)
            self.U4=self.gate_fun('D1',self.U2)
            self.V=self.pad(self.B(self.U2,[x*x%M for x in range(M)]),N-3)
    def pad(self,f,k):return [v*pow(x,k,self.M)%self.M for x,v in enumerate(f)]
    def B(self,f,g):
        M=self.M
        return [x*(f[x]*g[x]+(f[(x+1)%M]-f[x])*(g[(x+1)%M]-g[x])
                   +(f[(x-1)%M]-f[x])*(g[(x-1)%M]-g[x]))%M for x in range(M)]
    def gate_fun(self,g,f):
        if g=='V':return [x*y%self.M for x,y in zip(f,self.V)]
        aux={'D1':'U0','A':'U2','A4':'U4'}[g]
        return self.pad(self.B(f,getattr(self,aux)),self.N-1)
    def control(self,f):
        return (f[1%self.M]%self.two,f[1%self.M]%self.three,f[2%self.M]%self.three,
                f[1%self.M]%self.n,f[2%self.M]%self.n)
    def apply(self,g,k=1):
        x,y,z,mat=self.maps[g];h2,h3,g3,h,gg=self.state
        mat=mp(mat,k,self.n)
        self.state=(h2*pow(x,k,self.two)%self.two,h3*pow(y,k,self.three)%self.three,
                    g3*pow(z,k,self.three)%self.three,(mat[0]*h+mat[1]*gg)%self.n,
                    (mat[2]*h+mat[3]*gg)%self.n)
        if self.full:
            for _ in range(k):self.f=self.gate_fun(g,self.f)
            require(self.control(self.f)==self.state,'full polynomial gate versus controlled response')
            require(self.f[2%self.M]%self.two==0 and self.f[3%self.M]%self.three==0,'reset invariant')
        self.primitive_calls+=k
        if k:self.last=g
    def inverse(self,g,k=1):
        self.inverse_calls+=1;self.apply(g,(-k)%self.orders[g])
    def R(self):self.apply('A');self.inverse('D1',2)
    def E(self):self.apply('A4');self.inverse('D1',3)
    def S(self,k):
        old=self.state;self.S_calls+=1
        if self.full:
            for _ in range(self.r*k):
                self.apply('V');self.inverse('A');self.apply('D1',2);self.apply('A')
        else:
            h2,h3,g3,h,g=self.state
            self.state=(h2,h3,g3,(h+self.tau*k*g)%self.n,g)
            if k:self.last='A'
        require(self.state==old[:3]+((old[3]+self.tau*k*old[4])%self.n,old[4]),'positive shear realization')
    def seed(self,z):
        v=0
        while v<self.b and z%(3**(v+1))==0:v+=1
        modulus=self.two*3**(self.b-v)
        unit=next(u for u in range(modulus) if u%self.two==(z*pow(3**v,-1,self.two))%self.two
                  and u%(3**(self.b-v))==(z//(3**v))%(3**(self.b-v)))
        k,l,e=unit_words(modulus)[unit]
        if not k:k+=order(5,modulus)
        h=pow(3,v,self.M)*pow(5,k,self.M)*pow(17,l,self.M)*pow(23,e,self.M)%self.M
        g=pow(12,v,self.M)*pow(32,k,self.M)*pow(224,l,self.M)*pow(296,e,self.M)%self.M
        self.state=(h%self.two,h%self.three,g%self.three,h%self.n,g%self.n)
        require(self.state[0]==z%self.two and self.state[1]==z%self.three,'anchor seed')
        require(gcd(gcd(h,g),self.n)==1,'primitive odd seed')
        if self.n%5==0:require(h%5==0 and g%5!=0,'forced exceptional-prime sidecar')
        if self.full:
            X=list(range(self.M));seeds=[self.B(X,X),self.B([x*x%self.M for x in X],X),
                                        self.B([pow(x,4,self.M) for x in X],X),
                                        self.B([pow(x,3,self.M) for x in X],[x*x%self.M for x in X])]
            self.f=self.pad([pow(seeds[0][x],v,self.M)*pow(seeds[1][x],k,self.M)*pow(seeds[2][x],l,self.M)*pow(seeds[3][x],e,self.M)%self.M for x in X],self.N)
            require(self.control(self.f)==self.state,'full polynomial seed')
    def land(self,z):
        self.seed(z);h,g=self.state[3:]
        k=next(k for k in range(self.n0) if gcd(h+self.tau*k*g,self.n0)==1)
        self.S(k);self.R();require(gcd(self.state[4],self.n0)==1,'initial unit sidecar')
        if self.n%5==0:
            start=tuple(x%5 for x in self.state[3:]);queue=deque([start]);words={start:''};target=None
            while queue:
                p=queue.popleft()
                if p[0]==z%5 and p[1]:target=p;break
                for letter,q in [('R',(p[1],p[0])),('E',((2*p[0]+4*p[1])%5,(3*p[0]+3*p[1])%5))]:
                    if q not in words:words[q]=words[p]+letter;queue.append(q)
            require(target is not None,'exceptional-prime word exists')
            for letter in words[target]:
                h,g=self.state[3:];goal=0 if letter=='E' else g
                k=next(k for k in range(self.n0) if (h+self.tau*k*g-goal)%self.n0==0)
                self.S(k)
                if letter=='E':self.E()
                else:self.R()
                require(gcd(self.state[4],self.n0)==1,'safe E13 and R sidecars')
        require(gcd(self.state[4],self.n)==1,'final unit sidecar')
        h,g=self.state[3:];k=next(k for k in range(self.n) if (h+self.tau*k*g-z)%self.n==0)
        self.S(k+self.n)
        require(self.state[0]==z%self.two and self.state[1]==z%self.three and self.state[3]==z%self.n,'simultaneous CRT landing')
        require(self.last=='A','positive final cyclic context for strong extraction')
        if self.full:require(self.f[1%self.M]==z%self.M,'whole positive polynomial recipe')

def run():
    foundation=0;bank=[matrix(n,c) for n in range(1,4) for c in range(1<<(n*(n-1)//2))]
    profiles={a:profile(a) for a in bank}
    for a,b in product(bank,repeat=2):
        p,q=profiles[a],profiles[b];j=profile(combine(a,b,False));r=profile(combine(a,b,True))
        for x in range(-2,7):
            F=lambda t:evaluate(p,t);G=lambda t:evaluate(q,t)
            require(evaluate(j,x)==F(x)*G(x),'literal ordinal multiplication')
            require(evaluate(r,x)==x*(F(x)*G(x)+(F(x+1)-F(x))*(G(x+1)-G(x))+(F(x-1)-F(x))*(G(x-1)-G(x))),'literal bilinear cyclic identity')
        foundation+=1
    for r,s,expected in [(1,1,(3,12)),(2,1,(5,32)),(4,1,(17,224)),(3,2,(23,296))]:
        p=profile(combine(transitive(r),transitive(s),True))
        require((evaluate(p,1),evaluate(p,2))==expected,'literal seed coordinates')
    require(len(unit_words(72))==24,'full unit base modulo72')
    require((2*1+4*1)%5==1 and (3*1+3*1)%5==1,'fixed diagonal hostile')
    require(((11*2-3*3)*pow(8,-1,13))%13==0 and (3*2-2*3)%13==0,'E13 kernel hostile')
    cases=full_cases=primitive=inverses=0;digest=sha256()
    extra_full={(35,0),(35,13),(65,0),(65,13),(120,1),(120,119)}
    for M in list(range(1,97))+[120,125,169,175,343,455,1001,2340]:
        targets=range(M) if M<=96 else sorted({0,1,2,3,5,13,19,23,31,M//2,M-1})
        for z in targets:
            if (z-1)%gcd(2,M):continue
            C=Controller(M,M<=18 or (M,z) in extra_full);C.land(z)
            cases+=1;full_cases+=C.full
            if C.full:primitive+=C.primitive_calls;inverses+=C.inverse_calls
            digest.update(str((M,z,C.N,C.orders,C.state)).encode())
    print(json.dumps({'scope':'independent literal set-partition profiles and positive modular controller;all compatible targets M1..96 plus selected mixed/prime-power moduli;full polynomial recipes M1..18 plus six mixed cases at35,65,120',
                      'literal_factor_pairs':foundation,'literal_seed_checks':4,'simultaneous_landings':cases,
                      'full_polynomial_recipes':full_cases,'primitive_gate_applications_in_full_replays':primitive,'positive_inverse_calls_in_full_replays':inverses,
                      'semantic_sha256':digest.hexdigest(),
                      'verdict':'PASS bilinear foundation, anchors, positive inverses, resets, mod5 and13 hostiles, simultaneous landing and final strong-extraction marker'},sort_keys=True,indent=2))

if __name__=='__main__':run()
