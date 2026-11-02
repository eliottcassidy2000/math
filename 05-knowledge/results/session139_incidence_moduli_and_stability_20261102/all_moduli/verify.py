#!/usr/bin/env python3
"""Session139: exact contextual clocks, prime-function states, actual gates.
Standard library. No inherited program imports and no graph atlas inputs.
"""
from collections import deque
from array import array
from itertools import combinations
from math import factorial, gcd, lcm
import hashlib, json


def require(test, message):
    if not test:
        raise RuntimeError(message)


def factor(n):
    out = []
    p = 2
    while p*p <= n:
        if n % p == 0:
            a = 0
            while n % p == 0:
                n //= p; a += 1
            out.append((p,a))
        p += 1
    if n > 1: out.append((n,1))
    return out


def clock(M):
    if M == 1: return 1,1
    fs = factor(M); A=max(a for p,a in fs); L=1
    for p,a in fs:
        if p != 2: ell=(p-1)*p**(a-1)
        elif a <= 2: ell=1
        else: ell=2**max(1,a-3)
        L=lcm(L,ell)
    return A,L


def cutoff(M):
    c=1;d=1
    while d % M:
        c += 1;d *= c*(c-1)
    return c


def transitive_states(M, count):
    k=cutoff(M); s=[0]*k;s[0]=1
    ds=[0]+[factorial(c)*factorial(c-1)%M for c in range(1,k)]
    out=[]
    for n in range(1,count+1):
        s=[0]+[(c*s[c]+s[c-1])%M for c in range(1,k)]
        out.append(tuple(ds[c]*s[c]%M for c in range(1,k)))
    return out


def check_clocks():
    records=[];checks=0
    for M in list(range(1,513))+[625,1024,1331,2048,3125,6561,14641]:
        A,L=clock(M);ss=transitive_states(M,A+2*L)
        seen={};first=None
        for n,s in enumerate(ss,1):
            if s in seen:
                first=(seen[s],n-seen[s]);break
            seen[s]=n
        require(first==(A,L),('clock',M,first,(A,L)))
        for n in range(A-1,len(ss)-L):
            require(ss[n]==ss[n+L],('period',M,n));checks+=1
        records.append((M,A,L))
    # Wrong dyadic Euler-period and missing-preperiod controls.
    require(clock(16)==(4,2),'dyadic special period')
    require(transitive_states(25,3)[0]!=transitive_states(25,22)[20], 'preperiod hostile L1/L21')
    return {'moduli':len(records),'period_checks':checks,'selected':[r for r in records if r[0] in [1,8,16,25,55,121,625,1331,14641]],'digest':hashlib.sha256(json.dumps(records).encode()).hexdigest()}


def gate(v,par,p,op):
    d=len(v)
    def f(x):
        x %= p
        return 0 if x==0 else v[x-1] if x<=d else par*v[p-x-1]
    if op=='P': return tuple(x*f(x)%p for x in range(1,d+1)),-par
    if op=='C': return tuple((x*x*f(x)+x*(f(x+1)-f(x-1)))%p for x in range(1,d+1)),par
    if op=='D': return tuple(((x**3-2*x)*f(x)+(2*x*x+x)*f(x+1)+(x-2*x*x)*f(x-1))%p for x in range(1,d+1)),-par
    raise ValueError(op)


def determinant_mod(mat,p):
    a=[r[:] for r in mat]; ans=1
    for i in range(len(a)):
        j=next((j for j in range(i,len(a)) if a[j][i]%p),None)
        if j is None:return 0
        if j!=i:a[i],a[j]=a[j],a[i];ans=-ans
        z=a[i][i]%p;ans=ans*z%p;iz=pow(z,-1,p)
        for j in range(i+1,len(a)):
            t=a[j][i]*iz%p
            for k in range(i,len(a)):a[j][k]=(a[j][k]-t*a[i][k])%p
    return ans%p


def orbit(p):
    d=(p-1)//2;N=p**d;size=2*N
    def code(v,par):
        n=0
        for z in reversed(v):n=n*p+z
        return n+(N if par==-1 else 0)
    def decode(n):
        par=-1 if n>=N else 1;n%=N;v=[]
        for _ in range(d):v.append(n%p);n//=p
        return tuple(v),par
    start=code(tuple(range(1,d+1)),-1)
    parent=array('i',[-2])*size;ops=bytearray(size);parent[start]=-1
    q=array('i',[start]);i=0;dets={}
    for par in [1,-1]:
        for op in 'PCD':
            columns=[gate(tuple(int(j==k) for j in range(d)),par,p,op)[0] for k in range(d)]
            dets[f'{op}:{par}']=determinant_mod([[columns[j][i] for j in range(d)] for i in range(d)],p)
    while i<len(q):
        n=q[i];i+=1;v,par=decode(n)
        for j,op in enumerate('PCD'):
            w,ep=gate(v,par,p,op);m=code(w,ep)
            if parent[m]==-2:parent[m]=n;ops[m]=j;q.append(m)
    def word(v,par):
        n=code(v,par);require(parent[n]!=-2,('missing',p,v,par));w=[]
        while parent[n]>=0:w.append('PCD'[ops[n]]);n=parent[n]
        return ''.join(reversed(w))
    target_count=2*N if p in [3,7] else 2*N-2
    require(len(q)==target_count,('orbit_count',p,len(q),target_count))
    if p in [5,11]:
        require(all(dets.values()),('singular',p))
        require(parent[0]==parent[N]==-2,('unexpected unary zero',p))
        e1=(1,)+(0,)*(d-1);e2=(0,1)+(0,)*(d-2)
        wa,wb=word(e1,-1),word(e2,-1)
        for w,target in [(wa,e1),(wb,e2)]:
            poly=[0,1]
            for op in w:poly=poly_gate(poly,op)
            require(tuple(evaluate(poly,x)%p for x in range(1,d+1))==target,('full polynomial replay',p,w))
        zero={'left':wa,'right':wb,'left_order':1+sum('PCD'.index(x)+1 for x in wa),'right_order':1+sum('PCD'.index(x)+1 for x in wb),'operation':'J','profile_function':'zero'}
    else:zero={'word':word((0,)*d,-1)}
    return {'p':p,'tagged_unary_states':len(q),'contextual_classes':2*N-1,'determinants':dets,'zero_witness':zero,'visited_sha256':hashlib.sha256(parent.tobytes()).hexdigest()}


def add(a,b):
    c=[0]*max(len(a),len(b))
    for i,v in enumerate(a):c[i]+=v
    for i,v in enumerate(b):c[i]+=v
    return c

def mul(a,b):
    c=[0]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):c[i+j]+=x*y
    return c

def scale(a,n):return [n*x for x in a]

def shift(a,t):
    from math import comb
    c=[0]*len(a)
    for i,v in enumerate(a):
        for j in range(i+1):c[j]+=v*comb(i,j)*t**(i-j)
    return c

def evaluate(a,x):
    z=0
    for c in reversed(a):z=z*x+c
    return z

def poly_gate(a,op):
    if op=='P':return [0]+a
    plus,minus=shift(a,1),shift(a,-1)
    if op=='C':return add([0,0]+a,[0]+add(plus,scale(minus,-1)))
    return add(add(mul([0,-2,0,1],a),mul([0,1,2],plus)),mul([0,1,-2],minus))


def profiles(adj):
    n=len(adj);N=1<<n
    dp=[[0]*n for _ in range(N)];h=[0]*N;h[0]=1
    for v in range(n):dp[1<<v][v]=1
    for mask in range(1,N):
        for v in range(n):
            z=dp[mask][v]
            if not z:continue
            bits=adj[v]&~mask
            while bits:
                bit=bits&-bits;bits-=bit;w=bit.bit_length()-1;dp[mask|bit][w]+=z
        h[mask]=sum(dp[mask])
    pc=[[0]*(n+1) for _ in range(N)];pc[0][0]=1
    for mask in range(1,N):
        root=mask&-mask;part=mask
        while part:
            if part&root:
                for c in range(1,n+1):pc[mask][c]+=h[part]*pc[mask^part][c-1]
            part=(part-1)&mask
    return pc[-1]


def falling_eval(prof,x):
    z=0;t=1
    for c in range(1,len(prof)):
        t*=x-c+1;z+=prof[c]*t
    return z


def graph_gate(a,op):
    n=len(a)
    if op=='P':return [((1<<(n+1))-1)^1]+[x<<1 for x in a]
    # Cyclic blocks A, L_k, singleton; A -> L_k -> z -> A.
    k=1 if op=='C' else 2;z=n+k;out=a[:]+[0]*(k+1)
    for i in range(n):out[i]|=((1<<k)-1)<<n
    for i in range(k):out[n+i]|=(1<<z)|sum(1<<(n+j) for j in range(i+1,k))
    out[z]=(1<<n)-1
    return out


def actual_checks():
    classes={p:{} for p in [3,5,7,11]};function_classes={M:{} for M in range(1,33)};count=0;gate_count=0
    for n in range(1,6):
        edges=list(combinations(range(n),2))
        for mask in range(1<<len(edges)):
            a=[0]*n
            for k,(i,j) in enumerate(edges):
                if mask>>k&1:a[i]|=1<<j
                else:a[j]|=1<<i
            pc=profiles(a);count+=1
            for M,fc in function_classes.items():
                k=cutoff(M)
                sig=tuple(factorial(c)*factorial(c-1)*(pc[c] if c<len(pc) else 0)%M for c in range(1,k))
                transformed=[0]+[factorial(c-1)*pc[c] for c in range(1,len(pc))]
                values=tuple(falling_eval(transformed,x)%M for x in range(M))
                old=fc.get(values)
                require(old is None or old==sig,('Psi false merge',M))
                fc[values]=sig
            for p in classes:
                sig=tuple(factorial(c)*factorial(c-1)*pc[c]%p for c in range(1,min(p,len(pc))))+(0,)*max(0,p-len(pc))
                values=tuple(falling_eval(pc,x)%p for x in range(p))
                old=classes[p].get(values)
                require(old is None or old==sig,('prime function false merge',p))
                classes[p][values]=sig
                require(all(values[-x%p]==(-1)**n*values[x]%p for x in range(p)),('parity',p,n))
            if n<=4:
                # Literal subset/partition DP checks against inherited operator forms.
                poly=[0]*(n+1)
                fall=[1]
                for c in range(1,n+1):
                    fall=mul(fall,[-(c-1),1]);poly=add(poly,scale(fall,pc[c]))
                for op in 'PCD':
                    out=profiles(graph_gate(a,op));pred=poly_gate(poly,op)
                    require(all(evaluate(pred,x)==falling_eval(out,x) for x in range(len(pred)+1)),('actual gate',n,mask,op));gate_count+=1
    for p,m in classes.items():require(len(set(m.values()))==len(m),('false split',p))
    for M,m in function_classes.items():require(len(set(m.values()))==len(m),('Psi false split',M))
    return {'labelled_tournaments':count,'actual_gate_identities':gate_count,'prime_moduli':[3,5,7,11],'factorial_function_moduli':32,'factorial_function_tournament_checks':32*count}


def arithmetic_bridge():
    W5=[0,0,4,0,1]
    W11=[0,2,0,8,0,1]
    require(poly_gate(poly_gate([0,1],'P'),'C')==W5,'W5 polynomial')
    require(poly_gate(poly_gate([0,1],'P'),'D')==W11,'W11 polynomial')
    require([evaluate(W5,x)%5 for x in range(5)]==[0,0,2,2,0],'W5 supports')
    require([evaluate(W11,x)%5 for x in range(5)]==[0,1,0,0,4],'W11 supports')
    product=mul(W5,W11)
    require(all(evaluate(product,x)%5==0 for x in range(5)),'H55 zero function')
    require(evaluate(product,1)==55,'H55 actual count')
    e1=mul(W11,W11);e2=mul(mul(W5,W5),mul(W5,W5))
    for x in range(5):
        a=evaluate(e1,x)%5;b=evaluate(e2,x)%5
        require(a*a%5==a and b*b%5==b and a*b%5==0, 'projector')
        require((a+b)%5==int(x!=0),'complementary projectors')
    require(transitive_states(16,6)[3]==transitive_states(16,6)[5], 'ordinary F hostile contextual equality')
    require(pow(3,4,16)!=pow(3,6,16),'ordinary F hostile distinction')
    return {'W5_order':4,'W11_order':5,'zero_join_order':9,'zero_join_H':55,'W5_F2':32,'W11_F2':100,'ordinary_function_hostile':{'M':16,'orders':[4,6],'at_x3':[1,9]}}


def main():
    result={'status':'PASS','clock':check_clocks(),'actual':actual_checks(),'arithmetic_bridge':arithmetic_bridge(),'orbits':[orbit(p) for p in [3,5,7,11]]}
    print(json.dumps(result,sort_keys=True,indent=2))
if __name__=='__main__':main()
