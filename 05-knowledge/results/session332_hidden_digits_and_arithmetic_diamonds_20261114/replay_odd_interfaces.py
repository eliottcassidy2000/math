from collections import defaultdict
from itertools import permutations,product

def req(ok,s):
    if not ok:raise RuntimeError(s)

def g(m):
    req(m%3!=0,"not a unit")
    k=next(k for k in range(4) if ((1<<k)*m)%9 in (2,5))
    return k,((m<<k)+1)//3

def pred(u):
    n=3*u-1
    d=2 if u%3==1 else 4
    return [n>>k for k in range(d) if n%2**k==0]

def loop(m,r):
    for _ in range(r):
        a,m=g(m);b,m=g(m)
        if (a,b)!=(3,0):return False
    return True

def arithmetic():
    P=defaultdict(set)
    for m in range(1,5000):
        if m%3:P[g(m)[1]].add(m)
    for u in range(1,1300):
        if u%3:req(set(pred(u))==P[u],"preimage menu disagrees")
    req([2*x-1 for x in pred(4)]==[21],"wrong first hole predecessor")
    req([2*x-1 for x in pred(11)]==[63,31,15,7],"wrong second hole predecessors")
    for r in range(1,5):
        d=9**r
        for m in range(1,3*d+1):
            if m%3:
                req(loop(m,r)==((m-4)%d==0 and (pow(8,r-1,3)*((m-4)//d))%3 in (0,2)),
                    "frozen loop residue error")
    for a in range(5,18,2):
        M=5*(1+(1<<a))
        req([2*u-1 for u in pred((M+1)//2)]==
          [3*M,7+30*2**(a-2),3*(1+20*2**(a-4)),1+15*2**(a-3)],"menu arithmetic wrong")
    print("greedy affine conjugacy, exact predecessor menu, four repeat cylinders and physical-menu arithmetic passed")

def strong(A):
    q=1
    for _ in A:
        q|=sum(1<<i for i in range(len(A)) if any(q>>j&1 and A[j]>>i&1 for j in range(len(A))))
    if q!=(1<<len(A))-1:return False
    q=1
    for _ in A:
        q|=sum(1<<i for i in range(len(A)) if A[i]&q)
    return q==(1<<len(A))-1

def extend(A,c):
    return tuple(r|(((c>>i&1)^1)<<len(A)) for i,r in enumerate(A))+(c,)

def canon(A):
    G=defaultdict(list)
    for i,x in enumerate(A):G[x.bit_count()].append(i)
    return min(tuple(sum((A[i]>>j&1)<<k for k,j in enumerate(p)) for i in p)
       for blocks in product(*(tuple(permutations(G[k])) for k in sorted(G)))
       for p in [sum(blocks,())])

def visit(A,v):
    W=tuple(x for x in range(len(A)) if x!=v)
    F=[0]*(len(A)+1)
    for P in permutations(W):
        if not P:return (1,1)
        ok=int(bool(A[v]>>P[0]&1))+(int(bool(A[P[-1]]>>v&1))<<len(P))
        C=0
        for h,(x,y) in enumerate(zip(P,P[1:]),1):
            if not A[x]>>y&1:C|=1<<h
            if A[x]>>v&1 and A[v]>>y&1:ok|=1<<h
        if C&~ok:continue
        O=ok&~C
        t=O
        while True:
            F[(C|t).bit_count()]+=1
            if not t:break
            t=(t-1)&O
    return F

def hp(A):
    M=1<<len(A)
    E=[[0]*len(A) for _ in range(M)]
    for i in range(len(A)):E[1<<i][i]=1
    for s in range(1,M):
        for i,z in enumerate(E[s]):
            if z:
                b=A[i]&~s
                while b:
                    w=b&-b;b-=w
                    E[s|w][w.bit_length()-1]+=z
    return sum(E[-1])

def inflate(A,v,a):
    n=len(A)
    p=tuple(i for i in range(n) if i!=v)+(v,)*a
    return tuple(sum(1<<j for j,y in enumerate(p) if (x!=y and A[x]>>y&1) or (x==y==v and i<j))
                  for i,x in enumerate(p))

R=[
  ((2,12,9,17,7),0,1,6),
  ((2,12,17,37,43,71,31),0,5,32),
  ((2,4,17,7,43,15),0,5,5),
  ((2,8,3,37,15,87,31),0,5,10),
  ((2,12,25,97,11,23,55),0,7,30),
  ((2,4,9,51,39,71,31),6,1,20),
  ((2,4,41,19,7,91,31),3,1,15),
]
def physical():
    for A,v,C,D in R:
        F=visit(A,v)
        req(F[1]==C+2*D and F[2]==D and not any(F[3:]),"wrong visit profile")
        for a in range(1,5):
            X=inflate(A,v,a)
            req(strong(X) and hp(X)==C+D*(1<<a),"literal inflation failed")
    X={(0,)};O={1:set(),5:set(),7:set()}
    sizes=[1]
    for n in range(2,8):
        X={canon(extend(A,c)) for A in X for c in range(1<<(n-1))}
        sizes.append(len(X))
        for A in X:
            if strong(A):
                for v in range(n):
                    F=visit(A,v);D=F[2];C=F[1]-2*D
                    if D and C in O and not any(F[3:]):O[C].add(D//(D&-D))
    req(sizes==[1,1,2,4,12,56,456],"tournament census count")
    req(3 not in O[7],"offset 7 slope 3 physical through 7")
    full=[D for D in sorted(O[5]&O[1]) if D%3 and 3*D in O[7]&O[1]]
    req(full==[5],"unexpected new full branch ray")
    print("two-visit frames, literal strong inflation and complete order-7 atlas passed;",O)

def apex(w):
    A=tuple(sum(1<<j for j in range(i+1,len(w))) | ((w[i]=='1')<<len(w))
            for i in range(len(w)))
    return A+(sum(1<<i for i,b in enumerate(w) if b=='0'),)

def apex_formula(w):
    return (w[0]=='0')+(w[-1]=='1')+sum(a+b=='10' for a,b in zip(w,w[1:]))+sum(
        1<<(i-j-1) for i,b in enumerate(w) for j,a in enumerate(w[:i]) if a+b=='01')

def one_apex():
    for m in range(1,4):
        for q in range(2,5):
            for b in range(1,4):
                w='0'+'1'*m+'0'*b+'1'+'0'*(q-1)+'1'
                h=((1<<m)+1)*((1<<q)+1)*(1<<b)+(1<<m)-(1<<(q-1))+1
                req(h==apex_formula(w)==hp(apex(w)) and strong(apex(w)),"archetype formula")
    D=135
    output=[]
    for d in range(1,D+1,2):
        if D%d:continue
        U=''.join('0' if c=='1' else '1' for c in (bin(d-1)[2:] if d>1 else ''))
        V=bin(D//d)[2:][::-1]
        w=U+'0'+V
        req(hp(apex(w))==apex_formula(w),"literal 135-slope archetype")
        output.append(apex_formula(w)-2*D)
    req(sorted(output)==[-117,-27,-17,-5,5,17,35,85],"135 offset atlas")
    print("one-apex archetype, slope-45 return and all eight slope-135 intercepts passed")

if __name__=="__main__":
    arithmetic()
    physical()
    one_apex()
