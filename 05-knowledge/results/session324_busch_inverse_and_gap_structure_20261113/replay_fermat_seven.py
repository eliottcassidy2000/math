from replay_gap_structure import ends,check,strong,req,A,add

I=(4,1,10,67,15,31,55)
W=(8,1,3,6)
S=81;cut=13

def params(T,u):
    n=len(T);N=(1<<n)-1
    E=ends(T);P=ends(tuple(N^v^(1<<i) for i,v in enumerate(T)))
    h=[1]+[sum(r) for r in E[1:]]
    s=[sum(x for i,x in enumerate(t) if u>>i&1) for t in P]
    e=[sum(x for i,x in enumerate(t) if u>>i&1) for t in E]
    b=[sum(x for i,x in enumerate(t) if not(u>>i&1)) for t in E]
    alpha=sum(h[N^m]*s[m] for m in range(1,N+1))
    lam=sum(e[m]*s[N^m] for m in range(1,N))
    beta=h[N]+sum(b[m]*h[N^m] for m in range(1,N+1))
    d=sum(ends(add(T,u))[-1])
    return h[N],alpha,lam,beta,d

h,a,l,_,_=params(I,S);q,_,_,b,d=params(W,cut)
req((h,a,l,q,d,b)==(29,256,74,5,13,14))
req(strong(I) and strong(W) and a*b==3584 and h*d-q*l==7)

def tt(n):
    return tuple(((1<<n)-1)^((1<<(i+1))-1) for i in range(n))

def join(X,Y):
    return tuple(v|(((1<<len(Y))-1)<<len(X)) for v in X)+tuple(v<<len(X) for v in Y)

for m in (1,2):
    T=add(join(join(I,tt(m)),W),S|(cut<<(len(I)+m)))
    check(T,7*(2**(m+9)+1))

for s,T in enumerate([
    (54,64,90,83,2,94,17),A,
    (132,73,10,65,111,79,5,126),
    (150,176,194,215,164,77,147,32),
    (48,173,41,209,38,200,151,21),
    (440,273,483,486,460,146,291,322,32),
    (450,16,107,163,109,387,42,342,94),
    (282,680,699,32,682,897,447,777,542,73),
],2):check(T,7*(2**s+1))
print('strong seven-Fermat bridge s=2..9 and two literal rail lengths pass')
