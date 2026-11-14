from replay_ready_overlaps import end,converse,strong,profile,req

def join(A,B):
    return tuple(a|(((1<<len(B))-1)<<len(A)) for a in A)+tuple(b<<len(A) for b in B)

def Bmake(C,y):
    V=join((0,),C);n=len(V)
    return tuple(r|((i==1+y)<<n) for i,r in enumerate(V))+((((1<<n)-1)^(1<<(y+1))),)

def glue(A,v,B,w):
    n=len(A);p=[i for i in range(len(B)) if i!=w]
    N=n+len(p);R=[0]*N
    for i,r in enumerate(A):
        R[i]=r
        if i!=v:R[i]|=((1<<len(p))-1)<<n
        else:R[i]|=sum(1<<(n+j) for j,u in enumerate(p) if B[w]>>u&1)
    for j,u in enumerate(p):
        R[n+j]=sum(1<<(n+k) for k,t in enumerate(p) if B[u]>>t&1)
        if B[u]>>w&1:R[n+j]|=1<<v
    return tuple(R)

C=(10,28,25,16,1);letters=[Bmake(C,y) for y in (0,3)]
A=(8,1,3,6);H,h,L=profile(A,1)
for i,B in enumerate(letters):
    req(strong(B) and profile(B,6)==((31,55)[i],9,157))
    T=glue(A,1,B,6)
    req(strong(T) and profile(T,1)==(9*H+(22+24*i)*L,9*h,168*L-11*h))

for k in range(1,15):
    h0,L0=h,L;Ls=[]
    for _ in range(k):Ls.append(L0);h0,L0=9*h0,168*L0-11*h0
    w=[9**(k-1-i)*x for i,x in enumerate(Ls)]
    req(all(w[i]>sum(w[:i]) for i in range(k)))
print('literal binary seams and superincreasing inverse controls pass')
