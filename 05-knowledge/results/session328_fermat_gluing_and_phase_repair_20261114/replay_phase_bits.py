import itertools
from math import gcd,lcm,prod

def req(x):
    if not x: raise RuntimeError("replay mismatch")

def h(k,j): return 2048*4**k-1026*2**k+119+12*j

def q(k): return 354294*9**k-242340*6**k+50336*4**k-2052*2**k+2

def v2(x):
    return (abs(x)&-abs(x)).bit_length()-1

def physical(k,l,side):
    return h(k,side)*h(l,1-side),q(k)*q(l)

def safe(k): return k%10 not in (1,4)

def generators(D,e=4):
    M=2**e*D
    edges=[b*pow(a,-1,M)%M for a,b in ((73,85),(5,19),(13,41),(23,79))]
    # Successive cyclic quotients give one physical tuple at every base residue.
    reps={1:()}
    lengths=[]
    for g in edges:
        x=g; d=1
        while x not in reps:
            x=x*g%M; d+=1
        lengths.append(d-1)
        reps={x*pow(g,j,M)%M: coeff+(j,) for x,coeff in reps.items() for j in range(d)}
    H0=prod(a**n for a,n in zip((73,5,13,23),lengths))
    Q0=prod(q**n for q,n in zip((1088,256,416,1184),lengths))
    req(len(reps)==sum(gcd(i,M)==1 for i in range(M)))
    return M,H0,Q0,edges,reps

def bitplan(e,safe_only,shift,physical_mult):
    plan=[]
    k=0
    while k+3+shift<e:
        L=k+1
        while safe_only and not safe(L):
            L+=1
        a,b=(physical(k,L,j)[0]**physical_mult for j in (0,1))
        bits=min(L-k,e-(k+3+shift))
        copies=(1<<bits)-1
        plan.append((k+3+shift,bits,a,b,copies,(q(k)*q(L))**physical_mult))
        req(v2(a-b)==k+3+shift)
        k=L
    return plan

def dyadic_orbit(e,safe_only):
    M=2**e
    if e<3: return
    H={a*b%M for a in (73,85) for b in (373,379)}
    for i,bits,a,b,copies,Q in bitplan(e,safe_only,0,1):
        H={s*a**(copies-j)*b**j%M for s in H for j in range(copies+1)}
        req(len({x%(2**(i+bits)) for x in H})==2**(i+bits-1))
    req(len(H)==M//2)

def mixed_orbit(e,D,lam):
    M=2**e*D
    C,H0,Q0,gs,reps=generators(D)
    G0=2**4*D
    req(C==G0 and e>=4 and v2(lam)==1)
    bits=bitplan(e,True,1,lam)
    req(all((a-b)%D==0 for i,w,a,b,r,Q in bits))
    # The 4 controls are all literal common-clock switches; exhaust residues at
    # the small mixed base and sample lifts by exact descending digit targets.
    rules=[(73,85),(5,19),(13,41),(23,79)]
    wanted=[i for i in range(M) if gcd(i,M)==1]
    if len(wanted)>2500: wanted=wanted[::(1+len(wanted)//2500)]
    B=prod(pow(a,r,M) for i,w,a,b,r,Q in bits)%M
    for z in wanted:
        coeff=reps[z*pow((H0*B)%C,-1,C)%C]
        # Number of copies is the largest used exponent for a coset coordinate.
        n=[max(w[i] for w in reps.values()) for i in range(4)]
        x=prod(pow(a,t-j,M)*pow(b,j,M) for (a,b),t,j in zip(rules,n,coeff))%M
        x=x*B%M
        req(x%C==z%C)
        for i,w,a,b,r,Q in bits:
            d=2**(i+w)*D
            r0=b*pow(a,-1,M)%M
            t=next(t for t in range(r+1) if x*pow(r0,t,d)%d==z%d)
            x=x*pow(r0,t,M)%M
        req(x==z)

def slot_graph(A,z):
    n=len(A); N=1<<n
    end=[[0]*n for _ in range(N)]; H=[1]+[0]*(N-1)
    for S in range(1,N):
        for j in range(n):
            if S>>j&1:
                R=S^(1<<j)
                end[S][j]=1 if not R else sum(end[R][i] for i in range(n) if A[i]>>j&1)
        H[S]=sum(end[S])
    C=[1]+[0]*(N-1)
    for _ in range(z):
        D=[]
        for S in range(N):
            s=0; R=S
            while True:
                s+=H[R]*C[S^R]
                if not R:break
                R=(R-1)&S
            D.append(s)
        C=D
    return C

def R(A,v,z,sign):
    N=(1<<len(A))-1
    old=[j for j in range(len(A)) if j!=v]
    Q=slot_graph([sum(((A[i]>>j)&1)<<k for k,j in enumerate(old)) for i in old],z)
    L=1<<len(old); n=len(old)
    P=[[0]*n for _ in range(L)]
    for k,j in enumerate(old):
        P[1<<k][k]=(A[v]>>j)&1
    for S in range(1,L):
        for k,i in enumerate(old):
            if not P[S][k]: continue
            for u,j in enumerate(old):
                if not S>>u&1 and A[i]>>j&1:
                    P[S|(1<<u)][u]+=P[S][k]
    return Q[-1]+sum(sum(P[S])*sign**S.bit_count()*Q[(L-1)^S] for S in range(1,L))

def attach(A,v,j):
    n=len(A); N=n+j+1
    out=[0]*N
    for i in range(n):out[i]=A[i]|(((1<<(j+1))-1)<<n)
    out[v]^=1<<(N-1)
    for i in range(n,N-1):out[i]=((1<<N)-1)^((1<<(i+1))-1)
    out[-1]=1<<v
    return out

def cut_formula(w,z):
    n=len(w)
    S=z**(n+1)
    for i in range(n):
        for j in range(i+1,n):
            if w[i] and not w[j]:
                gap=j-i-1
                S+=z**(n-gap-1)*((z+1)**gap+(z-1)**gap)
    return S

def converse(A):
    return [((1<<len(A))-1)^r^(1<<i) for i,r in enumerate(A)]

def attach_general(A,B):
    # Root is last in each argument. Put merged root last in output.
    n,m=len(A)-1,len(B)-1
    full=(1<<m)-1
    out=[(A[i]&((1<<n)-1))+(full<<n)+(((A[i]>>n)&1)<<(n+m)) for i in range(n)]
    out += [((B[i]&full)<<n)+(((B[i]>>m)&1)<<(n+m)) for i in range(m)]
    out += [A[-1]+(B[-1]<<n)]
    return out

if __name__=='__main__':
    for A,H,Q in (((2,12,49,69,11,27,55),73,1088),
                  ((16,5,65,7,46,15,59),85,1088),
                  ((82,28,153,209,64,31,166,51),373,5896),
                  ((196,101,216,51,227,197,8,74),379,5896)):
        req(slot_graph(A,1)[-1]==H and slot_graph(A,2)[-1]==Q)
        for v in range(len(A)):
            S=1<<v
            for _ in A:
                for j,r in enumerate(A):
                    if S>>j&1:S|=r
            req(S==(1<<len(A))-1)
    for k in range(15):
        for j,base in enumerate(('10101010100','10111011100')):
            w=list(map(int,'1'*k+base+'0'*k))
            req(cut_formula(w,1)==h(k,j) and cut_formula(w,2)==q(k))
            req((gcd(h(k,j),99)==1)==(j or safe(k)))
            n=len(w)
            rp=2**n+sum(w[i]*2**i*3**(n-1-i) for i in range(n))
            rm=2**n-sum(w[i]*2**i for i in range(n))
            for L in range(1,9):
                out=w+[1]*L+[0]
                req(cut_formula(out,2)==2**(L+1)*q(k)+2*(3**L*rp-rm))
    for n in range(1,14):
        dyadic_orbit(n,False)
        dyadic_orbit(n,True)
    mixed_orbit(7,3,2)
    mixed_orbit(7,9,6)
    mixed_orbit(7,11,10)
    mixed_orbit(8,99,30)
    for A in [(2,4,1),(8,1,3,6),(8,5,9,50,7,23)]:
        for v in range(len(A)):
            for z in (1,2,3):
                a,b=R(A,v,z,1),R(A,v,z,-1)
                q0=slot_graph(A,z)[-1]
                for L in (1,2,3):
                    T=attach(A,v,L)
                    q1=slot_graph(T,z)[-1]
                    old=list(range(v))+list(range(v+1,len(A)))
                    C=[sum(((A[i]>>j)&1)<<k for k,j in enumerate(old)) for i in old]
                    qb=slot_graph(C,z)[-1]
                    if z==2 and L<=2:
                        req(R(T,v,z,1)==(z+1)**(L+1)*a-z**L*qb)
                        req(R(T,v,z,-1)==(z-1)**(L+1)*b+z**L*qb)
                    if A[v].bit_count()==1:
                        req(q0==z*(a+b-qb))
                    req(q1==z**(L+1)*q0+z*((z+1)**L*a-(z-1)**L*b))
    pair=((24,417,307,486,394,273,439,293,1),(36,273,314,259,361,394,303,351,1))
    for A,h,H,L in zip(pair,(525,561),(263,287),(568,607)):
        req(A[-1]==1)
        C=tuple(x&255 for x in A[:-1])
        req((slot_graph(A,1)[-1],slot_graph(C,1)[-1],R(A,8,1,1))==(h,H,L))
        req((slot_graph(A,2)[-1],R(A,8,2,1),R(A,8,2,-1),slot_graph(C,2)[-1])==(13152,7270,3922,4616))
        req(gcd(h+4*L,99)==1)
    second=((168,441,379,368,449,336,387,300,1),(82,316,353,469,484,329,258,359,1))
    orbit=[]
    for A,h,L in zip(second,(829,841),(877,892)):
        C=tuple(x&255 for x in A[:-1])
        req((slot_graph(A,1)[-1],slot_graph(C,1)[-1],R(A,8,1,1))==(h,349,L))
        req((slot_graph(A,2)[-1],R(A,8,2,1),R(A,8,2,-1),slot_graph(C,2)[-1])==(17732,9583,4771,5488))
        for u in range(len(A)):
            S=1<<u
            for _ in A:
                for i,row in enumerate(A):
                    if S>>i&1:S|=row
            req(S==(1<<len(A))-1)
        q,p,m,b=17732,9583,4771,5488
        for j in (1,1,16,2):
            h,L=h+(1<<j)*L,(2<<j)*L-349
            q,p,m,b=(2<<j)*q+2*(3**j*p-m),3**(j+1)*p-(1<<j)*b,m+(1<<j)*b,(2<<j)*b
        req(q==1909836390404740)
        orbit.append((h,L))
    req(tuple(h for h,L in orbit)==(7247174993,7388732915))
    seen=set()
    for t in range(990):
        if t%10==0:
            a,b=(h for h,L in orbit)
            req((b-a)%396==198)
            seen.add(a%99)
        orbit=[((h+4*L)%396,(8*L-349)%396) for h,L in orbit]
    req(len(seen)==99 and tuple((h,L) for h,L in orbit)==((377,71),(179,311)))
    from heapq import heappop,heappush
    for restricted,bound,answer,count0 in [
        (False,1090389448,((1,2,5),812755,828397),92),
        (True,4682656948,((1,4,2,1),2067209,2106611),80)]:
        heap=[(17732,(),829,877,841,892,9583,4771,5488)]
        count=0
        while heap:
            q,w,h,l,H,L,p,s,b=heappop(heap); count+=1
            if (H-h)%396==198 and gcd(H*h,99)==1:
                req(q==bound and count==count0 and (w,h,H)==answer)
                break
            for j in itertools.count(1):
                u=1<<j
                qq=2*u*q+2*(3**j*p-s)
                if qq>bound: break
                if restricted and j&(j-1): continue
                heappush(heap,(qq,w+(j,),h+u*l,2*u*l-349,H+u*L,2*u*L-349,
                               3**(j+1)*p-u*b,s+u*b,2*u*b))
        else:
            raise RuntimeError('missing initialized phase')
    for j in range(1,22):
        L=j+2
        a,b=525+568*2**j,561+607*2**j
        c,d=525+568*2**L,561+607*2**L
        req(v2(a*d-b*c)==j)
        req((gcd(a,99)==gcd(b,99)==1)==(j%10!=1))
    for A in ((2,4,1),(8,1,3,6)):
        for B in ((2,4,1),(16,9,3,5,14)):
            T=attach_general(A,B); C=tuple(x&((1<<(len(A)-1))-1) for x in A[:-1])
            D=tuple(x&((1<<(len(B)-1))-1) for x in B[:-1])
            for z in (1,2,3):
                qa,qb,qc,qd=[slot_graph(S,z)[-1] for S in (A,B,C,D)]
                qq=qa*qd+qc*qb-z*qc*qd
                for e in (-1,1):
                    ra=R(A,len(A)-1,z,e)
                    rb=R(B,len(B)-1,z,e)
                    eb=R(converse(B),len(B)-1,z,e)
                    req(R(T,len(T)-1,z,e)==(ra-qc)*slot_graph(D,z+e)[-1]+qc*rb)
                    qq+=z*(ra-qc)*(eb-qd)
                req(slot_graph(T,z)[-1]==qq)
    print("same-Q bit controls, safe 99 orbit, and Fermat boundary transfer verified")
