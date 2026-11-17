from itertools import product,permutations

# Each mask tags exactly the old vertices beaten by q.
BANK=(
((108,73,34,52,71,82,12),82,7,137,657,645,59,61,-53,-79),
((80,105,251,209,194,89,128,35),136,18,89,344,508,18,0,57,39),
((74,44,41,112,7,17,54),50,44,151,653,779,73,59,-71,-75),
((142,116,24,226,105,5,165,54),129,65,523,1835,1489,109,137,311,323),
((120,145,75,162,172,198,26,69),35,22,659,3150,2952,234,256,-93,-99),
((158,164,56,242,98,1,167,52),144,6,389,1168,1376,102,88,217,225),
((10,52,81,70,105,141,163,31),40,20,521,1662,1662,124,124,331,267),
)

def end_profiles(a):
    n=len(a);u=1<<n;dp=[[0]*n for _ in range(u)];h=[0]*u;h[0]=1
    for i in range(n):dp[1<<i][i]=1
    for m in range(1,u):
        h[m]=sum(dp[m])
        for i,v in enumerate(dp[m]):
            if not v:continue
            s=a[i]&~m
            while s:
                b=s&-s;s-=b;dp[m|b][b.bit_length()-1]+=v
    return h,dp

def profile(a,q):
    h,ends=end_profiles(a)
    _,starts=end_profiles(tuple(sum((a[j]>>i&1)<<j for j in range(len(a))) for i in range(len(a))))
    e=[sum(v for j,v in enumerate(row) if q>>j&1) for row in ends]
    s=[sum(v for j,v in enumerate(row) if q>>j&1) for row in starts]
    f=len(h)-1
    J=h[f]
    P=sum(h[x]*h[f^x] for x in range(f+1))
    W=sum(s[x]*h[f^x] for x in range(f+1))
    Z=sum(e[x]*h[f^x] for x in range(f+1))
    U=sum(s[x]*e[f^x] for x in range(f+1))
    return J,P,W,Z,e[f],s[f],J-U

def append(profile0,A):
    J,W,K,E=profile0
    h,p,w,z,e,s,k=A
    return J*h,p*W+J*w,h*K-J*(h-k)-E*s-W*z,J*e

def external(a,q,k):
    # R first, then q, then TT_k in reverse order of its labels.
    n=len(a);D=tuple(x|((1-(q>>i&1))<<n) for i,x in enumerate(a))+(q,)
    t=list(D)
    for j in range(k):
        outgoing=(1<<n)|(((1<<j)-1)<<(n+1))
        t=[x|((1-(outgoing>>i&1))<<len(t)) for i,x in enumerate(t)]+[outgoing]
    return t

def ordinal(a,b):
    return tuple(x|(((1<<len(b))-1)<<len(a)) for x in a)+tuple(x<<len(a) for x in b)

def strong(a):
    for i in range(len(a)):
        seen=1<<i
        for _ in a:
            for j in range(len(a)):
                if seen>>j&1:seen|=a[j]
        if seen!=(1<<len(a))-1:return False
    return True

def verify():
    for q,k in [(70,-35),(69,-21)]:
        base=(56,97,51,54,98,64,13)
        p=profile(base,q)
        if (p[0],p[2],p[-1])!=(101,553,k):raise ValueError(p)
        for t in (1,2,3):
            h=end_profiles(external(base,q,t))[0][-1]
            if h!=553*2**t+k or not strong(external(base,q,t)):raise ValueError(h)
        print("twin",q,p)
    A=[]
    for rows,q,q1,h,w,z,e,s,k,k1 in BANK:
        p,p1=profile(rows,q),profile(rows,q1)
        if p[:-1]!=p1[:-1] or (p[0],*p[2:])!=(h,w,z,e,s,k) or p1[-1]!=k1:raise ValueError((p,p1))
        print("local",p,p1[-1])
        A.append((p,p1))
    return A

def state(A,order,word):
    out=(1,0,1,0)
    for i,v in zip(order,word):out=append(out,A[i][v])
    return out

def crt(A,order,J,W,K0,H,k):
    r=H-W*2**k-K0
    bits=[]
    for i in order:
        h=A[i][0][0]
        dr=(A[i][1][-1]-A[i][0][-1])*(J//h)
        a=r%h
        if a==0:bits.append(0)
        elif a==dr%h:bits.append(1)
        else:return None
    if sum((A[i][1][-1]-A[i][0][-1])*(J//A[i][0][0])*v for i,v in zip(order,bits))!=r:return None
    return bits

def main():
    A=verify()
    chosen=None
    for n in (2,3,4,5,6,7):
        best=(0,None)
        for order in permutations(range(len(A)),n):
            J,W,K0,_=state(A,order,(0,)*n)
            if W%7:continue
            residues=[(K0+sum((A[i][1][-1]-A[i][0][-1])*(J//A[i][0][0])*b for i,b in zip(order,v)),v) for v in product(range(2),repeat=n)]
            seven=[(K,v) for K,v in residues if K%7==0]
            if len(seven)>best[0]:
                best=(len(seven),(order,J,W,K0,seven))
        if best[1]:
            print("best",n,best[0],best[1][0])
            chosen=best[1]
    order,J,W,K0,seven=chosen
    for K,v in seven:
        a=crt(A,order,J,W,K0,W*2**4+K,4)
        if a!=list(v):raise ValueError(a)
        rows=();q=0
        for i,b in zip(order,v):
            q|=BANK[i][1+b]<<len(rows);rows=ordinal(rows,BANK[i][0])
        if not strong(external(rows,q,1)):raise ValueError('not strong')
    if {K//7%16 for K,_ in seven}!=set(range(1,16,2)):raise ValueError('mod16')
    if {K//7%9 for K,_ in seven}!=set(range(9)):raise ValueError('mod9')
    print("mod16",sorted({K//7%16 for K,_ in seven}),"mod9",sorted({K//7%9 for K,_ in seven}))
    print("joint",order,J,W,K0,"seven",len(seven),"offset range",min(t[0] for t in seven),max(t[0] for t in seven))
    # Independent direct literal DP on two smaller tagged blocks combined.
    order=(0,1);u=(1,0);d,x=(),0
    for i,b in zip(order,u):
        rows=BANK[i][0];q=BANK[i][1+b]
        x|=q<<len(d);d=ordinal(d,rows)
    ref=profile(d,x)
    J,W,K,E=state(A,order,u)
    if (ref[0],ref[2],ref[-1],ref[4])!=(J,W,K,E):raise ValueError(ref)
    print("direct-two-block",J,W,K,E)

if __name__=='__main__':main()
