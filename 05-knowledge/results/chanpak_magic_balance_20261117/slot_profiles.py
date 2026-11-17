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

