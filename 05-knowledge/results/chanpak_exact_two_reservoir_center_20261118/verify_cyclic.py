from math import isqrt

A=(2,4,9,19,39,15)
B=(0,1,3,23,39,15)
C=(128,5,9,19,39,79,31,126)
D=(4,9,26,81,163,15,55,111)
M79=(2,20,17,7,105,15,175,63)
M85=(16,5,65,7,46,15,59,127)
M91=(8,5,17,22,67,31,175,63)

def TT(n):return tuple((1<<n)-(1<<(i+1)) for i in range(n))
def ord(A,B):return tuple(x|(((1<<len(B))-1)<<len(A)) for x in A)+tuple(x<<len(A) for x in B)
def host(A,B,m=0,n=0):
    A=ord(A,TT(m));B=ord(B,TT(n));a,b=len(A),len(B)
    return tuple(x|(((1<<b)-1)<<a) for x in A)+tuple((x<<a)|(1<<(a+b)) for x in B)+((1<<a)-1,)
def counts(rows):
    n=len(rows);A=[[0]*n for _ in range(1<<n)]
    for i in range(n):A[1<<i][i]=1
    for s in range(1,len(A)):
        for i,w in enumerate(A[s]):
            mask=rows[i]&~s
            while mask and w:
                j=mask&-mask;mask-=j;A[s|j][j.bit_length()-1]+=w
    return [1]+[sum(q) for q in A[1:]]
def profile(rows):
    H=counts(rows);M=len(H)-1
    return H[-1],sum(v*H[M^s] for s,v in enumerate(H))
def H(X,Y,m,n):
    h,a=X;k,b=Y
    return (a*2**m-h)*(b*2**n-k)+2*h*k
def inverse(X,Y,N):
    h,a=X;k,b=Y
    M=N-2*h*k
    if M<=0:return
    for d in range(1,isqrt(M)+1):
        if M%d:continue
        ds=(d,M//d)
        for x,y in dict.fromkeys((ds,ds[::-1])):
            u,r=divmod(x+h,a);v,s=divmod(y+k,b)
            if r or s or not u or not v or u&(u-1) or v&(v-1):continue
            yield u.bit_length()-1,v.bit_length()-1

if __name__=='__main__':
    assert [profile(X) for X in (A,B,C,D,M79,M85,M91)]==[(17,224),(3,96),(241,3760),(195,3720),(79,2176),(85,2176),(91,2176)]
    for m in range(5):
        for X,Y,ans in ((A,TT(2),672*2**m-17),(B,TT(3),672*2**m-15),(TT(4),D,56400*2**m-3135),(C,TT(4),56400*2**m-3133)):
            p,q=profile(X),profile(Y)
            assert counts(host(X,Y,m))[-1]==ans==H(p,q,m,0)
            assert (m,0) in tuple(inverse(p,q,ans))
    assert counts(host(M85,M85))[-1]-counts(host(M79,M91))[-1]==108
    for m in range(10):
        assert H((85,2176),(85,2176),m,m)-H((79,2176),(91,2176),m,m)==108
    for X,Y in (((1,2),(3,12)),((17,224),(195,3720))):
        for m in range(5):
            for n in range(5):
                N=H(X,Y,m,n)
                assert (m,n) in tuple(inverse(X,Y,N))
                d=min(m,n)
                assert (N-3*X[0]*Y[0])%(1<<d)==0
    print('PASS seven literal profiles, 20 literal dyadic unit endpoints, literal fixed 108 gap, and exact height inverses')
