"""Cross-check the target examples, exact phase carry, and shared-height law."""
from math import isqrt
from replay_independent import H45, contexts, paths, value, g, solve
import replay_diagonal as diag

DEF=(45,432,129,258,15,30,137)
# An independently movable window outside every old pure complementary chart.
def old_window(N):
    for k in range(4,(N+1).bit_length()):
        M=2**k;u=isqrt(N)//M
        B=M*M*u*(u+1)+M*M//4-u+2*g(u)+(u&1)-3
        lo,hi=0,M//2-1-(u==0)
        while lo<=hi:
            q=(lo+hi)//2
            v=B+2*g(q)-q*q-2*q+(q&1)-2*(q>>(k-2))
            if v>N:lo=q+1
            elif v<N:hi=q-1
            else:return True
    return False

def prime(n):return n>=2 and all(n%d for d in range(2,isqrt(n)+1))
for i,N in enumerate([7605,7607,7609]):
    assert value(1,12,DEF)+2*(i-1)==N
    assert (1,12,i) in set(solve(N,DEF))
    assert sum(paths(contexts(H45,1,12,i))[-1])==N and not old_window(N)
assert prime(7607) and prime(1087) and 7609==7*1087
assert not any(old_window(N) for N in range(5953,5962,2))

# q carries at any dyadic tilt R; width bound makes the shifted tails legal.
for old,U in diag.tilt7:
    h,f,P,Q,S,T,K=diag.profile(old,U)
    for M in (4,8,16,32,64):
        for u in (1,2,7):
            X=M*u
            for R in (2,4,8):
                base=f*X*(R*X-1)+Q*X+P*(R*X-1)+h*(g(X)+g(R*X-1)-2)+K+S
                for q in range(-M//2,M//2):
                    t=q if q>=0 else -q-1
                    d=2*g(t)-t+(t&1)-2
                    x,y=X+q,R*(X-q)-1
                    exact=f*x*y+Q*x+P*y+h*(g(x)+g(y)-2)+K+T*(x&1)+S*(y&1)
                    assert exact==base-f*q*(R*q+1)+q*(Q-R*P)+h*d+T*(q&1)

# A single shared baseline across *different* source F2 and outer widths.
synchronised=[
 (4,8,(28,25,10,16,4,31),45),
 (4,8,(28,25,10,16,4,31),46),
 (2,8,(124,121,58,0,8,88,28),74),
 (4,8,(28,25,10,16,4,31),54),
 (4,8,(28,1,2,6,14,31),43),
]
def val(u,a,b,h,D,C,delta):
    return D*u*u+C*u+h*(g(u)+g(u-1)+(u&1)-2)+delta
for i,(a,b,rows,U) in enumerate(synchronised):
    h,f,P,Q,S,T,K=diag.profile(rows,U)
    delta=K-P+S
    D,C=a*b*f,b*P+a*(Q-f+h)-h
    assert (h,D,C,delta)==(9,5376,627,-59+2*i)
    assert diag.count(diag.physical(rows,U,a,b-1))==5953+2*i
    for u in (1,2,3,17,128):
        x,y=a*u,b*u-1
        assert val(u,a,b,h,D,C,delta)==value(x,y,(h,f,P,Q,S,T,K))
        assert val(u+1,a,b,h,D,C,delta)>val(u,a,b,h,D,C,delta)
print("PASS: distinct exact target windows, general phase carry, and varied-width common anchor")

def tilted_windows(N):
    specs=(
      (diag.banks[0][1],1,1), (diag.banks[1][1],1,1),
      (diag.tilt,1,2), (diag.tilt7,1,4),
    )
    for bank,L,R in specs:
        for rows,U in bank:
            h,f,P,Q,S,T,K=diag.profile(rows,U)
            M=2
            while f*L*R*M*M<=2*N:
                a,b=L*M,R*M
                D,C=a*b*f,b*P+a*(Q-f+h)-h
                lo,hi=1,isqrt(2*N//D)
                while lo<=hi:
                    u=(lo+hi)//2
                    w=val(u,a,b,h,D,C,K-P+S)
                    if w==N:
                        yield a,b,u,diag.physical(rows,U,a*u,b*u-1)
                        break
                    if w<N:lo=u+1
                    else:hi=u-1
                M*=2

assert any((a,b,u)==(2,8,1) for a,b,u,_ in tilted_windows(5953))
if __name__=='__main__':
    import sys
    if len(sys.argv)>1:
        N=int(sys.argv[1])
        for a,b,u,rows in tilted_windows(N):
            print("TILTED WITNESS",N,a,b,u,rows)
