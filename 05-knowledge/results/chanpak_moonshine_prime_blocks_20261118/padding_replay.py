"""Replay the switch identity and the equal/unequal padding certificates."""
from replay import dp, transpose, pad, profile, prime, R, wrap
from frame_replay import join, add, cut

cases=(
 ((62,40,26,48,34,4),4,1,1,(57,75,3696)),
 ((110,88,114,116,65,82,0),48,2,2,(709,529,107584)),
 ((128,241,211,119,161,197,17,72),55,3,3,(15029,15461,6443776)),
 ((42,48,27,18,1,28),26,2,3,(905,725,138368)),
 ((22,244,40,99,204,81,5,109),10,3,2,(7013,8813,2032160)),
 ((138,120,211,52,161,197,25,74),33,-1,-1,(69975,69579,11628288)),
)
triangle=(2,4,1)

def slots(A,k):
    h=[1]+[sum(w) for w in dp(A)[1:]];f=len(h)-1
    bits=[1]+[0]*f
    result=[bits]
    for j in range(1,k+1):
        nxt=[]
        for s in range(f+1):
            u=s;t=0
            while True:
                t+=h[u]*bits[s^u]
                if not u:break
                u=(u-1)&s
            nxt.append(t)
        result.append(nxt);bits=nxt
    return result

def rooted(U,I):
    E=dp(U);S=dp(transpose(U));f=len(E)-1
    h=slots(U,3)
    signed=[]
    for q in range(4):
        a=[[0,0],[0,0]]
        for x in range(1,f+1):
            v=h[q][f^x]
            for iy,Y in enumerate((E,S)):
                t=v*sum(z for i,z in enumerate(Y[x]) if I>>i&1)
                a[iy][0]+=t;a[iy][1]+=t*(-1)**x.bit_count()
        signed.append(a)
    for q in range(1,4):
        assert signed[q][0][1]==-signed[q-1][1][0]
        assert signed[q][1][1]==-signed[q-1][0][0]
    return [tuple(x[0] for x in a) for a in signed]

def materialize(U,I,a,b):
    L=triangle if a<0 else pad((),a)
    R=triangle if b<0 else pad((),b)
    X=join(L,join(U,R));T=add(X,I<<len(L));Tp=cut(T,1<<len(X))
    return (L,R,T,Tp)

def check():
    for U,I,a,b,answer in cases:
        L,P,T,Tp=materialize(U,I,a,b)
        v=rooted(U,I);e,s=v[1];ee,ss=v[2]
        X=slots(L,3);Y=slots(P,3)
        l=[z[-1] for z in X];r=[z[-1] for z in Y]
        h,B=profile(T);hh,BB=profile(Tp)
        assert (h,hh,B)==answer and B==BB
        assert hh-h==r[1]*l[2]*e-l[1]*r[2]*s
        assert BB-B==2*(r[2]*(l[3]*ee-l[1]*s)-l[2]*(r[3]*ss-r[1]*e))
        print(a,b, (e,s,ee,ss),(h,hh,B))
    U,I,a,b,(h,hh,B)=cases[4]
    _,_,T,Tp=materialize(U,I,a,b)
    for X,y in [(T,h),(Tp,hh)]:
        p=6*B*2**16+y
        assert prime(p)
        for scale in (1,7):
            D,c,_,_=R[(1,scale)]
            assert len(wrap(D,c,pad(X,16)))==len(D)+len(X)+17
        print('same-lattice skew-pad prime and partner',p,7*p)

if __name__=='__main__':check()
