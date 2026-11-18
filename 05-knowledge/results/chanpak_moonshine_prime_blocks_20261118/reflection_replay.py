from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent
sys.path[:0]=[str(HERE.parent/'ans_inverse'),str(HERE.parent/'amm12592_selector')]
from fixtures import dp,rev,state,H,forward
from benchmark import BASE,FROZEN,union,primes_to

def rooted_last(A,r):
    J=tuple(i for i in range(len(A)) if i!=r)+(r,)
    return tuple(sum((A[i]>>j&1)<<k for k,j in enumerate(J)) for i in J)

def arm(A,r):
    T=rooted_last(A,r)
    f=(1<<(len(T)-1))-1
    Q=tuple(x&f for x in T[:-1])
    E,F=dp(Q),dp(rev(Q))
    h=[1]+[sum(row) for row in E[1:]]
    b=sum(h[X]*h[f^X] for X in range(f+1))
    U=V=Ub=Vb=0
    for S in range(1,f+1):
        sign=(-1)**S.bit_count()
        a=sum(z for j,z in enumerate(F[S]) if T[-1]>>j&1)*h[f^S]
        c=sum(z for j,z in enumerate(E[S]) if not (T[-1]>>j&1))*h[f^S]
        U+=a;V+=sign*a;Ub+=c;Vb+=sign*c
    return H(T),h[f],b,U,V,Ub,Vb

def glue(A,B):
    n,m=len(A)-1,len(B)-1;f=(1<<n)-1;g=(1<<m)-1
    return (tuple((r&f)|(g<<n)|((r>>n)<<(n+m)) for r in A[:-1])
            +tuple(((r&g)<<n)|((r>>m)<<(n+m)) for r in B[:-1])
            +(A[-1]|(B[-1]<<n),))

def transitive_arm(I):
    m=I.bit_length();f=(1<<m)-1
    return tuple((f^((1<<(i+1))-1))|(((I>>i)&1)<<m) for i in range(m))+(f^I,)

def apex_increment(I):
    return sum(1<<max(1,j-i-1) for j in range(I.bit_length()) if I>>j&1
               for i in range(j) if not (I>>i&1))

def count(A,I):
    h,g,b,U,V,ub,vb=A
    return h+g*apex_increment(I)+U*I-V*(I&1)

A0=((126,28,56,112,96,66,6),3)
A1=((126,28,40,112,100,66,6),4)
S=(41,15,240,63,-9,99,-15)
def two_digits():
    C0=((62,124,24,112,96,68,5),2)
    C1=((16,65,3,71,14,31,53),1)
    C=(65,17,224,72,-6,135,-11)
    assert arm(*C0)==C and arm(*C1)==(69,*C[1:])
    for i,A in enumerate((A0,A1)):
        for j,B in enumerate((C0,C1)):
            T=glue(rooted_last(*A),rooted_last(*B))
            v=len(T)-1;h=10021+34*i+60*j
            assert arm(T,v)==(h,255,53760,15192,-90,34083,-255)
            assert state(T,v)==(h,165,15447,255)
    print('independent signed digits',10021,10055,10081,10115)

def main():
    two_digits()
    assert arm(*A0)==S and arm(*A1)==(43,*S[1:])
    A,B=(rooted_last(*w) for w in (A0,A1))
    for I in range(1,33):
        R=transitive_arm(I)
        X,Y=glue(A,R),glue(B,R);r=len(X)-1
        u,v=arm(X,r),arm(Y,r)
        assert u[0]==count(S,I) and v==(u[0]+2,*u[1:])
        assert state(forward(Y,r,(0,1,0,1)),r)==tuple(
            n+(2 if j==0 else 0) for j,n in enumerate(state(forward(X,r,(0,1,0,1)),r)))
    p=primes_to(262144)
    for lo,hi in ((4096,16384),(16384,65536),(65536,262144)):
        X=union(BASE+FROZEN,lo,hi)
        rows=[]
        for I in range(1,(hi-S[0])//S[3]+1):
            j=count(S,I)
            if lo<=j<hi-2 and (j in X)!=(j+2 in X):
                missing=j+2 if j in X else j
                rows.append((I,j,missing))
        odds={q[2] for q in rows};prime={z for z in odds if p[z]}
        print(lo,hi,'matched transfers',len(odds),'primes',len(prime),'first',rows[:4],
                'first prime',next((q for q in rows if q[2] in prime),None))
        assert (len(odds),len(prime))=={4096:(37,10),16384:(121,41),65536:(552,152)}[lo]
        if lo==65536:
            I,J,miss=next(q for q in rows if q[2] in prime)
            R=transitive_arm(I)
            for T,j in ((glue(A,R),J),(glue(B,R),J+2)):
                r=len(T)-1
                assert H(T)==j
                assert H(forward(T,r,(0,)))==j
                assert state(T,r)==(j,15,76107,15)
                print('literal heldout',I,j,state(T,r))
if __name__=='__main__': main()
