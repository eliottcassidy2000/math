"""Integer proof of matching-width and color-branch bounds."""
from fractions import Fraction as F
from math import factorial


def floorroot(x,k):
    hi=1<<((x.bit_length()+k-1)//k)
    lo=0
    while hi-lo>1:
        b=(lo+hi)//2
        if b**k<=x: lo=b
        else: hi=b
    return lo


# upper bound a/x**(119/125), without floating-point
def upper(a,x,Q=10**8):
    rhs=(a*Q)**125
    lhs=x**119
    if rhs<=lhs: return F(1,Q)
    assert rhs<=lhs*Q**125
    i,j=1,Q
    while j-i>1:
        t=(i+j)//2
        if lhs*t**125>=rhs: j=t
        else: i=t
    return F(j,Q)


def compute():
    N=256
    facts=[factorial(n) for n in range(N+1)]
    omega=[0]+[floorroot(facts[k+1]**2//(2*(k+1)),k+1) for k in range(1,N-1)]
    g=[1]+[max([facts[t]]+[facts[t]//facts[t-k]*omega[k]**(t-k)
                           for k in range(1,t+1)]) for t in range(1,N-2)]
    z=F(7,8)
    assert 16**12500*15**4760*5**11781 < 35**12500*256**4041
    tail=3*z**256*(256/(1-z)+z/(1-z)**2)
    assert tail<F(1,10**9)
    upperK={}
    for m,goal in ((3,255),(4,83),(5,49),(6,34)):
        K=tail
        for h in range(2,N):
            t=h-2
            G=max([m**(2*t),(m**4+4*m*m)**(t//2)*m**(2*(t%2))]+
                  [m**(2*(t-a)+1)*g[a] for a in range(1,t+1)])
            future=max(factorial(m+h-1)//factorial(m-1),m*(m+1)*G)
            K+=upper(facts[h]-facts[h-1],future)
        assert K < F(goal,1000),(m,K,float(K))
        upperK[m]=F(goal,1000)
        print('K_%d < %s (compiled bound %.9f)'%(m,float(upperK[m]),float(K)))
    zm={m:upper(1,m) for m in range(3,7)}
    U={3:F(99,100),4:F(49,100),5:F(35,100),6:F(3,10),
       7:F(3,10),8:F(3,10)}
    for m in (6,5,4,3):
        x=zm[m]
        A=x+x**3+x**4/(1-2*x)
        B=A+x*U[m+1]/(1-U[m+1])+upperK[m]*U[m+2]/(1-U[m+2])
        assert B<U[m]
        print('S_%d inductive upper bound %.9f < %.2f'%(m,float(B),float(U[m])))
    return upperK


if __name__=='__main__':
    compute()
