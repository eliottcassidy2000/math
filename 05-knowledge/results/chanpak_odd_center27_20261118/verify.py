from itertools import permutations
minus=((2,4,1,7,15),4)
plus=((2,20,9,3,13),13)
small={(1,-1):(8,1,67,38,15,23,59),(1,1):(4,65,18,7,11,31,61),
       (2,-1):(4,17,162,7,13,27,63,123),(2,1):(66,48,3,7,13,29,62)}

def cols(a):
    n=len(a);f=[[0]*n for _ in range(1<<n)]
    for j in range(n):f[1<<j][j]=1
    for S in range(1,1<<n):
        for j,v in enumerate(f[S]):
            U=a[j]&~S
            while U and v:
                b=U&-U;U-=b;f[S|b][b.bit_length()-1]+=v
    return [1]+[sum(q) for q in f[1:]],f
def apex(a,q):
    n=len(a)
    return tuple(v|((1-(q>>i&1))<<n) for i,v in enumerate(a))+(q,)
def TT(n):return tuple((1<<n)-(1<<(i+1)) for i in range(n))
def ord(a,b):return tuple(v|(((1<<len(b))-1)<<len(a)) for v in a)+tuple(v<<len(a) for v in b)
def sig(a,q):
    h,f=cols(a);M=len(h)-1
    R=sum(h[M^S]*sum(v for j,v in enumerate(row) if not(q>>j&1)) for S,row in enumerate(f))
    return h[-1],R,cols(apex(a,q))[0][-1]
def build(m,sign):
    if m<1 or sign not in (-1,1):raise ValueError('m≥1 and sign ±1')
    if m<3:return small[m,sign]
    a,q=minus if sign<0 else plus
    return apex(ord(TT(m),a),3|(q<<m))
def strong(T):
    n=len(T)
    for v in range(n):
        R=1<<v
        for _ in T:
            for j in range(n):
                if R>>j&1:R|=T[j]
        if R!=(1<<n)-1:return False
    return True

from fractions import Fraction as F

def ck(b,*msg):
    if not b:
        raise RuntimeError(msg)

def steps(a):
    high=a
    n=0
    while a!=1:
        a=3*a+1 if a&1 else a//2
        n+=1
        high=max(high,a)
    return n,high

def digits(a):
    return [a%3,a//3%3,a//9]
def code(a):
    return sum((v%3)*3**i for i,v in enumerate(a))
def sub(a,b):
    return code([x-y for x,y in zip(digits(a),digits(b))])
def mul(a,b):
    r=[0]*5
    for i,x in enumerate(digits(a)):
        for j,y in enumerate(digits(b)):
            r[i+j]+=x*y
    for i in [4,3]:
        r[i-3]+=r[i]
        r[i-2]+=r[i]
    return code(r[:3])

def had():
    Q={mul(x,x) for x in range(1,27)}
    ck(len(Q)==13 and 0 not in Q)
    H=[[1]*28]+[[-1]+[1 if x==y or sub(y,x) in Q else -1
                       for y in range(27)] for x in range(27)]
    for i,a in enumerate(H):
        for j,b in enumerate(H):
            ck(sum(x*y for x,y in zip(a,b))==28*(i==j),i,j)
            if i!=j:
                ck(H[i][j]==-H[j][i],i,j)

def atan(a):
    s=sum(F((-1)**i,(2*i+1)*a**(2*i+1)) for i in range(12))
    t=s+F(1,25*a**25)
    return s,t

def main():
    ck(sig(*minus)==(3,33,5) and sig(*plus)==(11,25,23),'profile')
    for m in range(1,10):
        for sign in (-1,1):
            T=build(m,sign)
            ck(strong(T) and cols(T)[0][-1]==27*2**m+sign,m,sign)
    for a,b in [(69,13),(71,107),(73,55)]:
        z=3*a+1
        ck(z//(z&-z)==b,a)
    ck((3*73-1)//2==109 and steps(27)==(111,9232),'Collatz')
    ck(pow(2,18,27)==1 and all(pow(2,j,27)!=1 for j in range(1,18)))
    for e in range(1,19):
        c=[(3*x+1)*pow(2,-e,27)%27 for x in range(27)]
        ck(len(set(c))==9 and all(c.count(x)==3 for x in set(c)))
    for d in range(1,6):
        f=[j for j in range(7**(d-1)) if (27*pow(8,j,7**d)+1)%7**d==0]
        ck(len(f)==1,d)
        if d==2:ck(f==[4])
    had()
    a,b=atan(5);c,d=atan(239)
    lo,hi=16*a-4*d,16*b-4*c
    def ratio(m):
        C=27*2**m
        return F(355*(C-1),113*(C+1))
    ck(ratio(19)<lo<hi<ratio(20))
    print("PASS: 18 literal centered hosts; core, Collatz and 7-adic checks; 28 Hadamard rows; Machin threshold")

if __name__=='__main__':
    main()
