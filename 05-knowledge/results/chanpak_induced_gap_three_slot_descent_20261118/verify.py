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


from random import Random
from math import isqrt, factorial

def ck(ok,*tags):
    if not ok:raise RuntimeError(tags)
def prime(n):
    return n>=2 and all(n%d for d in range(2,isqrt(n)+1))
def slot_profile(B,J=3):
    H=cols(B)[0];a=[0]*len(H)
    a[0]=1
    vec=[]
    for _ in range(J):
        z=[0]*len(H)
        for S in range(len(H)):
            t=S
            while True:
                z[S]+=a[S^t]*H[t]
                if not t:break
                t=(t-1)&S
        vec.append(z[-1]);a=z
    return vec

minus17=((30,156,8,128,12,671,959,277,575,159),31)

def host(B,m,kind):
    D,q=[plus,minus,minus17][kind]
    a=len(D)
    # Uniform old clone-core orientations, then a last clone as the mark.
    E=tuple(v if q>>i&1 else v|(((1<<len(B))-1)<<a)
            for i,v in enumerate(D))
    E+=tuple((row<<a)|q for row in B)
    return apex(ord(TT(m),E),3|(q<<m))

def read(B,m,kind):
    h,p,q=slot_profile(B)
    if kind==0:return h+27*2**m*q
    if kind==1:return 2**m*(18*q+12*p-3*h)-h
    return 17*27*2**(m+1)*q-h

def sample(n,r):
    d=[0]*n
    for i in range(n):
        for j in range(i):
            x,y=(i,j) if r.randrange(2) else (j,i)
            d[x]|=1<<y
    return tuple(d)
def main():
    rng=Random(26)
    found={'p':((2,4,9,3),4,50549),
           '7p':((28,5,32,86,38,75,23),3,1905869)}
    for label,(B,m,N) in found.items():
        T=host(B,m,0)
        ck(cols(T)[0][-1]==N and strong(T),label)
        ck(prime(N) if label=='p' else (N%7==0 and prime(N//7)))
        print(label,B,m,slot_profile(B),N)
    for n in range(1,6):
        for _ in range(4):
            B=sample(n,rng)
            for m in (3,4):
                for kind in range(3 if n<=3 and m==3 else 2):
                    T=host(B,m,kind)
                    ck(strong(T) and cols(T)[0][-1]==read(B,m,kind),B,m,kind)
    P=[((8,9,3,100,15,23,55),(144,5,161,7,78,27,47,122)),
       ((32,97,19,7,11,28,61),(112,97,67,7,14,28,56,127))]
    v=[[x*y for x,y in zip(slot_profile(a,4),slot_profile(b,4))] for a,b in P]
    ck(v==[[17325,5105184,317017071,8535609600],
           [17325,5105184,317017071,8540125440]],v)
    # Every arbitrary induced inclusion of the labeled 5-vertex tournaments.
    n=5
    edge=[(i,j) for i in range(n) for j in range(i)]
    pairs=set()
    for i in range(1<<len(edge)):
        d=[0]*n
        for j,(a,b) in enumerate(edge):
            if i>>j&1:d[a]|=1<<b
            else:d[b]|=1<<a
        H=cols(d)[0];full=H[-1]
        for y in H:
            delta=full-y
            if delta==2:pairs.add((y,full))
            if delta>0:
                g=delta//2
                ck(delta%2==0 and full<=g*factorial(2*g+2))
    ck(pairs=={(1,3),(3,5)},pairs)
    print('PASS: three module reads, prime-axis witnesses, induced gap on all labeled order-five inclusions.')
if __name__=='__main__':main()
