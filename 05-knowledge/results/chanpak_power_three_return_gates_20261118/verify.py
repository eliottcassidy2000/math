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

rails = [
 (81,3,1,(2,4,9,3,15,31),4,(5,108,11,2)),
 (81,3,-1,(2,24,3,21,37,15,63),47,(25,108,49,2)),
 (243,3,-1,(50,32,3,7,14,28),47,(37,81,73,0)),
 (243,3,1,(30,100,96,102,110,1,33),100,(47,324,95,2)),
 (729,3,1,(60,113,74,34,12,20,57),63,(109,243,219,0)),
 (729,3,-1,(18,44,89,65,106,77,3),95,(109,243,217,0)),
 (2187,9,-1,(38,4,104,3,79,26,43),23,(81,486,323,1)),
 (2187,9,1,(20,109,120,1,42,9,57),11,(59,486,237,1)),
 (81,9,-1,(2,4,17,7,43,15),18,(15,144,59,4)),
]
boundary = {
 161:(2,160,3,71,15,29,55,125),
 163:(2,96,3,135,15,29,61,119),
 487:(122,104,147,100,170,132,52,75),
 973:(284,437,272,118,288,261,183,61,200),
 4373:(978,376,403,373,256,853,788,122,640,158),
}

def ck(ok,*args):
    if not ok:
        raise RuntimeError(args)
def F(a):
    if not a:
        return 1
    b=a//2
    return F(b)+(0 if a&1 else b+(b&1))
def shape(D,q,a,t):
    w=bin(a)[2:]+'0'*t
    return apex(ord(TT(len(w)),D),sum(int(x)<<i for i,x in enumerate(w))
                | (q<<len(w)))
def signed_R(D,q):
    H,Z=cols(D)
    return sum((-1)**i.bit_count()*H[(len(H)-1)^i]*
               sum(n for j,n in enumerate(r) if not q>>j&1)
               for i,r in enumerate(Z))

def check():
    for N,a,s,D,q,(h,W,K,e) in rails:
        hh,R,KK=sig(D,q)
        ck((hh,hh+R,KK)==(h,W,K) and a*W==N*2**e
           and K-h*(a-F(a))==s,N,a,s)
        for t in range(4):
            T=shape(D,q,a,t)
            if t:
                expected=N*2**(e+t)+s
            else:
                expected=N*2**e+s-h-signed_R(D,q)
            ck(strong(T) and cols(T)[0][-1]==expected,N,a,s,t)
    for N,T in boundary.items():
        ck(strong(T) and cols(T)[0][-1]==N,N)
    m35=(8,37,1,6,15,29)
    m35=tuple((a<<1) if i==0 else ((a&(1<<0))*3)|((a&~1)<<1)
              for i,a in enumerate(m35))
    m35=(m35[0]|2,m35[0])+m35[1:]
    ck(cols(m35)[0][-1]==35)
    ck(cols((2,4,9,3))[0][-1]==5 and 35*5**3==4375)
    for D,q,sgn in [(*plus,1),(*minus,-1)]:
        for k in range(8):
            h,R,K=sig(D,q)
            W=36*3**k if sgn==1 else 24*3**k+16*2**k-4
            ck(h+R==W and K-2*h==sgn,k,sgn)
            if sgn==1:
                ck(h==12*2**k-1 and strong(apex(D,q)))
                ck(cols(shape(D,q,3,1))[0][-1]==3**(k+3)*8+1,k)
            D=apex(D,q)
    print("PASS: marked signatures, both signs through 2187, quotient-product boundary, and clone reads")
if __name__=='__main__':
    check()
