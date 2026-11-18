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

from math import prod
def ck(ok,*context):
    if not ok:raise RuntimeError(context)
def trial():
    seed=(30,156,8,128,12,671,959,277,575,159)
    q=31
    D=seed
    ck([(i,j) for i in range(10) for j in range(10)
        if q>>i&1 and not q>>j&1 and D[i]>>j&1]==[(1,7),(3,7)])
    for k in range(6):
        h,R,K=sig(D,q)
        ck(h==90*2**k+1 and h+R==1224*3**k and K-2*h==-1,k)
        E=apex(D,q)
        reach=1<<(len(E)-1)
        for _ in E:
            for i,r in enumerate(E):
                if r&reach:reach|=1<<i
        ck(reach==(1<<len(E))-1)
        if k<2:
            for m in (3,4):
                L=TT(m)
                T=apex(ord(L,D),3|(q<<m))
                C=17*3**(k+3)*2**(m+1)
                ck(strong(T) and cols(T)[0][-1]==C-1)
                ck((3*(C-1)+1)//2==3*C//2-1 and (3*(C-1)+1)%4==2)
        D=E
    ck(7343==7*1049 and all(1049%d for d in range(2,33)))
F = {
 3:(2,4,1), 5:(4,9,10,1), 11:(2,12,25,17,3),
 113:(96,45,73,17,103,12,42), 73:(94,92,112,4,40,11,56),
 313:(72,205,241,196,75,27,32,113),337:(156,237,72,208,102,141,33,84),
 239:(54,228,32,23,166,136,189,13),577:(108,197,184,178,131,18,60,97),
}

def run():
    trial()
    for v,a in F.items():
        ck(all(not row>>i&1 and row < 1<<len(a) for i,row in enumerate(a)))
        ck(all((a[i]>>j&1)+(a[j]>>i&1)==1 for i in range(len(a)) for j in range(i)))
        ck(cols(a)[0][-1]==v,v)
    for v,seq in {57121:(239,239),57123:(3,3,11,577),
                  114243:(3,113,337),114245:(5,73,313)}.items():
        ck(all(i in F for i in seq) and prod(seq)==v,v)
    ck(239**2+1==57122==2*13**4 and 2*57122==338**2)
    ck((239-1)//2==7*17 and 119**2+120**2==169**2)
    x=y=1
    for a in (1,7,41,239):
        ck(x==a and x*x-2*y*y==-1)
        x,y=3*x+4*y,2*x+3*y
    ck(complex(5,1)**4==complex(476,480)==2*complex(1,1)*complex(239,1))
    M=[[1-2*(i==j) for j in range(4)] for i in range(4)]
    J=[0,2,3,1]
    M=[[M[i][j]*(1 if j else -1) for j in J] for i in range(4)]
    ck(all(M[i][i]==1 and all(M[i][j]==-M[j][i] for j in range(i))
           for i in range(4)))
    ck(all(sum(x*y for x,y in zip(a,b))==4*(i==j)
           for i,a in enumerate(M) for j,b in enumerate(M)))
    print("PASS: exact minus growth, literal factor witnesses, Pell and circulant order-4 controls")
if __name__=='__main__':
    run()
