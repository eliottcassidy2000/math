from build import TT,D4,htable,ordinal,host,transpose,cyclic
from itertools import product
from math import comb

def profiles(A,k):
    h=htable(A); M=len(h)-1; prev=[1]+[0]*M; ans=[prev[M]]
    for j in range(k):
        z=[0]*(M+1)
        for s in range(M+1):
            t=s
            while True:
                z[s]+=h[t]*prev[s^t]
                if t==0:break
                t=(t-1)&s
        prev=z;ans.append(z[M])
    return ans

def Q(P,j):return sum((-1)**(j-i)*comb(j,i)*P[i] for i in range(j+1))
def expr(B):
    b=len(B);q=profiles(B,b);f=[]
    for j in range(3,3*b+3):
      for start in range(3):
        counts=[0]*3
        for i in range(j):counts[(i+start)%3]+=1
        a,c,z=counts
        if not z or z>b:continue
        f.append((a,c,Q(q,z)))
    return f

def evaluate(j,m,f,mod):
    P=[0]+[pow(PD[i],j,mod)*pow(i,m+1,mod)%mod for i in range(1,7)]
    return sum(Q(P,a)*Q(P,c)*z for a,c,z in f)%mod
PD=profiles(D4,6)
Z3=(2,4,1)
# Keyed by H mod 8. Slots are indexed by H mod 7, from 0 through 6.
rows={
 1:(tuple(TT(3)),[(0,5),(1,3),(1,1),(1,5),(0,3),(0,2),(2,0)]),
 3:((14,20,24,18,1),[(0,2),(0,1),(0,4),(1,3),(0,0),(4,0),(0,3)]),
 5:(tuple(TT(2)),[(0,3),(2,2),(1,3),(0,1),(0,5),(0,0),(1,1)]),
 7:(Z3,[(1,3),(0,0),(0,2),(1,4),(0,5),(0,1),(0,4)])}

def decode(k,target):
    if k<1 or not target%2:raise ValueError('k>=1 and odd target')
    B, pairs=rows[target%8]; f=expr(B); j,m=pairs[target%7]
    step=6; mod=7
    for e in range(1,k):
        m,=[m+step*d for d in range(7) if evaluate(j,m+step*d,f,7*mod)==target%(7*mod)]
        step*=7;mod*=7
    return B,j,m

if __name__=='__main__':
    for r,(Z,pairs) in rows.items():
        if sum(profiles(Z,2)[1:])%8!=r:raise ArithmeticError('8 phase')
        zz=len(Z)
        if any((Z[i]>>j&1)!=(Z[zz-1-j]>>(zz-1-i)&1) for i in range(zz) for j in range(zz)):raise ArithmeticError('center reflection')
        for s,(j,m) in enumerate(pairs):
            h=evaluate(j,m,expr(Z),49);hs=evaluate(j,m+6,expr(Z),49)
            if h%7!=s or (hs-h)%49==0:raise ArithmeticError('bad seed')
        for j,m in product(range(2),repeat=2):
            X=ordinal(TT(1),TT(m))
            if j:X=ordinal(D4,X)
            A=cyclic(X,transpose(X),Z)
            if htable(A)[-1]%1009!=evaluate(j,m,expr(Z),1009):raise ArithmeticError('unfolded counts disagree')
    tested=0
    for k in range(1,5):
      for target in range(1,8*7**min(k,2),2):
        target+=8*(target*233)*7**max(k-2,0)
        B,j,m=decode(k,target)
        if evaluate(j,m,expr(B),7**k)!=target%7**k:raise ArithmeticError('lift')
        if sum(profiles(B,2)[1:])%8!=target%8:raise ArithmeticError('lift changed 8 phase')
        tested+=1
    print(tested,'CRT lifts, 28 validated seed profiles, and 16 independent full-DP checks')
