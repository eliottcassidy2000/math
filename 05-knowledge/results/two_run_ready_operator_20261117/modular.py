from replay import hp,move,rooted

def append(A,r):
    n=len(A)
    return tuple(x|((i!=r)<<n) for i,x in enumerate(A))+(1<<r,)

def initial():
    A=move((6,4,0));r=6
    return append(append(A,r),r)

def evaluate(A,k,modulus):
    h=hp(A);N=len(h);P=[1]+[0]*(N-1)
    a=[P[-1]]
    for _ in range(k):
        Q=[]
        for s in range(N):
            q=0;i=s
            while True:
                q+=h[i]*P[s^i]
                if not i:break
                i=(i-1)&s
            Q.append(q%modulus if modulus else q)
        P=Q;a.append(P[-1])
    return tuple(a)

def c(q):
    return q**5+2*q**3,2*q**4+3*q**3+6*q*q+3*q,-2*q**4+3*q**3-6*q*q+3*q

def transform(P,m):
    return tuple(sum(a*P[b%m] for a,b in zip(c(q),(q,q+1,q-1)))%m for q in range(m))

if __name__=='__main__':
    B=initial()
    assert hp(B)[-1]==511
    assert rooted(B,6)==(511,17,264,281)
    s0=evaluate(B,8,9)
    print("seed",B,"mod 9 state",s0)
    S={}
    cur=s0;ph=[]
    while cur not in S:
        S[cur]=len(ph)
        qmod9=(4*cur[1])%9
        ph.append(qmod9 if qmod9%2 else qmod9+9)
        cur=transform(cur,9)
    print('p=H/7 mod 18; preperiod',ph[:S[cur]],'cycle',ph[S[cur]:],flush=True)
    assert S[cur]==0 and ph==[1,17,17,7,1,5,17,1,1,11,17,13]
    Q=evaluate(B,4,None);print("exact seed path-color values",Q)
    for n in range(3):
        print("exact count",n,Q[1])
        Q=tuple(sum(a*Q[j] for a,j in zip(c(q),(q,q+1,q-1))) if q else 0 for q in range(len(Q)-1))
    # Direct evaluate a small predecessor and compare the full shifted operator.
    for t in ((0,), (2,0),(2,4,1),(2,4,9,3)):
        for m in (5,7,9):
            assert evaluate(move(t),m-1,m)==transform(evaluate(t,m-1,m),m)
