from itertools import combinations, permutations
from math import comb, factorial

def tournament(n,m):
    A=[0]*n
    for b,(i,j) in enumerate(combinations(range(n),2)):
        x,y=(i,j) if m>>b&1 else (j,i)
        A[x]|=1<<y
    return A

def steps(n,m):
    A=tournament(n,m)
    d=[[[0]*n for _ in range(n)] for _ in range(1<<n)]
    for i in range(n):d[1<<i][i][0]=1
    for S in range(1,1<<n):
        for i,poly in enumerate(d[S]):
            if any(poly):
                for j in range(n):
                    if not S>>j&1:
                        b=int(not(A[i]>>j&1))
                        q=d[S|1<<j][j]
                        for k,v in enumerate(poly[:n-b]):q[k+b]+=v
    W=[sum(c[k] for c in d[-1]) for k in range(n)]
    assert W==W[::-1] and sum(W)==factorial(n)
    return W

def choose(n,r):
    return comb(n,r) if 0<=r<=n else 0

def shells(n,m):
    W=steps(n,m);D=comb(n,2)
    return [sum(w*choose(D+1-n,r-k) for k,w in enumerate(W)) for r in range(D+1)]

def kval(D,r,k):
    return sum((-1)**b*choose(k,b)*choose(D-k,r-b) for b in range(r+1))

def verify_field(n,m,p):
    D=comb(n,2);S=shells(n,m)
    q=[k for k in range(0,n,2) if n-k<p]
    B=[[kval(D,r,k)%p for k in q]+[S[r]%p] for r in range(len(q))]
    for j in range(len(q)):
        k=next(k for k in range(j,len(q)) if B[k][j])
        B[j],B[k]=B[k],B[j]
        g=pow(B[j][j],-1,p)
        B[j]=[x*g%p for x in B[j]]
        for i in range(len(q)):
            if i!=j:
                g=B[i][j]
                B[i]=[(x-g*y)%p for x,y in zip(B[i],B[j])]
    for r in range(D+1):
        assert S[r]%p==sum(B[i][-1]*kval(D,r,k) for i,k in enumerate(q))%p

def digits(D,p):
    R=[];r=1
    while r<=D:R.append(r);r*=p
    signs={tuple(kval(D,r,k)%p for r in R) for k in range(D+1)}
    assert len(signs)==D+1

def main():
    for n,m,h in [(4,2,3),(8,153661494,315)]:
        D=comb(n,2);W=steps(n,m)
        assert W==[h*comb(n-1,r) for r in range(n)]
        # Independent brute-force vertex-order count for these exact fixtures.
        A=tournament(n,m);W1=[0]*n
        for seq in permutations(range(n)):
            W1[sum(not(A[i]>>j&1) for i,j in zip(seq,seq[1:]))]+=1
        assert W1==W
        assert shells(n,m)==[h*comb(D,r) for r in range(D+1)]
        print('perfect ball',n,m,h)
    a,b=[shells(9,m) for m in (31070013745,45517909057)]
    assert a[:4]==[1215,44094,779358,8933310] and b[:4]==[1589,56458,975994,10935850]
    assert [x%11 for x in a[:3]]==[x%11 for x in b[:3]] and a[3]%11!=b[3]%11
    S=shells(7,(1<<21)-1)
    assert S[:4]==[1,135,3096,33336]
    assert [S[j]%3 for j in range(4)]==[1,0,0,0] and comb(21,3)%3==1
    for p in (3,5,7,11,13,19):
        digits(36,p)
        for n,m in [(7,153662),(8,153661494),(9,31070013745)]:
            verify_field(n,m,p)
    mm=lambda fs:10^sum(1<<e for e in fs)
    assert [steps(5,mm(fs))[0] for fs in ([],[7],[6],[7,6],[5],[7,5])]==[9,11,13,15,13,11]
    assert [steps(7,158750^x)[0] for x in (0,1<<16,1<<14,(1<<16)^(1<<14))]==[131,133,135,137]
    print('modular shell forecasts and radius-two collision verified')

if __name__=='__main__':main()
