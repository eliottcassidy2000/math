"""Independent Hamilton subset DP and literal leaf operations."""
from random import Random

def leaf(rows,q):
    n=len(rows)
    return tuple(v|(0 if j==q else 1<<n) for j,v in enumerate(rows))+(1<<q,)

def profiles(rows):
    n=len(rows); f=(1<<n)-1
    a=[[0]*n for s in range(f+1)];h=[1]+[0]*f
    for s in range(1,f+1):
        for i in range(n):
            if not s>>i&1:continue
            t=s^(1<<i)
            a[s][i]=sum(a[t][j] for j in range(n) if t&rows[i]&(1<<j)) if t else 1
        h[s]=sum(a[s])
    out=[]
    for i in range(n):
        S=d=0
        for s in range(1,f+1):
            v=a[s][i]*h[f^s]; S+=v;d+=v if s.bit_count()%2 else -v
        K=h[f^(1<<i)]
        out.append((S,d,K,S-d-h[-1]))
    return h[-1],out

def delta_fan(a,b,L,K,d,eta):
    return L*2**a*(2*3**b-1-2**(b+1))+(K-eta)*2**b+(K+eta)*(b==0)+(d-K)*(a==0)+eta

def main():
    G=Random(19721)
    seen=0
    for n in range(1,7):
        for _ in range(25):
            r=[0]*n
            for j in range(n):
                for i in range(j):
                    if G.randrange(2):r[i]|=1<<j
                    else:r[j]|=1<<i
            q=G.randrange(n)
            H,A=profiles(r)
            S,d,K,delta=A[q]; L=S-K;eta=delta+2*(d-K)
            for k in range(1,6):
                H0,A0=profiles(r)
                r=leaf(r,q); h,B=profiles(r)
                if h!=2**k*L-eta or B[q][1:]!=(K,K,eta) or B[-1][3]!=2*A0[q][1]-H0:
                    raise AssertionError((r,q,h,B[q]))
                for i in range(k):
                    want=delta_fan(i,k-i-1,L,K,d,eta)
                    if B[n+i][3]!=want:raise AssertionError((r,q,k,i,B[n+i],want))
                    seen+=1
    print('Endpoint DP checks for exact individual fan-vertex defects:',seen)
    R=(256,509,505,497,481,449,129,257,64)
    h,A=profiles(R)
    X=profiles(leaf(R,0));Y=profiles(leaf(leaf(R,0),0))
    if (h,A[0],X[0],X[1][0][3],Y[0])!=(5,(6,2,3,-1),9,-3,15):raise AssertionError(X)
    print('Unstable original minus-one port: counts 5,9,15; defect changes -1 to -3')
    R=(6,12,8,1)
    h,A=profiles(R)
    if A[0][1:]!=(1,1,1):raise AssertionError((h,A))
    for k in range(6):
        h,A=profiles(R)
        if h!=6*2**k-1 or A[0][3]!=1:raise AssertionError(A)
        R=leaf(R,0)
    print('Stable plus-one fan: first six counts 5,11,23,47,95,191')
if __name__=='__main__':main()
