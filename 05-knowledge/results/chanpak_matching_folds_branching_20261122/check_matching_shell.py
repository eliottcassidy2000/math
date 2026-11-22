from itertools import permutations
from math import comb, factorial

def count(adj, cycle=False):
    n=len(adj)
    if not n: return 1
    ds=[{} for _ in range(1<<n)]
    for i in (range(1) if cycle else range(n)): ds[1<<i][i]=1
    for s, d in enumerate(ds):
        for i, val in d.items():
            for j in range(n):
                if not (s>>j & 1) and adj[i][j]:
                    de=ds[s | 1<<j]
                    de[j]=de.get(j,0)+val
    return sum(v for k,v in ds[-1].items() if not cycle or adj[k][0])

def orient(n, f):
    return [[False if i==j else f(i,j) if i<j else not f(j,i) for j in range(n)] for i in range(n)]

def matrix_gadget(m, k, mask):
    def f(i,j):
        if i//m==0:
            return bool(mask >> (i*m+j-m) & 1) if j//m==1 else j<m
        return True
    return orient(2*m+k,f)

def surj(k,r):
    return sum((-1)**(r-j)*comb(r,j)*j**k for j in range(r+1))

def recover(hs):
    a=[]
    for k,h in enumerate(hs,1):
        y=h-sum(z*surj(k,r) for r,z in enumerate(a,1))
        if y%factorial(k): raise ValueError('nonintegral')
        a.append(y//factorial(k))
    return a

def indep_perm(m, mask):
    return sum(all(mask>>(i*m+j)&1 for i,j in enumerate(p)) for p in permutations(range(m)))

def k2(A):
    n=len(A)
    def induced(vertices): return [[A[i][j] for j in vertices] for i in vertices]
    return sum(count(induced([i for i in range(n) if s>>i&1]))
               *count(induced([i for i in range(n) if not(s>>i&1)])) for s in range(1<<n))

def shell(A,k):
    n=len(A)
    def f(i,j):
        if (i,j)==(0,2*n+k+1): return False
        if 1<=i<j<=n: return A[i-1][j-1]
        if n+k+1<=i<j<=2*n+k: return A[j-(n+k+1)][i-(n+k+1)]
        return True
    return orient(2*n+k+2,f)

def main():
    npairs=0
    for m in range(1,4):
        for mask in range(1<<(m*m)):
            p=indep_perm(m,mask)
            for cyclic in (False,True):
                r=m if cyclic else m+1
                a=recover([count(matrix_gadget(m,k,mask),cyclic) for k in range(1,r+1)])
                if a[-1] != factorial(r-1)*p: raise ValueError((m,mask,a,p))
            npairs+=1
    nshells=0
    for n in range(1,5):
        for mask in range(1<<comb(n,2)):
            ix={(i,j):q for q,(i,j) in enumerate((i,j) for i in range(n) for j in range(i+1,n))}
            A=orient(n,lambda i,j:bool(mask>>ix[i,j]&1))
            h=count(A)
            K=k2(A)
            for k in (1,2,3):
                if count(shell(A,k)) != h*h+2**k*K*K: raise ValueError((n,mask,k))
                nshells+=1
    print(f'matrix instances {npairs}, each tested for path and cycle; symmetry shells {nshells}')
if __name__=='__main__': main()
