"""E[N_j] for a uniformly random word of shape (jm, jx): N_j = number of r in [0,m) with all j windows fair.
Exact: E[N_j] = m * C(m,x)^j / C(jm, jx). Asymptotics: sqrt(j) m^((3-j)/2) (2 pi rho(1-rho))^(-(j-1)/2).
Monte Carlo for P(N_j >= 1) at large m (j = 3, 4) and the exact expectation check against brute force (small)."""
import random, sys
from math import comb, sqrt, pi
from itertools import combinations
def N_j(w,j):
    K=len(w); m=K//j; x=sum(w)//j; ww=w+w; cnt=0
    # prefix sums
    F=[0]
    for b in ww: F.append(F[-1]+b)
    for r in range(m):
        if all(F[r+(i+1)*m]-F[r+i*m]==x for i in range(j)): cnt+=1
    return cnt
# exact expectation check on all words (not necklaces) for small shapes
for (K,X,j) in [(12,6,3),(12,4,4),(15,6,3),(16,8,4),(18,9,3)]:
    m=K//j; x=X//j; tot=0; nw=0; atleast=0
    for ones in combinations(range(K),X):
        w=[0]*K
        for p in ones: w[p]=1
        n=N_j(w,j); tot+=n; nw+=1; atleast+=(n>0)
    exact=m*comb(m,x)**j/comb(K,X)
    print(f"shape ({K},{X}) j={j}: mean N over all {nw} words = {tot/nw:.6f}, formula m*C(m,x)^j/C(K,X) = {exact:.6f}; P(N>=1) = {atleast/nw:.4f}")
# Monte Carlo at large m
random.seed(20260929)
print("\nMonte Carlo, uniformly random words of shape (jm, jx):")
for j in (3,4):
    for rho_lab,(xm_num,xm_den) in (("1/2",(1,2)),("~log_3 2 = 12/19",(12,19))):
        for m in (19,38,76,152,304,608,1216):
            x=m*xm_num//xm_den
            if x==0 or x==m: continue
            K=j*m; X=j*x; T=4000 if m<=304 else 1500
            hits=0; s=0
            for _ in range(T):
                w=[1]*X+[0]*(K-X); random.shuffle(w)
                n=N_j(w,j); hits+=(n>0); s+=n
            rho=x/m; asym=sqrt(j)*m**((3-j)/2)*(2*pi*rho*(1-rho))**(-(j-1)/2)
            print(f" j={j} rho={rho_lab} m={m}: P(N>=1)={hits/T:.4f}  E[N]={s/T:.4f} (exact {m*comb(m,x)**j/comb(K,X):.4f}, asym {asym:.4f})")
