"""Perfect-power clocks: |2^K - q^X| = m^r with r>=2 (or =1), q odd, 1<=X<=K.
Each such (q,K,X) makes the mixed shape (K,X) free for the 2-adic map y->y/2,(qy+d)/2 with d=+-m^r
(THM-4484 shift criterion), and each is a Pillai / Fermat-Catalan-type identity 2^K -+ m^r = q^X."""
import gmpy2, sys
from gmpy2 import mpz
qmax=int(sys.argv[1]) if len(sys.argv)>1 else 201
Kmax=int(sys.argv[2]) if len(sys.argv)>2 else 400
FC=[ "1^m+2^3=3^2","2^5+7^2=3^4","7^3+13^2=2^9","2^7+17^3=71^2","3^5+11^4=122^2",
     "17^7+76271^3=21063928^2","1414^3+2213459^2=65^7","9262^3+15312283^2=113^7",
     "43^8+96222^3=30042907^2","33^8+1549034^2=15613^3"]
def power_data(n):
    """largest r with n = m^r, r>=2; returns (m,r) or None"""
    n=mpz(n)
    if n<=1: return None
    if not gmpy2.is_power(n): return None
    for r in range(64,1,-1):
        m,exact=gmpy2.iroot(n,r)
        if exact: return (int(m),r)
    return None
found=[]
for q in range(3,qmax+1,2):
    P=[mpz(1)]
    for X in range(1,Kmax+1): P.append(P[-1]*q)
    for K in range(1,Kmax+1):
        two=mpz(1)<<K
        for X in range(1,K+1):
            g=two-P[X]
            a=abs(g)
            if a==1:
                found.append((q,K,X,int(g),1,1)); continue
            pd=power_data(a)
            if pd: found.append((q,K,X,int(g) if a<10**40 else ('sign',1 if g>0 else -1),pd[0],pd[1]))
print("perfect-power clocks |2^K - q^X| = m^r (r>=2) or =1, q odd <=%d, K<=%d:"%(qmax,Kmax))
for f in found:
    q,K,X,g,m,r=f
    typ='mixed' if 1<=X<=K-1 else 'all-odd'
    print(" q=%d (K,X)=(%d,%d) gap=%s = %s%d^%d  [%s; density X/K=%.4f vs log_q 2=%.4f]"%(q,K,X,g,'-' if (isinstance(g,int) and g<0) or (isinstance(g,tuple) and g[1]<0) else '',m,r,typ,X/K,1/ (float(gmpy2.log(q))/float(gmpy2.log(2)))))
print("\nKnown Fermat-Catalan solutions (10):"); print("\n".join(" "+s for s in FC))
