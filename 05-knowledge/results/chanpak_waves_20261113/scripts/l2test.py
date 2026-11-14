"""Wave 24, P55: L^2 local-limit test. For each dyadic window of odd N, model(N) = mean * f9(N mod 9) * f5(N mod 5)
* f7(N mod 7), with the class factors measured in that window. V = mean((r/model - 1)^2); Poisson noise alone gives
about 1/model. Reports V, V*mean (overdispersion), and the share of N with r < model/2."""
import numpy as np, sys
for f in sys.argv[1:]:
    r=np.fromfile(f,dtype=np.uint16).astype(np.float64)   # r[N>>1]
    print(f)
    print(" window   mean     V        V*mean   share(r<model/2)  max|r/model-1|")
    K=int(round(np.log2(len(r))))+1     # r has 2^(K-1) entries (odd N < 2^K)
    for k in range(14,K):
        a,b=(1<<k)+1,(1<<(k+1))-1
        N=np.arange(a,b+1,2); x=r[N>>1]
        model=np.full_like(x,x.mean())
        for m in (9,5,7):
            cls=N%m; fac=np.array([x[cls==c].mean() for c in range(m)])/x.mean(); model*=fac[cls]
        z=x/model-1; V=(z**2).mean()
        print(" 2^%-2d %8.2f %9.5f %8.3f %12.6f %10.3f"%(k,x.mean(),V,V*x.mean(),(x<model/2).mean(),np.abs(z).max()))
