import math, numpy as np, sys
g=math.log(2)/math.log(3)
H=-g*math.log2(g)-(1-g)*math.log2(1-g)
rho=2**(H-1)
print(f"theory: survivors fraction ~ C K^-3/2 rho^K, rho = 2^(H(log_3 2)-1) = 2^-{1-H:.5f} = {rho:.6f}")
def load(fn):
    d={}
    for l in open(fn):
        if l.startswith('#'): continue
        p=l.split(); d[int(p[0])]=int(p[1])
    return d
sets={'descent':'out_desc_K36.txt','join':'out_join_K32.txt','pred1':'out_pred1_K32.txt','pred1+join':'out_pred1_join_K30.txt','branch(max2adic)':sys.argv[1] if len(sys.argv)>1 else 'out_branch_mergeafter_K30.txt'}
for name,fn in sets.items():
    d=load(fn); Ks=[K for K in sorted(d) if K>=14]
    y=np.array([math.log(d[K]/2**K) for K in Ks]); K=np.array(Ks,float)
    # free fit log f = c + b log K + K log r
    A=np.vstack([np.ones_like(K),np.log(K),K]).T; coef,*_=np.linalg.lstsq(A,y,rcond=None)
    # fixed theoretical rho and exponent -3/2: fit only constant -> residual trend
    resid=y-(-1.5*np.log(K)+K*math.log(rho)); C=np.exp(resid)
    print(f"{name:18s} K={Ks[0]}..{Ks[-1]}: free fit rho={math.exp(coef[2]):.5f}, K-exponent {coef[1]:+.2f};  C(K)=f*K^1.5/rho^K at K={Ks[-3]},{Ks[-2]},{Ks[-1]}: {C[-3]:.3f} {C[-2]:.3f} {C[-1]:.3f}")
