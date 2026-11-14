import sys, time
from raa import phi, best_d
def raa_profile(bits):
    T=len(bits); k=next(i for i,b in enumerate(bits) if b)
    alpha=int(''.join(str(b) for b in reversed(bits)),2)
    f=(0,2); g=(1<<k,1); maxdev=0; arg=0; jumps={}; prevlog=None
    while k<T-1:
        k+=1; mod=1<<(k+1); am=alpha&(mod-1)
        if (am*g[1]-g[0])%mod==0: f=(2*f[0],2*f[1])
        elif phi(g)<phi(f):
            h=best_d(f,g); f=(2*g[0],2*g[1]); g=h
        else:
            h=best_d(g,f); g=h; f=(2*f[0],2*f[1])
        lg=phi(g).bit_length()
        if prevlog is not None and lg>prevlog: jumps[lg-prevlog]=jumps.get(lg-prevlog,0)+1
        prevlog=lg
        dev=abs(2*lg-(k+1))
        if dev>maxdev: maxdev,arg=dev,k+1
    return g, lg, maxdev, arg, jumps
if __name__ == "__main__":
    T=int(sys.argv[1]); raw=open(sys.argv[2],'rb').read()[:T]; bits=list(raw)
    t0=time.time(); g,lg,md,arg,j=raa_profile(bits)
    print(f"T={T}: final 2-adic complexity log2 phi = {lg} (T/2={T//2}); max |2 log2phi(n) - n| = {md} at n={arg}; jump histogram {dict(sorted(j.items()))}; {time.time()-t0:.1f}s")
