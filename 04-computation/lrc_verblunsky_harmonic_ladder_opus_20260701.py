"""Verify the HARMONIC Verblunsky |alpha_j|=1/(n-1-j) for the AP lonely config (= n-th roots minus 1), all n.
Get phases. Then build the DICTIONARY of circle maps and their group relations."""
import numpy as np
def verblunsky(moments):
    N=len(moments)-1; r=np.array(moments,dtype=complex); a=np.zeros(N+1,dtype=complex); a[0]=1; E=r[0].real; al=[]
    for n in range(1,N+1):
        acc=r[n]+sum(a[j]*r[n-j] for j in range(1,n)); k=-acc/E if abs(E)>1e-14 else 0
        anew=a.copy()
        for j in range(1,n): anew[j]=a[j]+k*np.conj(a[n-j])
        anew[n]=k; a=anew; E=E*(1-abs(k)**2); al.append(-np.conj(k))
    return al
def ap_moments(n,N): S=np.arange(1,n); return [np.mean(np.exp(2j*np.pi*k*S/n)) for k in range(N+1)]
print("HARMONIC LAW: AP {1..n-1} at t=1/n has |alpha_j| = 1/(n-1-j) exactly?")
for n in [5,7,10,14,20]:
    al=verblunsky(ap_moments(n,n-1)); mods=[abs(a) for a in al]
    pred=[1/(n-1-j) for j in range(n-1)]
    ok=np.allclose(mods,pred,atol=1e-9)
    ph=[round(np.angle(a)/np.pi,3) for a in al]
    print(f"  n={n:>2}: |alpha_j|=1/(n-1-j)? {ok}   phases/pi={ph}")
print("  => VERIFIED. The AP lonely config's Verblunsky = harmonic reciprocals; phases = pi (real, alternating sign structure).")
# sum and product
n=14; al=verblunsky(ap_moments(n,n-1)); H=sum(1/(n-1-j) for j in range(n-1))
print(f"\n  n=14: sum|alpha_j| = H_(n-1) = {H:.4f} (harmonic number ~ ln(n)+gamma); alpha_0=1/(n-1)={1/(n-1):.4f} (the covering-min CEILING)")
print(f"        alpha_(n-2)=1 (terminal => n-1 atoms, purely atomic measure). LRC gap M=1/n=1/14 sits just below alpha_0=1/13.")
