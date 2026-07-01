"""VERBLUNSKY / OPUC lens on the LRC. A measure mu on the unit circle <-> Verblunsky coeffs {alpha_n} (Simon).
Runner config at time t for set S = uniform on {e^{2pi i v t}: v in S} (points on a loop). Compute alpha_n via
Levinson-Durbin on the moments m_k = (1/|S|) sum_v e^{2pi i k v t}. Look for LRC structure (tight vs not, the gap)."""
import numpy as np
def verblunsky(moments):
    # moments m[0..N] (m[0]=1, complex); Levinson-Durbin reflection coeffs = -conj(alpha)
    N=len(moments)-1; r=np.array(moments,dtype=complex)
    a=np.zeros(N+1,dtype=complex); a[0]=1.0; E=r[0].real; alphas=[]
    for n in range(1,N+1):
        acc=r[n]+sum(a[j]*r[n-j] for j in range(1,n))
        k=-acc/E if abs(E)>1e-14 else 0
        anew=a.copy()
        for j in range(1,n): anew[j]=a[j]+k*np.conj(a[n-j])
        anew[n]=k; a=anew; E=E*(1-abs(k)**2)
        alphas.append(-np.conj(k))   # Verblunsky alpha_{n-1}
    return alphas
def config_moments(S,t,N):
    S=np.array(S); return [np.mean(np.exp(2j*np.pi*k*S*t)) for k in range(N+1)]
def uniform_roots_moments(n,N):  # uniform on all n-th roots
    return [1.0 if k%n==0 else 0.0 for k in range(N+1)]
n=14
print("Verblunsky coeffs |alpha_n| for LRC configs (runners on the loop):")
# (1) uniform n-th roots (AP lonely config)
vr=verblunsky(uniform_roots_moments(14,14))
print(f"  uniform 14-th roots: |alpha|={[round(abs(a),3) for a in vr]}  (expect 0..0,1 -- atomic/Szego)")
# (2) AP {1..13} runner config at the LONELY time t=1/14
ap=verblunsky(config_moments(list(range(1,14)),1/14,13))
print(f"  AP {{1..13}} at t=1/14 (lonely): |alpha|={[round(abs(a),3) for a in ap]}")
# (3) AP at a GENERIC (non-lonely) time
apg=verblunsky(config_moments(list(range(1,14)),0.137,13))
print(f"  AP {{1..13}} at t=0.137 (generic): |alpha|={[round(abs(a),3) for a in apg]}")
# (4) construction {1..12,182} at its lonely time 14/183
con=verblunsky(config_moments(list(range(1,13))+[182],14/183,13))
print(f"  construction at t=14/183: |alpha|={[round(abs(a),3) for a in con]}")
# (5) a spread/far set -> equidistribution -> Verblunsky should DECAY (Szego)
np.random.seed(0); spread=[1]+sorted(np.random.choice(range(2,400),12,replace=False).tolist())
sp=verblunsky(config_moments(spread,0.234,13))
print(f"  spread set at generic t: |alpha|={[round(abs(a),3) for a in sp]}  (decay = equidistribution/Szego)")
print("\n=> reading: atomic/structured (AP,construction) => |alpha| stay O(1) near the boundary; spread => decay.")
print("   The Szego/Verblunsky decay = the STRUCTURED-vs-SPREAD two-regime (HYP-3790) in OPUC language.")
