# The GMC<->LRC structural bridge, made concrete.
# LRC good-set G'~v is a UNION OF ARCS on the circle. Its v-grid autocorrelation discrepancy is
#    disc_v = (1/v) sum_j A(j/v) - |G|^2  =  sum_{m!=0} |c_{mv}|^2   (Wiener-Khinchin / Parseval).
# CLAIM: the Fourier coeff of a union of arcs is an EXPONENTIAL SUM OVER ARC ENDPOINTS / m --
# EXACTLY the GMC reconstruction "jumps": c_m = (1/(2 pi i m)) sum_j (+/-) e^{-2pi i m x_j},
# x_j = arc endpoints (landings), +/- = jump orientation, 1/m = soft-Weyl / fold decay.
import numpy as np
# union of arcs [a_j,b_j) in [0,1)
arcs=[(0.05,0.23),(0.31,0.52),(0.66,0.70),(0.80,0.97)]
endpoints=[]  # (x, sign): +1 at left end (indicator jumps up), -1 at right end (jumps down)
for a,b in arcs: endpoints += [(a,+1),(b,-1)]
G=sum(b-a for a,b in arcs)
def chat(m):  # exact Fourier coeff c_m = int_G e^{-2pi i m x} dx
    if m==0: return G
    return sum(s*np.exp(-2j*np.pi*m*x) for x,s in endpoints)/(2j*np.pi*m)
# (A) verify the jump/exponential-sum form: c_m * (2 pi i m) = S_m := sum_j sign_j e^{-2pi i m x_j}
for m in [1,3,7,20]:
    Sm=sum(s*np.exp(-2j*np.pi*m*x) for x,s in endpoints)
    print(f"  m={m:2d}: |c_m|={abs(chat(m)):.5f}, |S_m|={abs(Sm):.4f}, |c_m*2pi i m - S_m|={abs(chat(m)*2j*np.pi*m-Sm):.2e}  (c_m = S_m/(2pi i m))")
# (B) verify disc_v two ways: spatial autocorrelation vs Parseval sum of |c_{mv}|^2
def disc_spatial(v,Ngrid=200000):
    x=np.linspace(0,1,Ngrid,endpoint=False)
    ind=np.zeros(Ngrid)
    for a,b in arcs: ind[(x>=a)&(x<b)]=1.0
    # A(tau)=int ind(x)ind(x+tau)dx via FFT autocorrelation
    F=np.fft.fft(ind); A=np.fft.ifft(F*np.conj(F)).real/Ngrid
    taus=np.arange(v)/v
    idx=(taus*Ngrid).astype(int)
    return A[idx].mean()-G**2
def disc_parseval(v,M=4000):
    return sum(abs(chat(m*v))**2 for m in range(-M,M+1) if m!=0)
for v in [5,7,13]:
    print(f"  v={v}: disc_v spatial={disc_spatial(v):.6f}, Parseval sum|c_mv|^2={disc_parseval(v):.6f}  (equal => Wiener-Khinchin)")
print("\nSTRUCTURE: disc_v = sum_{m!=0} |S_{mv}|^2/(2pi m v)^2, S_m = sum_j sign_j e^{-2pi i m x_j}")
print("  = the GMC reconstruction: S_m is exactly GMC's sum_j beta_j e^{-i m theta_j} (jumps=arc endpoints),")
print("  the 1/m is the fold/soft-Weyl decay, and disc_v is its Parseval energy (positive-definite).")
print("  GMC thin-tie closure used Cesaro |S_m|^2 -> sum|beta_j|^2; LRC needs the 1/m^2-WEIGHTED sum bounded.")
