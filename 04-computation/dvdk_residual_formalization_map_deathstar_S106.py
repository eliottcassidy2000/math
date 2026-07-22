# STATUS CORRECTION — HYP-8935 / MISTAKE-241: these floating checks illustrate
# asymptotics; they do not verify formal-log identities, unramified Hensel
# descent, the rational local/global small-root product, or DvdK1.
import numpy as np
np.set_printoptions(suppress=True, precision=6)

# ---------- Laurent polynomial as dict {exp: coeff} ----------
def laurent_pow_CT(fd, m):
    # CT(f^m) via convolution on a shifted coeff array
    lo = min(fd); hi = max(fd)
    arr = np.zeros(hi-lo+1, dtype=complex)
    for e,c in fd.items(): arr[e-lo]=c
    P = np.array([1.0+0j]); 
    for _ in range(m): P = np.convolve(P,arr)
    # exponent of P index k is k + m*lo ; want exponent 0
    idx = -m*lo
    return P[idx] if 0<=idx<len(P) else 0.0

def to_RM(fd):
    # f = u^{-M} R(u), R ordinary poly (coeffs c0..cd), M = -min exponent
    lo=min(fd); M=-lo; d=max(fd)-lo
    R=np.zeros(d+1,dtype=complex)
    for e,c in fd.items(): R[e-lo]=c
    return R, M, d

def D_m(R,M,m):
    P=np.array([1.0+0j])
    for _ in range(m): P=np.convolve(P,R)
    idx=M*m
    return P[idx] if idx<len(P) else 0.0

# ---------- test faces ----------
faces = {
 'S100 hard {-2,-1,1,2} c=(1,1,-1,-1)': {-2:1,-1:1,1:-1,2:-1},
 'THM2070 dihedral u^2+u+u^-1-u^-2':    {2:1,1:1,-1:1,-2:-1},
 'unique-cycle {-1,1,2} (DvdK-free)':   {-1:1.3,1:0.7,2:-2.1},
 'generic complex {-2,-1,1,3}':         {-2:0.4+0.9j,-1:1.1-0.3j,1:-0.7+0.5j,3:2.0+0.1j},
}
print("="*70)
print("CHECK A: CT(f^m) == D_m=[u^{Mm}]R^m, and the cancellation pattern")
for name,fd in faces.items():
    R,M,d=to_RM(fd)
    cts=[laurent_pow_CT(fd,m) for m in range(1,9)]
    dms=[D_m(R,M,m) for m in range(1,9)]
    ok=np.allclose(cts,dms)
    firstnz=next((m for m in range(1,30) if abs(laurent_pow_CT(fd,m))>1e-9),None)
    print(f"\n{name}\n  M={M} d={d}  CT(f^1..8)={np.round(np.array(cts).real,3)}")
    print(f"  D_m match={ok}   first nonzero CT at m={firstnz}")

print("\n"+"="*70)
print("CHECK B: elementary THM-1550 re-derivation")
print("  claim: log[ (-t r_d)(-1)^N * prod(large roots) ] = -sum_{m>=1} D_m t^m/m")
print("  and    prod(small roots) Pi(t) = (-1)^(N+d+1) r_0 t  <=>  all D_m vanish")
for name,fd in faces.items():
    R,M,d=to_RM(fd); N=d-M; r0=R[0]; rd=R[-1]
    # Phi(X)=X^M - t R(X): coeff array in X, for numeric t
    def phi_roots(t):
        coeff=np.zeros(d+1,dtype=complex)      # index i = X^i
        coeff += -t*R
        coeff[M]+=1.0
        # np.roots wants highest-degree first
        return np.roots(coeff[::-1])
    t=1e-3
    rts=phi_roots(t); rts=rts[np.argsort(np.abs(rts))]
    small=rts[:M]; large=rts[M:]
    Pi=np.prod(small)
    const=(-t*rd)*((-1)**N)*np.prod(large)
    # RHS series
    Dms=[D_m(R,M,m) for m in range(1,40)]
    rhs=-sum(Dms[m-1]*t**m/m for m in range(1,40))
    lhs=np.log(const)
    # Pi vs c t
    c=((-1)**(N+d+1))*r0
    print(f"\n{name}")
    print(f"  log(const)={lhs:.6f}   -sum D_m t^m/m={rhs:.6f}   match={np.isclose(lhs,rhs,atol=1e-6)}")
    print(f"  Pi(t)={Pi:.6e}   c*t={c*t:.6e}   Pi==ct? {np.isclose(Pi,c*t,rtol=1e-3)}  (false unless all D_m vanish)")

print("\n"+"="*70)
print("CHECK C: small-root product via UNRAMIFIED Hensel  Z^M = R(sZ),  t=s^M")
print("  the M unit-roots Z_i(s) in C[[s]] (lift of Z^M=r0); Pi(t)=s^M*prod Z_i = t*prod Z_i")
def hensel_unit_roots(R,M,d,s_val):
    # solve psi(Z)=Z^M - R(s*Z)=0 numerically at s=s_val, pick M roots near the M-th roots of r0
    r0=R[0]
    # psi as polynomial in Z: Z^M - sum_i R[i] s^i Z^i
    deg=max(M,d)
    coeff=np.zeros(deg+1,dtype=complex)
    coeff[M]+=1.0
    for i in range(d+1): coeff[i]-= R[i]*(s_val**i)
    Z=np.roots(coeff[::-1])
    # unit roots: |Z| ~ |r0|^{1/M}, O(1); others blow up ~ s^{-...}
    target=abs(r0)**(1.0/M)
    Z=Z[np.argsort(np.abs(np.abs(Z)-target))]
    return Z[:M]
for name,fd in faces.items():
    R,M,d=to_RM(fd); N=d-M
    t=1e-3; s=t**(1.0/M)
    Zi=hensel_unit_roots(R,M,d,s)
    Pi_hensel = (s**M)*np.prod(Zi)
    # direct
    coeff=np.zeros(d+1,dtype=complex); coeff+=-t*R; coeff[M]+=1.0
    rts=np.roots(coeff[::-1]); rts=rts[np.argsort(np.abs(rts))]
    Pi_direct=np.prod(rts[:M])
    print(f"\n{name}\n  Pi_direct={Pi_direct:.6e}   Pi_hensel={Pi_hensel:.6e}   match={np.isclose(Pi_direct,Pi_hensel,rtol=1e-2)}")

print("\n"+"="*70)
print("CHECK D: Vieta C_Phi=prod(all roots)=(-1)^d r0/rd (val_t=0);  Pi(t) unramified (int powers)")
for name,fd in faces.items():
    R,M,d=to_RM(fd); N=d-M; r0=R[0]; rd=R[-1]
    t=1e-3
    coeff=np.zeros(d+1,dtype=complex); coeff+=-t*R; coeff[M]+=1.0
    rts=np.roots(coeff[::-1])
    Cphi=np.prod(rts); vieta=((-1)**d)*r0/rd
    # Pi(t) unramified test: Pi(t)/(c t) should -> exp(sum D_m t^m/m), a power series in t (integer powers).
    # test: Pi(t) for t and for t*e^{2pi i/M} (rotate t) : if unramified, Pi rotates by same factor to leading order
    c=((-1)**(N+d+1))*r0
    def Pi_of(tt):
        cf=np.zeros(d+1,dtype=complex); cf+=-tt*R; cf[M]+=1.0
        rr=np.roots(cf[::-1]); rr=rr[np.argsort(np.abs(rr))]; return np.prod(rr[:M])
    # ratio Pi(t)/(c t) at two small t; should agree to O(t) (a genuine power series in t)
    r1=Pi_of(1e-4)/(c*1e-4); r2=Pi_of(2e-4)/(c*2e-4)
    print(f"\n{name}\n  C_Phi={Cphi:.5f}  (-1)^d r0/rd={vieta:.5f}  Vieta match={np.isclose(Cphi,vieta,rtol=1e-3)}")
    print(f"  Pi/(ct) at t=1e-4:{r1:.5f}  t=2e-4:{r2:.5f}  ->1 as t->0 (unramified power series in t): {np.isclose(r1,1,atol=2e-3) and np.isclose(r2,1,atol=3e-3)}")
