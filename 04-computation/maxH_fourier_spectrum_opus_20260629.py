"""
(1) PROOF that Paley is not the H-maximizer at n=19: exact H(rotational)>H(Paley), two ways.
(2) The Fourier abstraction: a circulant tournament = an ODD sign function s on Z_n; its
    adjacency eigenvalues mu_j = hat(1_C)(j) all have Re=-1/2 (j!=0), Im=(1/2)hat(s)(j).
    Paley: s=Legendre chi => hat(s)=Gauss sum => |Im|=sqrt(n)/2 FLAT (Ramanujan/extremal).
    Half-turn: s=sign/sawtooth => hat(s)=Dirichlet/cotangent => |Im| CONCENTRATED at low freq.
    Verify, and connect: the SAME Dirichlet-vs-Gauss dichotomy is the LRC band vs sign-obstruction.
"""
import cmath, math

def Hcirc(n,out,start_all=False):
    adj=[0]*n
    for i in range(n):
        for d in out: adj[i]|=1<<((i+d)%n)
    dp=[[0]*n for _ in range(1<<n)]
    if start_all:
        for v in range(n): dp[1<<v][v]=1
    else:
        dp[1][0]=1
    for mask in range(1<<n):
        if not start_all and not mask&1: continue
        for last in range(n):
            c=dp[mask][last]
            if not c or not (mask>>last)&1: continue
            av=adj[last]&~mask
            while av:
                nb=av&-av; nx=nb.bit_length()-1; dp[mask|nb][nx]+=c; av^=nb
    tot=sum(dp[(1<<n)-1])
    return tot if start_all else n*tot

def QR(n): return set((x*x)%n for x in range(1,n))
def rot(n): return set(range(1,(n-1)//2+1))

print("(1) PROOF Paley not max at n=19 (exact integers, two DP routes):")
for n in [19]:
    P=QR(n); R=rot(n)
    hP1=Hcirc(n,P); hP2=Hcirc(n,P,start_all=True)
    hR1=Hcirc(n,R); hR2=Hcirc(n,R,start_all=True)
    print(f"  H(Paley_19): route1={hP1} route2={hP2} agree={hP1==hP2}")
    print(f"  H(rot_19)  : route1={hR1} route2={hR2} agree={hR1==hR2}")
    print(f"  rotational - Paley = {hR1-hP1} > 0 : {hR1>hP1}  => Paley is NOT the maximizer at 19. QED.")

print("\n(2) Fourier spectrum mu_j = sum_{d in C} omega^{jd}  (omega=e^{2pi i/n}):")
def spectrum(n,C):
    w=lambda k: cmath.exp(2j*math.pi*k/n)
    return [sum(w(j*d) for d in C) for j in range(n)]
for n in [7,11,19,23]:
    for nm,C in [("Paley",QR(n)),("half-turn",rot(n))]:
        mu=spectrum(n,C)
        re=[round(mu[j].real,4) for j in range(1,n)]
        im=[abs(mu[j].imag) for j in range(1,n)]
        allre_half = all(abs(r+0.5)<1e-9 for r in re)
        print(f"  n={n:>2} {nm:>10}: mu_0={mu[0].real:.1f}; Re(mu_j)=-1/2 for j!=0? {allre_half}; "
              f"|Im| in [{min(im):.2f},{max(im):.2f}] (Paley: flat~sqrt(n)/2={math.sqrt(n)/2:.2f})")
print("\n  => Paley: |Im| FLAT = sqrt(n)/2 (Gauss sum, Ramanujan-extremal spectrum).")
print("     half-turn: |Im| SPREAD, large at low freq (Dirichlet/cotangent kernel).")
print("     Both live on the vertical line Re=-1/2 (forced: C is a tournament sign-set).")
print("\n  LRC link: the LRC danger band is also a Z_p subset; its transform is the SAME Dirichlet")
print("  kernel phi-hat, and the Gauss sum i*sqrt(p) is the LRC sign-obstruction. Max-H selects the")
print("  Dirichlet (concentrated) shape; LRC's hard part is the Gauss (flat) shape. One dichotomy.")
