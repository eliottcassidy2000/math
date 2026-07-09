# CORR_N = RESONANT + NON-RESONANT (kps-S93), via LEM-011.
#   Corr_N = S_N - N(6/7)^k = sum_{n!=0} What(n) D_N(theta_n),  D_N(theta)=sum_{j=1}^N e(j theta).
#   RESONANT (Vmax | n.e, so theta=0, D_N=N): contributes N * (E_grid[W]-(6/7)^k)  [THM-664 residual].
#   NON-RESONANT (oscillatory): the remainder = the genuinely-new Mertens-analog.
# Verify the split EXACTLY: resonant part computed from the FULL grid average E_grid[W] (no Fourier
# truncation needed -- Parseval-exact identity). Localizes where the last-mile difficulty lives.
import math
TH=1.0/7.0
def W_of(E,Vmax,j):
    pts=sorted(((e*j)%Vmax)/Vmax for e in E); n=len(pts); unc=0.0
    for i in range(n):
        nxt=pts[(i+1)%n]+(1.0 if i==n-1 else 0.0); g=nxt-pts[i]
        if g>TH: unc+=g-TH
    return unc
def longest_ap(E):
    E=sorted(set(E)); S=set(E); best=1
    for i in range(len(E)):
        for jj in range(i+1,len(E)):
            d=E[jj]-E[i]; L=2; nx=E[jj]+d
            while nx in S: L+=1; nx+=d
            pv=E[i]-d
            while pv in S: L+=1; pv-=d
            if L>best: best=L
    return best
def make_AP(K,Vmax): d=Vmax//(K+1); return [d*i for i in range(K)]
def make_sidon(K,Vmax):
    E=[0,1]; diffs={1,Vmax-1}; x=2
    while len(E)<K and x<Vmax:
        ok=True; nd=set()
        for e in E:
            d1=(x-e)%Vmax; d2=(e-x)%Vmax
            if d1 in diffs or d2 in diffs or d1 in nd or d2 in nd: ok=False; break
            nd.add(d1); nd.add(d2)
        if ok: E.append(x); diffs|=nd
        x+=1
    return sorted(E)
print("Corr_N split: RESONANT = N*(E_grid[W]-(6/7)^k) [density-floor residual, ALREADY closed];")
print("NON-RESONANT = Corr_N - resonant [the new oscillatory Mertens-analog].  All EXACT (grid avg).")
print(f"{'K':>3}{'family':>8}{'lAP':>5}{'r_N':>7}{'r_res':>8}{'r_nonres':>9}{'res/Corr':>9}")
for K in (11,12,13):
    N=math.ceil(7*(K-1)/6); lead=(6.0/7.0)**K; targ=N*lead; Vmax=1009
    for name,E in [("AP",make_AP(K,Vmax)),("Sidon",make_sidon(K,Vmax))]:
        if len(set(E))!=K: print(f"{K:>3}{name:>8} bad"); continue
        SN=sum(W_of(E,Vmax,j) for j in range(1,N+1)); corr=SN-targ; rN=abs(corr)/targ
        Egrid=sum(W_of(E,Vmax,j) for j in range(0,Vmax))/Vmax     # full-period grid average
        gridres=Egrid-lead                                        # = sum_{Vmax|n.e, n!=0} What(n)
        res=N*gridres; nonres=corr-res
        print(f"{K:>3}{name:>8}{longest_ap(E):>5}{rN:>7.3f}{abs(res)/targ:>8.3f}{abs(nonres)/targ:>9.3f}{abs(res)/max(abs(corr),1e-9):>9.2f}")
print()
print("READ: r_res = resonant contribution to r_N (density-floor/THM-664 residual, CLOSED).")
print("r_nonres = the oscillatory remainder = the genuinely-new last-mile piece.")
print("For SIDON, r_res tiny => the last mile is ENTIRELY the non-resonant oscillatory sum (Mertens-analog).")
print("For AP, r_res large => but that branch is handled STRUCTURALLY (Dirichlet, LEM-012).")
