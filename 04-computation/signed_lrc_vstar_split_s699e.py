"""DECISIVE: the signed structure SPLITS the n=14 worry-set. V*=(1..11,13,24) is an exact floor row
(M=1/14, repo Res_27 certificate) WITH shell-partner 3+24=27=C — refuting 'tight⟹no shell-partner'.
AP has none (max sum 25<27). The V* shell-partner comes from doubling 24=2·12 (prime-2 × C=27=3³).
opus-2026-06-06-S699e."""
from itertools import combinations
def dist_f(x):
    x-=int(x)
    if x<0:x+=1
    return min(x,1-x)
def M_float(V):
    ds=set(V)
    for a,b in combinations(V,2): ds.add(a+b); ds.add(abs(a-b))
    ds.discard(0); best=-1.0
    for d in ds:
        for m in range(d):
            v=min(dist_f(x*m/d) for x in V)
            if v>best: best=v
    return best
def shell_partners(V,C): return [(V[i],V[j]) for i in range(len(V)) for j in range(i+1,len(V)) if (V[i]+V[j])%C==0]
def main():
    n=14; C=2*n-1
    AP=tuple(range(1,14))
    Vstar=(1,2,3,4,5,6,7,8,9,10,11,13,24)
    AP2=tuple(2*v for v in range(1,14))   # 2*AP (nonprimitive)
    for name,V in [("AP",AP),("V*",Vstar),("2*AP",AP2)]:
        m=M_float(V); sp=shell_partners(V,C)
        print(f"{name:5s} V={V}")
        print(f"      M≈{m:.6f} (1/14={1/14:.6f}, tight={abs(m-1/14)<1e-6}); shell-partners (sum≡0 mod {C}) = {sp if sp else 'NONE'}")
    print()
    print("V* shell-partner 3+24=27=C ; note 24 = 2*12 (DOUBLING) — V* replaces AP's 12 by 24.")
    print("So the worry-set SPLITS by signed structure:")
    print("  AP-type (n≤7 all tight rows, AP at all n): NO shell-partner (max sum 2n-3 < 2n-1) — no signed zero-clock.")
    print("  V*-type (first at n=14): HAS shell-partner via the doubling 24=2·12 hitting C=27=3³ — a signed zero-clock.")
    print("The signed LRC is thus STRICTLY FINER than M: it separates the n=14 floor {AP, V*} that M cannot.")
    print("n=14 is exactly where C=2n-1=27 first admits a doubled-speed shell-partner (prime-2 × prime-3).")
if __name__=='__main__': main()
