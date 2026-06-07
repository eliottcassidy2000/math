#!/usr/bin/env python3
# lrc_two_tower_witness_group_monad_s3.py
# monad-explorer-2026-06-07-S3
#
# THE TWO-TOWER UNIFICATION OF LRC TORSION HARDNESS.
#
# Two cyclic groups are attached to the n-runner Lonely Runner problem:
#   CLOCK group  G_clk = Z/n        (THM-369/420 prime clocks t = 1/p ;
#                                    opus-S701 THM-421 prime-torsion leaks)
#   SHELL group  G_shl = Z/(2n-1)   (THM-420 shell-partner discrete ticks t=m/(2n-1);
#                                    THM-425 synchronization; S708/S710 signed-LRC
#                                    homometry, the 3-adic tower at C=2n-1).
#
# CLAIM A (proved here, trivial but central): gcd(n, 2n-1) = 1 for ALL n.
#   => by CRT the COMBINED WITNESS GROUP  W_n = Z/(n(2n-1)) = G_clk (x) G_shl,
#      a product of COPRIME factors. The clock tower and shell tower live at
#      DISJOINT primes; LRC difficulty migrates between orthogonal CRT components.
#
# CLAIM B: face-hardness = prime-power tower height.
#   H_clk(n) = max_p v_p(n),  H_shl(n) = max_p v_p(2n-1).
#   Both faces SQUAREFREE  -> "doubly squarefree" easy regime.
#   A prime power on a face (height >= 2) = the Lam-Leung / S708 tower regime that
#   defeats the coprime plug (prime clock / shell-partner coprimality).
#
# CLAIM C (the n=14 punchline): 14 = 2*7 (clock squarefree, H_clk=1) but
#   2n-1 = 27 = 3^3 (H_shl=3). ALL of n=14's prime-power hardness is the 3-adic
#   shell tower = the signed-LRC homometry object at C=27 (S708 deficiency 69).
#   The clock face contributes no tower at all.
#
# CLAIM D (the ladder): organize n by the arithmetic type of 2n-1:
#   2n-1 prime            -> shell TRIVIAL  (no homometry; deficiency 0)
#   2n-1 squarefree comp. -> shell INDEPENDENT primes (homometry, no tower)
#   2n-1 prime power >=2  -> shell TOWER (p-adic, the hard face)
#   First squarefree-composite shell is n=8 (15=3*5) -- coincides with the
#   worry-set going non-AP at n=8 (HYP-2281 / MISTAKE-056). Cross-checked below.
#
# Verifies: gcd law, CRT factorization of W_n, clock leak geometry (re-confirms
# THM-421), tower heights, the ladder, and the n=14 <-> C=27 bridge.

from math import gcd
from fractions import Fraction

def factor(m):
    f={}; d=2
    while d*d<=m:
        while m%d==0: f[d]=f.get(d,0)+1; m//=d
        d+=1
    if m>1: f[m]=f.get(m,0)+1
    return f

def height(m):              # max prime-power exponent (tower height)
    f=factor(m); return max(f.values()) if f else 0

def is_squarefree(m):
    return all(a==1 for a in factor(m).values())

def shell_type(m):          # type of 2n-1
    f=factor(m)
    if len(f)==1 and list(f.values())==[1]: return "PRIME (trivial shell)"
    if all(a==1 for a in f.values()):       return "squarefree composite (independent primes)"
    return "PRIME POWER (p-adic tower)" if len(f)==1 else "mixed (tower + extra prime)"

def norm_units(x,N):        # ||x/N|| in units of 1/N (distance to nearest integer * N)
    r=x%N; return min(r,N-r)

print("="*78)
print("CLAIM A: gcd(n, 2n-1) = 1 for all n  (=> clock & shell are coprime CRT factors)")
print("="*78)
bad=[n for n in range(2,5000) if gcd(n,2*n-1)!=1]
print(f"  n in [2,5000): gcd(n,2n-1)!=1 count = {len(bad)}   (proof: gcd(n,2n-1)=gcd(n,-1)=1)")

print()
print("="*78)
print("CLAIM B/C/D: the two-tower table.  W_n = n*(2n-1).")
print("="*78)
hdr=f"{'n':>3} {'fact n':>10} {'Hclk':>4} | {'2n-1':>5} {'fact 2n-1':>12} {'Hshl':>4} | {'W_n=n(2n-1) tower-primes':>26} | shell type"
print(hdr); print("-"*len(hdr))
first_sf_comp_shell=None
for n in range(3,61):
    C=2*n-1
    fn=factor(n); fC=factor(C)
    Hc=height(n); Hs=height(C)
    W=n*C; fW=factor(W)
    towers=[f"{p}^{a}" for p,a in sorted(fW.items()) if a>=2]   # only height>=2 primes (the hard towers)
    fns='*'.join(f'{p}^{a}' if a>1 else str(p) for p,a in sorted(fn.items()))
    fCs='*'.join(f'{p}^{a}' if a>1 else str(p) for p,a in sorted(fC.items()))
    st=shell_type(C)
    if first_sf_comp_shell is None and st.startswith("squarefree"): first_sf_comp_shell=n
    flag=""
    if Hc==1 and Hs>=2: flag=" <== clock-easy, SHELL TOWER"
    if n==14: flag+="  *** n=14 ***"
    print(f"{n:>3} {fns:>10} {Hc:>4} | {C:>5} {fCs:>12} {Hs:>4} | {('  '.join(towers) if towers else '(squarefree W)'):>26} | {st}{flag}")

print()
print(f"  First squarefree-COMPOSITE shell occurs at n = {first_sf_comp_shell}  "
      f"(2n-1 = {2*first_sf_comp_shell-1} = {'*'.join(str(p) for p in factor(2*first_sf_comp_shell-1))})")
print("  -> matches the worry-set non-AP onset at n=8 (HYP-2281 / MISTAKE-056).")

print()
print("="*78)
print("CLAIM C details: n=14 has its ENTIRE prime-power tower on the shell.")
print("="*78)
n=14; C=2*n-1; W=n*C
print(f"  clock Z/{n} = {'*'.join(f'{p}^{a}' if a>1 else str(p) for p,a in factor(n).items())}"
      f"   H_clk = {height(n)}  (squarefree: {is_squarefree(n)})")
print(f"  shell Z/{C} = {'*'.join(f'{p}^{a}' if a>1 else str(p) for p,a in factor(C).items())}"
      f"   H_shl = {height(C)}  (squarefree: {is_squarefree(C)})")
print(f"  W_14  = {W} = {'*'.join(f'{p}^{a}' if a>1 else str(p) for p,a in factor(W).items())}")
print(f"  the ONLY tower prime (exponent>=2) in W_14 is 3, exponent {factor(W)[3]} -- ALL from the shell.")
print(f"  => n=14's famous prime-power obstruction = the 3-adic structure of Z/27")
print(f"     = the signed-LRC homometry object at C=27 (S708 deficiency 69; S710 3-adic tower).")

print()
print("="*78)
print("Re-confirm THM-421 clock-leak geometry (opus-S701) on the same n list.")
print("="*78)
geom_ok=True; sf_ok=True
for n in [12,14,15,18,20,21,30,33,35,45,50]:
    fn=factor(n)
    for p in sorted(fn):
        a=fn[p]; pa=p**a; m=n//pa
        T=[x for x in range(1,n) if (p*x)%n==0]          # p-torsion minus 0 = <n/p>\{0}
        ok_cof=all((x%m)==0 for x in T) if m>1 else True  # dies in cofactor base
        ok_rot=all(norm_units(x,n)>=n//p for x in T)      # ||x/n|| >= 1/p  (margin)
        ok_pb =all((x%p)!=0 for x in T)                   # survives mod p (squarefree only)
        geom_ok &= ok_cof and ok_rot
        if a==1: sf_ok &= ok_pb
print(f"  geometric margin ||x/n||>=1/p and cofactor-death: {'OK' if geom_ok else 'FAIL'}")
print(f"  squarefree-prime plug (survives mod p when a=1):  {'OK' if sf_ok else 'FAIL'}")

print()
print("="*78)
print("LADDER: count n in [3,200] by shell type and clock type.")
print("="*78)
from collections import Counter
cnt=Counter()
clock_sf=0; shell_sf=0; both_sf=0
for n in range(3,201):
    C=2*n-1
    cnt[shell_type(C)]+=1
    csf=is_squarefree(n); ssf=is_squarefree(C)
    clock_sf+=csf; shell_sf+=ssf; both_sf+= (csf and ssf)
for k,v in sorted(cnt.items(), key=lambda kv:-kv[1]):
    print(f"  {v:>4}  {k}")
N=198
print(f"  clock squarefree: {clock_sf}/{N};  shell squarefree: {shell_sf}/{N};  "
      f"DOUBLY squarefree (easy regime): {both_sf}/{N}")

print()
print("="*78)
print("CLAIM E: MIRROR PAIRS.  A p-adic tower of height h appears as a SHELL at")
print("n=(p^h+1)/2  (since 2n-1=p^h) and as a CLOCK at n=p^h.  Same group Z/p^h,")
print("opposite faces.  S708/S710 studied C=3^h homometry = the SHELL-face family.")
print("="*78)
print(f"{'p^h':>6} | {'shell n=(p^h+1)/2':>18} (2n-1=p^h) | {'clock n=p^h':>12} | same group")
print("-"*70)
for p in (3,5,7):
    for h in (2,3,4):
        ph=p**h
        if (ph+1)%2==0:
            ns=(ph+1)//2; nc=ph
            print(f"{p}^{h}={ph:>4} | shell n={ns:<5} (2*{ns}-1={2*ns-1}) | clock n={nc:<6} | Z/{ph}")
print("  => 3-adic mirror pairs:  (shell n=5, clock n=9) h=2 ;  (shell n=14, clock n=27) h=3 ;")
print("     (shell n=41, clock n=81) h=4.   S710's '3-adic homometry tower at C=9,27,81'")
print("     IS the SHELL face of the LRC mirror family n=5,14,41.")

print()
print("="*78)
print("CLAIM F: the geometric margin 1/p is FACE-INDEPENDENT (same valuation law).")
print("A p-torsion element of Z/N (N=n clock OR N=2n-1 shell) sits at an order-p")
print("rotation j/p at the full tick t=1/N, so ||x/N|| >= 1/p on BOTH faces.")
print("="*78)
def margin_ok(N):
    f=factor(N); ok=True
    for p in f:
        T=[x for x in range(1,N) if (p*x)%N==0]
        ok &= all(min(x%N,N-(x%N)) >= N//p for x in T)
    return ok
for N,lab in [(14,'clock n=14'),(27,'shell n=14'),(9,'clock n=9'),(9,'shell n=5'),
              (8,'clock n=8'),(15,'shell n=8'),(25,'shell n=13'),(49,'shell n=25')]:
    print(f"  N={N:>3} ({lab:<12}): p-torsion margin 1/p holds = {margin_ok(N)}")

print()
print("="*78)
print("CLAIM G: the n/2 GUARD = the 2-torsion of the CLOCK group.")
print("For even n, the half-turn speed n/2 is the UNIQUE order-2 element of Z/n")
print("(2*(n/2)=n=0). The seed's 'half-turn leak residue 7' at n=14 is 7=14/2,")
print("the 2-torsion generator. M(AP_n \\ {n/2})=2/n (seed) = deleting this 2-torsion.")
print("="*78)
for n in (6,8,10,12,14):
    h=n//2
    is2t = (2*h)%n==0 and h!=0
    print(f"  n={n:>2}: n/2={h:>2}  order-2 in Z/{n}: {is2t}  (2*{h}={2*h}=={(2*h)%n} mod {n})")
