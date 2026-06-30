"""
UNIFICATION 2: the apex atom of LRC(2p) = 2+lambda_min(C_p) = 4sin^2(pi/2p) (the minimal odd p-cycle);
the genus of X_0(2p) is the obstruction multiplicity; LRC(2p) hard <=> genus>=1. And the dual minimal
cycles: tournament C_3 (3-cycle, OCF) <-> LRC C_p. Apex primes 3,7 = Mersenne n Heegner n 3mod4.
"""
import math
def genus_X0_2p(p):
    # genus of X_0(2p), 2p squarefree: g = 1 + N/12 - nu2/4 - nu3/3 - ninf/2 ; ninf=#cusps=4 (for 2p)
    N=2*p
    # nu2 = # elliptic pts of order 2: prod over p|N of (1+(-1|p)) ; here primes 2 and p
    def leg(a,q): 
        a%=q
        return 0 if a==0 else (1 if pow(a,(q-1)//2,q)==1 else -1)
    nu2 = (1+ (-1 if 2==2 else 0)) # special at 2: for N even, nu2=0
    nu2 = 0  # 4|N? no; for 2p with p odd, nu2 = prod (1+(-1/q)); at q=2 factor is 0 if 4|N else... use table
    # use known: genus X_0(6)=0,X_0(10)=0,X_0(14)=1,X_0(22)=2,X_0(26)=2
    return {3:0,5:0,7:1,11:2,13:2}.get(p,None)
print("LRC(2p): apex atom = 2+lambda_min(C_p) = 4sin^2(pi/2p)  vs  genus(X_0(2p)) (obstruction multiplicity):")
print(f"   {'p':>3} {'2p':>4} {'apex atom 4sin^2(pi/2p)':>24} {'= 2+lambda_min(C_p)':>20} {'genus':>6} {'status':>14}")
for p in [3,5,7,11,13]:
    atom=4*math.sin(math.pi/(2*p))**2
    lam=2+2*math.cos(math.pi*(p-1)/p)  # 2+lambda_min(C_p)
    g=genus_X0_2p(p)
    status={0:"easy (Eisenstein)",1:"FIRST HARD",2:"harder"}[g]
    print(f"   {p:>3} {2*p:>4} {atom:>24.5f} {lam:>20.5f} {g:>6} {status:>14}")
print("   => atom DECREASES smoothly (~pi^2/p^2); GENUS JUMPS 0->1 at p=7 (X_0(14)). LRC(14)=first genus-1.")
print("   p=3 (LRC6): atom=1, genus 0 -- the C_3 = the tournament's minimal odd cycle (OCF origin).")
print("   p=7 (LRC14): atom=0.198, genus 1 -- the C_7 = ONE cusp form. Apex primes {3,7}=Mersenne n Heegner n 3mod4.\n")
print("THE TWO-COLUMN UNIFICATION (the whole LRC(14) in one duality):")
rows=[("ATOM","the DOUBLET (2-core, min resonance)","the EMPTY TOOTH (comb gap, min hole)"),
      ("regime","OFF-CUSP (proper core)","CUSP (core=Z_7)"),
      ("carries","MEASURE (density rho_j>=0.198)","EXISTENCE (witness D(q*)>=1)"),
      ("value","4cos^2(3pi/7)=0.198","M=n/Phi_6(n)=14/183"),
      ("Heegner field","Q(sqrt-7) (apex/Mersenne)","Q(sqrt-3) (Eisenstein/hexagonal)"),
      ("denominator","the gap (cyclotomic Q(cos2pi/7))","Phi_6(n)=n^2-n+1 (Eisenstein norm)"),
      ("X_0(14) cusp","d=7 (- / apex / w_7=-1)","d=14 (S / Fricke / the comb)"),
      ("modular","f_14 leading coeff (vanishes at cusp)","the comb witness / rank-0 L(1)>0"),
      ("status","PROVED (THM-590, finite)","OPEN (covering-min >= 1/n)"),
      ("razor-thin","the MEASURE -> 0 (vanishes)","the WITNESS robust (exists)"),
      ("tournament","the 3-CYCLE (min odd cycle, OCF)","the empty tooth / max-packing"),]
print(f"   {'':>14} | {'MEASURE column':>38} | {'EXISTENCE column':>38}")
for r in rows:
    print(f"   {r[0]:>14} | {r[1]:>38} | {r[2]:>38}")
print("\n  => ONE duality: (measure/existence) = (doublet/empty-tooth) = (Q(-7)/Q(-3)) = (d=7/d=14 cusp)")
print("     = (proved/open) = (4cos^2/Phi_6). The two Heegner fields ARE the two columns.")
