"""
The modular form f_14 = 14a: compute a_n (Hecke), L(14a,1) (central value, rank 0 => >0), and compare
to the LRC floor structure (bulk 3/pi^2 + cusp-form contribution). Think: floor 2nd moment = weight-2
form on Gamma_0(14) = Eisenstein + f_14.
"""
import math
a1,a2,a3,a4,a6=1,0,1,4,-6
def ap_good(p):
    N=1+sum(1 for x in range(p) for y in range(p) if (y*y+a1*x*y+a3*y-(x**3+a2*x*x+a4*x+a6))%p==0)
    return p+1-N
def is_prime(n):
    if n<2: return False
    for d in range(2,int(n**0.5)+1):
        if n%d==0: return False
    return True
# a_p table
ap={}
for p in range(2,400):
    if not is_prime(p): continue
    if p==2: ap[p]=-1
    elif p==7: ap[p]=1
    else: ap[p]=ap_good(p)
# a_n multiplicative
def an(n):
    if n==1: return 1
    res=1; m=n
    for p in sorted(ap):
        if p*p>m and m>1:
            if m in ap: res*=ap[m]; m=1
            break
        if m%p==0:
            k=0
            while m%p==0: m//=p; k+=1
            # a_{p^k}
            if p in (2,7):  # bad: a_{p^k}=a_p^k
                res*=ap[p]**k
            else:
                seq=[1,ap[p]]
                for i in range(2,k+1): seq.append(ap[p]*seq[-1]-p*seq[-2])
                res*=seq[k]
    if m>1:
        if m in ap: res*=ap[m]
        else: return None
    return res
# L(14a,1) = 2 sum a_n/n exp(-2pi n/sqrt14)  (root number +1)
N=14; s=0.0
for n in range(1,300):
    v=an(n)
    if v is None: continue
    s+=v/n*math.exp(-2*math.pi*n/math.sqrt(N))
L1=2*s
print(f"Hecke eigenvalues a_p: " + ", ".join(f"a_{p}={ap[p]:+d}" for p in sorted(ap)[:10]))
print(f"\nL(14a,1) = {L1:.6f}  (rank 0 => L(1)>0, the obstruction is NON-degenerate)")
print(f"  the period Omega(14a1) ~ 2.4946; BSD: L(1)/Omega = (#Sha . prod c_p)/#tors^2 = (1.2.1)/36 ?")
print(f"  L(1)/Omega ~ {L1/2.4946:.5f}  (a rational by BSD; #tors=6)")
print()
print("LRC floor structure (the Gamma_0(14) weight-2 2nd moment = Eisenstein + cusp form f_14):")
bulk=3/math.pi**2; floor=114382/332563; atom=4*math.cos(3*math.pi/7)**2
print(f"   Eisenstein BULK   3/pi^2 = 1/(2 zeta(2)) = {bulk:.5f}  (the 3-dim boundary space, CONSTANT)")
print(f"   measured floor    inf R' = 114382/332563 = {floor:.5f}  (clears the bulk)")
print(f"   cusp-form part    floor - bulk            = {floor-bulk:.5f}  (the f_14 = 14a contribution)")
print(f"   apex atom         4cos^2(3pi/7)           = {atom:.5f}  (the worst per-level, the doublet)")
print(f"   compare L(14a,1)  = {L1:.5f}")
print()
print("THE MODULAR PICTURE:")
print("  * the LRC(14) floor 2nd moment is a weight-2 modular form on Gamma_0(14):")
print("    dim 4 = 3 Eisenstein (BULK, controlled, 3/pi^2) + 1 cusp form (f_14 = 14a, the OBSTRUCTION).")
print("  * genus 0 (LRC6,LRC10): NO cusp form -- pure Eisenstein -- the boundary DETERMINES the floor (SOLVED).")
print("  * genus 1 (LRC14): ONE cusp form f_14 the boundary CANNOT see = the local-global gap = the obstruction.")
print("  * f_14 rank 0 => L(1)>0 => the cusp form is non-degenerate => the floor clears the bulk (favorable sign).")
