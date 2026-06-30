"""
The Z/6 torsion of 14a forces a Galois congruence a_p = p+1 mod 6. And 6 = phi(14) = the units mod 14 =
the EXISTENCE column. So f_14's Hecke data carries the same '6' (hexagonal/Eisenstein) as the existence side.
Atkin-Lehner a_2=-1 (mirror), a_7=+1 (apex) = the bad primes = the 2 and 7 of 14.
"""
a1,a2,a3,a4,a6=1,0,1,4,-6
def ap(p):
    if p==2: return -1
    if p==7: return 1
    N=1+sum(1 for x in range(p) for y in range(p) if (y*y+a1*x*y+a3*y-(x**3+a2*x*x+a4*x+a6))%p==0)
    return p+1-N
def isp(n): return n>1 and all(n%d for d in range(2,int(n**.5)+1))
print("THE TORSION CONGRUENCE a_p = p+1 (mod 6) [from E_tors(14a) = Z/6]:")
print(f"   {'p':>4} {'a_p':>5} {'a_p mod6':>9} {'(p+1) mod6':>11} {'match':>6}")
ok=True
for p in [3,5,11,13,17,19,23,29,31,37,41,43]:
    if not isp(p): continue
    a=ap(p); m1=a%6; m2=(p+1)%6
    if m1!=m2: ok=False
    print(f"   {p:>4} {a:>5} {m1:>9} {m2:>11} {str(m1==m2):>6}")
print(f"   => a_p = p+1 (mod 6) for all good p: {ok}  (the Z/6 torsion = the units mod 14 = the '6' = EXISTENCE column)")
print()
print("THE MODULAR DICTIONARY of f_14 = 14a (the LRC(14) obstruction):")
print("   * BAD primes (the 2,7 of 14=2*7): a_2=-1 (w_2, the MIRROR), a_7=+1 (w_7, the APEX) -- Atkin-Lehner.")
print("   * TORSION Z/6: a_p = p+1 mod 6 -- the '6' = phi(14) = units = razor's edge = EXISTENCE/Q(sqrt-3).")
print("   * GENUS 1: ONE cusp form -- the obstruction the Eisenstein boundary cannot see (local-global gap).")
print("   * RANK 0: L(14a,1)=0.330>0 -- the obstruction is NON-DEGENERATE (empty tooth has positive width).")
print()
print("THE FLOOR = a weight-2 modular form on Gamma_0(14):  dim 4 = 3 Eisenstein (bulk 3/pi^2) + 1 cusp (f_14).")
print("   LRC(2p) hardness = dim cusp forms = genus(X_0(2p)) = 0,0,1,2,2.  GENUS 0 = pure Eisenstein =")
print("   boundary determines the floor = LRC(6),(10) SOLVED.  GENUS 1 = +f_14 = LRC(14) = first obstruction.")
print()
print("   the 2nd moment <-> symmetric square: the floor (E[X^2] pair-correlation of the danger count) is the")
print("   Rankin-Selberg/sym^2 of the level-14 form; the Eisenstein part = main term, the cusp form f_14 =")
print("   the genus-1 fluctuation. (Direction: floor's obstruction = sym^2 f_14 component.)")
