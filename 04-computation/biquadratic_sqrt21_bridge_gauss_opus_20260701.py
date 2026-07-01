"""Test the Q(sqrt-3, sqrt-7) bridge conjecture. Classical quadratic Gauss sum g(m)=sum_{a=0}^{m-1} e^{2pi i a^2/m}.
p=3 mod4 => g=i*sqrt(p) (iota-ODD, imaginary); m=1 mod4 => g=sqrt(m) (iota-EVEN, real). Test sqrt21 = g(21) =
the PRODUCT of the two iota-odd certificates i*sqrt3 (Eisenstein) and i*sqrt7 (heptagon)."""
import numpy as np
from math import gcd
def gauss(m): return sum(np.exp(2j*np.pi*(a*a)/m) for a in range(m))
for m in [3,7,21,61,183]:
    g=gauss(m); modp=m%4
    print(f"  g({m}) = {g.real:+.4f}{g.imag:+.4f}i   |g|={abs(g):.4f}  sqrt({m})={np.sqrt(m):.4f}  ({m}=%d mod4 => {'iota-ODD i*sqrt' if modp==3 else 'iota-EVEN sqrt'})"%modp)
g3,g7,g21=gauss(3),gauss(7),gauss(21)
print(f"\n  i*sqrt3 (Eisenstein, iota-odd) * i*sqrt7 (heptagon, iota-odd) = {(g3*g7).real:+.4f}{(g3*g7).imag:+.4f}i")
print(f"  vs g(21)=sqrt21 (iota-EVEN) = {g21.real:+.4f}  => product of two ODD certs = +-EVEN sqrt21? {np.isclose(abs(g3*g7),abs(g21))}")
print(f"  |g3*g7|={abs(g3*g7):.4f} = sqrt(21)={np.sqrt(21):.4f}  (two iota-odd i-Gauss => one iota-even real: cup product odd x odd -> even)")
print("\n  FIELD: Q(sqrt-3, sqrt-7) biquadratic, Gal=Z/2xZ/2, 3 quad subfields Q(sqrt-3),Q(sqrt-7),Q(sqrt21).")
print(f"  LRC-14: N=14=2*7 (heptagon=7 => sqrt-7); Phi6(14)={14*14-14+1}=3*61 (Eisenstein 3 => sqrt-3); |Aut(Paley_7)|=21=3*7 (shared 3).")
print(f"  So sqrt21=sqrt(3*7) BRIDGES the Eisenstein covering (3, from Phi6) and the heptagon (7) = the CROSS term.")
# is 61 also relevant? 61 mod 3, mod 7
print(f"  61 mod3={61%3} (=1 => 61 splits in Q(sqrt-3), Eisenstein); 61 mod7={61%7}; Phi6=3*61 both Eisenstein-friendly.")
