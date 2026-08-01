"""
FC(n)-homogeneous <=> the polynomial moment problem on the (n-1)-simplex:
   L(f^m) = (Dm+n-1)! int_{Delta_{n-1}} phi_D^m dsigma.
Arc vs area: n=2 (interval) -> 1-D arc -> phi=0 (FC(2) true). n=3 (triangle) -> 2-D area can have all
holomorphic moments zero (disk area measure) -> FC(3) borderline. S_3/C_3 form phi=alpha+omega beta+omega^2
gamma kills 2/3 of the moments (3 nmid m). Full counterexample needs a rotationally-balanced image.
"""
import mpmath as mp
mp.mp.dps = 30
w = mp.exp(2j*mp.pi/3)

def tri_int(g):  # int over 2-simplex {a,b>=0, a+b<=1}
    return mp.quad(lambda a: mp.quad(lambda b: g(a, b) if a+b <= 1 else 0, [0, 1-a]), [0, 1])

# S_3/C_3 generator: 3-fold symmetric equilateral image
def phic(a, b):
    c = 1-a-b
    return a + w*b + w*w*c
print("S_3 mechanism: phi=alpha+omega beta+omega^2 gamma on the triangle (equilateral image, 3-fold symmetric):")
for m in range(1, 7):
    I = tri_int(lambda a, b: phic(a, b)**m)
    tag = "~0 (3 nmid m)" if abs(I) < 1e-14 else "NONZERO (3|m)"
    print("  int_Delta phi^%d = %s   |.|=%s   %s" % (m, mp.nstr(I, 5), mp.nstr(abs(I), 4), tag))
print("=> 2/3 of moments vanish from the 2-D area-enclosing image -- impossible for a 1-D arc (FC(2)).")

# disk area measure: the clean all-moments-zero target for a full FC(3) counterexample
print("\ndisk area measure int_{|z|<=1} z^m dA (the all-moments-zero target):")
for m in [1, 2, 3, 6]:
    val = mp.quad(lambda rr: mp.quad(lambda th: (rr*mp.exp(1j*th))**m*rr, [0, 2*mp.pi]), [0, 1])
    print("  m=%d: int z^m dA = %s (=0 for all m>=1; a polynomial realizing this pushforward would refute FC(3))" % (m, mp.nstr(abs(val), 4)))
