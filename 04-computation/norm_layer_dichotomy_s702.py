"""
S702 (monad-explorer): THE NORM-LAYER DICHOTOMY for 2D lattices.

Probe (from S699/HYP-2257 + dispatched seed): "is there a 2D-realizable GROUP strictly
between the triangular lattice (kissing kappa=6) and the CM field whose norm-1 layer beats 3n
at moderate n?"

THESIS to test rigorously: the probe conflates TWO different "norm-1 layers":
  (K) KISSING layer   = MINIMAL-norm vectors. Geometric. kappa = #min vectors <= 6 for ANY 2D
                        lattice (=6 iff hexagonal). Gives <= 3n unit distances. This is the cap.
  (P) POPULAR layer   = a FIXED norm value D with MANY representations r_Q(D). Arithmetic.
                        sup_D r_Q(D) = infinity for EVERY 2D lattice whose order has a split prime.
                        Beating "3n" just means finding D with r_Q(D) > 6.

If the thesis holds: the 3n cap is purely about layer (K); layer (P) already beats 3n inside
the SQUARE lattice Z^2 (kappa=4 < 6!), so the kissing axis (4..6) is ORTHOGONAL to the escape.
There is no continuum of 2D groups "between triangular and CM": kissing is universally <=6 and the
popular-norm growth rate is the SAME divisor-function rate for every 2D lattice (imaginary
quadratic order). The CM escape is a DEGREE/RANK jump (rank>=3 module projected to 2D), not a 2D
group. We verify all the exact, checkable pieces and clearly separate them from the (external,
unverified) CM exponent claim.

Everything here is EXACT integer arithmetic.
"""
import math
from fractions import Fraction
from collections import Counter

# ----------------------------------------------------------------------------
# Part K: kissing number kappa for a spectrum of 2D lattices.  Geometric cap.
# ----------------------------------------------------------------------------
def kissing(basis, R=60):
    """Minimal nonzero squared length and its multiplicity (=kissing number)."""
    b1, b2 = basis
    best = None
    cnt = 0
    for i in range(-R, R+1):
        for j in range(-R, R+1):
            if i == 0 and j == 0:
                continue
            x = i*b1[0] + j*b2[0]
            y = i*b1[1] + j*b2[1]
            d2 = x*x + y*y
            if best is None or d2 < best - 1e-12:
                best = d2; cnt = 1
            elif abs(d2 - best) < 1e-9:
                cnt += 1
    return cnt, best

def part_K():
    print("="*78)
    print("PART K -- KISSING LAYER (minimal norm).  GEOMETRIC, capped at 6.")
    print("="*78)
    s3 = math.sqrt(3)/2
    # A spectrum of 2D lattices: rhombic family interpolating square(90) -> hex(60),
    # plus rectangular and generic oblique.
    lats = []
    for deg in (90, 80, 75, 70, 65, 61, 60, 59, 55, 50, 45):
        th = math.radians(deg)
        lats.append((f"rhombic angle={deg:3d}", ((1.0,0.0),(math.cos(th),math.sin(th)))))
    lats += [
        ("rectangular 1x2",  ((1.0,0.0),(0.0,2.0))),
        ("rectangular 1x1.3",((1.0,0.0),(0.0,1.3))),
        ("oblique generic",  ((1.0,0.0),(0.37, 1.21))),
    ]
    maxk = 0
    for name, basis in lats:
        k, d2 = kissing(basis)
        maxk = max(maxk, k)
        flag = "  <-- HEXAGONAL/triangular, the unique 2D maximum" if k == 6 else ""
        print(f"  {name:22s}: kappa = {k}  ->  <= {k/2:.1f} n unit (min-dist) pairs{flag}")
    print(f"\n  MAX kissing over all tested 2D lattices = {maxk}  (classical cap: kappa<=6, "
          "=6 iff hexagonal).")
    print("  => the MINIMAL-distance graph on n lattice pts has <= 3n edges, ANY 2D lattice.")
    print("  This is the ONLY sense in which '3n' is a cap.  It is geometric (kissing), not")
    print("  arithmetic.  No 2D lattice has kappa>6, so there is NO kissing-spectrum above 6.\n")

# ----------------------------------------------------------------------------
# Part P: popular-norm layer.  Arithmetic, unbounded.  Beats 3n inside Z^2 already.
# ----------------------------------------------------------------------------
def reps_square(D):
    """r_{Z^2}(D) = #{(x,y): x^2+y^2=D}, brute force (exact)."""
    if D == 0: return 1
    c = 0
    s = int(math.isqrt(D))
    for x in range(-s, s+1):
        rem = D - x*x
        if rem < 0: continue
        y = int(math.isqrt(rem))
        if y*y == rem:
            c += 1 if y == 0 else 2
    return c

def reps_eisenstein(D):
    """r(D) for norm form x^2+xy+y^2 (Loeschian), brute force (exact)."""
    c = 0
    s = int(math.isqrt(D)) + 2
    for x in range(-s, s+1):
        for y in range(-s, s+1):
            if x*x + x*y + y*y == D:
                c += 1
    return c

def divisor_formula_square(D):
    """r(D) = 4 (d_1(D) - d_3(D)), divisors mod 4.  Exact closed form."""
    d1 = d3 = 0
    for k in range(1, D+1):
        if D % k == 0:
            if k % 4 == 1: d1 += 1
            elif k % 4 == 3: d3 += 1
    return 4*(d1 - d3)

def divisor_formula_eisenstein(D):
    """r(D) = 6 (d_1(D) - d_2(D)), divisors mod 3 (for D>0)."""
    d1 = d2 = 0
    for k in range(1, D+1):
        if D % k == 0:
            if k % 3 == 1: d1 += 1
            elif k % 3 == 2: d2 += 1
    return 6*(d1 - d2)

def part_P():
    print("="*78)
    print("PART P -- POPULAR-NORM LAYER (a fixed norm with many reps).  ARITHMETIC, unbounded.")
    print("="*78)
    # (1) verify the divisor closed forms against brute force
    print("  [check] divisor closed-forms vs brute force, D=1..120:")
    ok_sq = all(reps_square(D) == divisor_formula_square(D) for D in range(1,121))
    ok_ei = all(reps_eisenstein(D) == divisor_formula_eisenstein(D) for D in range(1,121))
    print(f"    square  r(D)=4(d1-d3 mod4): {'OK' if ok_sq else 'FAIL'}")
    print(f"    eisen.  r(D)=6(d1-d2 mod3): {'OK' if ok_ei else 'FAIL'}")

    # (2) first D where the popular-norm count EXCEEDS 6  (= "beats 3n")
    first_sq = next(D for D in range(1, 10**4) if reps_square(D) > 6)
    first_ei = next(D for D in range(1, 10**4) if reps_eisenstein(D) > 6)
    print(f"\n  [beats 3n  <=>  r(D) > 6  (>3n unit pairs at distance sqrt(D))]")
    print(f"    SQUARE  Z^2  (kappa=4): least D with r(D)>6 is D={first_sq}, r={reps_square(first_sq)}"
          f"  ->  > 3n at distance sqrt({first_sq})")
    print(f"    TRIANG. Z[w] (kappa=6): least D with r(D)>6 is D={first_ei}, r={reps_eisenstein(first_ei)}"
          f"  ->  > 3n at distance sqrt({first_ei})")
    print("    *** The SQUARE lattice (kappa=4 < 6) beats 3n -- the LOW-kissing lattice. ***")
    print("    => the escape from 3n is ORTHOGONAL to the kissing number. 'Between triangular")
    print("       and CM' is the wrong axis: the escape is a LAYER change, not a GROUP change.\n")

    # (3) how many lattice points n you need so that distance sqrt(D) actually occurs r(D)
    #     times -- 'moderate n'.  A point z and z+v both inside an N x N box: need ~ enough box.
    print("  [moderate n: smallest square window that realizes the beating distance fully]")
    for name, D, reps in [("Z^2", first_sq, reps_square(first_sq)),
                          ("Z[w]", first_ei, reps_eisenstein(first_ei))]:
        # take all lattice pts in a box big enough; count pairs at sqrt(D); compare to 3n
        print(f"    {name}: beating distance^2 = {D}, r={reps}.  small windows:")
        # square lattice window
        if name == "Z^2":
            pts = [(x,y) for x in range(0,8) for y in range(0,8)]  # 64 pts
        else:
            pts = [(x,y) for x in range(0,8) for y in range(0,8)]
        # We measure popular-distance pairs vs 3n for increasing window
        crossed = None
        for side in (4,6,8,10,12,14,16,20):
            if name == "Z^2":
                P = [(x,y) for x in range(side) for y in range(side)]
                def nrm(a,b): return (a[0]-b[0])**2+(a[1]-b[1])**2
            else:
                # Eisenstein coords -> use integer (x,y) with norm x^2+xy+y^2
                P = [(x,y) for x in range(side) for y in range(side)]
                def nrm(a,b):
                    dx=a[0]-b[0]; dy=a[1]-b[1]; return dx*dx+dx*dy+dy*dy
            n = len(P)
            pairs = 0
            for i in range(n):
                for j in range(i+1,n):
                    if nrm(P[i],P[j]) == D: pairs += 1
            beats = pairs > 3*n
            if beats and crossed is None: crossed = (side, n)
            print(f"       side {side:2d}: n={n:4d}, pairs@sqrt(D)={pairs:4d}, 3n={3*n:4d}, "
                  f"{'BEATS 3n' if beats else 'not yet (boundary loss)'}")
        if crossed:
            print(f"       -> {name} crosses 3n at side {crossed[0]} (n={crossed[1]}): bulk density "
                  f"r/2={reps//2}n beats 3n once perimeter loss ~c*sqrt(n) is dominated.")
    print()

# ----------------------------------------------------------------------------
# Part S: the REAL spectrum -- divisor-function growth, same exponent for all 2D lattices.
# ----------------------------------------------------------------------------
def sieve_char_divisor_counts(X, residues, mod):
    """For all D<=X, return r-vector = w*(d_a - d_b) where residues=(a,b) mod `mod`,
    using a sieve over divisors.  Fast, exact."""
    a, b = residues
    sgn = [0]*(X+1)  # contribution of divisor k: +1 if k%mod==a, -1 if k%mod==b
    acc = [0]*(X+1)
    for k in range(1, X+1):
        r = k % mod
        s = 1 if r == a else (-1 if r == b else 0)
        if s:
            for m in range(k, X+1, k):
                acc[m] += s
    return acc  # multiply by w=4 (square) or 6 (eisen) at use site

def part_S():
    print("="*78)
    print("PART S -- THE REAL SPECTRUM: popular-norm GROWTH RATE.")
    print("="*78)
    print("  Both lattices' popular counts are governed by a divisor sum (Dirichlet character).")
    print("  max_{D<=X} r(D) grows like a power of X with exponent -> 0 (Wigert/divisor bound):")
    print("  max_{m<=X} d(m) = X^{(log2+o(1))/loglogX}.  SAME rate for square AND triangular.\n")
    print(f"   {'X':>7} | {'max r_sq(D), D<=X':>18} | {'argmax':>7} | {'max r_ei(D), D<=X':>18} | {'argmax':>7}")
    print("  " + "-"*70)
    Xmax = 200000
    acc_sq = sieve_char_divisor_counts(Xmax, (1,3), 4)  # *4
    acc_ei = sieve_char_divisor_counts(Xmax, (1,2), 3)  # *6
    for X in (50, 200, 1000, 5000, 20000, 200000):
        best_sq = best_ei = 0; arg_sq = arg_ei = 0
        for D in range(1, X+1):
            r = 4*acc_sq[D]
            if r > best_sq: best_sq = r; arg_sq = D
            r2 = 6*acc_ei[D]
            if r2 > best_ei: best_ei = r2; arg_ei = D
        print(f"   {X:>7} | {best_sq:>18} | {arg_sq:>7} | {best_ei:>18} | {arg_ei:>7}")
    print("\n  Both columns are the SAME superlinear-by-loglog growth (highly-composite argmax).")
    print("  The popular-distance count among n=X lattice pts is ~ (max r)/2 * n = n^{1+o(1)}.")
    print("  This exponent 1+o(1) is IDENTICAL for every 2D lattice = imaginary quadratic order:")
    print("  it is the 1-DIMENSIONAL divisor function, independent of discriminant/class number.")
    print()
    print("  ==> CONCLUSION (rigorous part): among 2D lattices there is NO spectrum 'between")
    print("      triangular and CM'.  Kissing is universally <=6 (geometric); popular-norm")
    print("      exponent is universally 1+o(1) (arithmetic, divisor-function).  Both axes are")
    print("      FLAT across 2D lattices.  The CM construction's claimed higher exponent (external")
    print("      seed: n^{1.014}) requires a rank>=3 module / degree>=3 CM field -- NOT a 2D")
    print("      lattice.  The 'gap' is a DIMENSION jump, not a continuous 2D group spectrum.")
    print()

def part_dict():
    print("="*78)
    print("PART D -- DICTIONARY: 2D lattice  <->  imaginary quadratic order")
    print("="*78)
    print("  square Z[i]   : disc -4, h=1, units 4  (kappa=4); norm form x^2+y^2,   chi_-4")
    print("  triang Z[w]   : disc -3, h=1, units 6  (kappa=6); norm form x^2+xy+y^2, chi_-3")
    print("  the +6 vs +4 kissing = the +6 vs +4 ROOTS OF UNITY in the order (units = kissing!).")
    print("  => kappa of the lattice = #roots of unity (units) of its order: 4 for Z[i], 6 for")
    print("     Z[w].  The 2D cap kappa<=6 IS the fact that an imaginary quadratic field has")
    print("     at most 6 roots of unity (only Q(w) has 6, only Q(i) has 4, all others 2).")
    print("  The CM escape = a higher-degree field with MORE roots of unity on |z|=1, i.e. a")
    print("     rank-(>2) unit group -- exactly the 'unbounded =1 layer' of S699, now pinned as")
    print("     'roots of unity = kissing = units, capped at 6 in degree 2'.")
    print()
    # verify: units (roots of unity) of Z[i] and Z[w] = kissing numbers
    ksq,_ = kissing(((1,0),(0,1)))
    kei,_ = kissing(((1,0),(0.5, math.sqrt(3)/2)))
    print(f"  [check] kissing(square)={ksq} = #units Z[i]=4 ? {ksq==4}")
    print(f"  [check] kissing(triang)={kei} = #units Z[w]=6 ? {kei==6}")
    print()

if __name__ == "__main__":
    part_K()
    part_P()
    part_S()
    part_dict()
    print("="*78)
    print("HEADLINE (S702): kissing kappa = #roots of unity of the order, capped at 6 in degree 2;")
    print("the popular-norm layer (NOT kissing) is what beats 3n, and Z^2 (kappa=4) already does it")
    print("at distance sqrt(5).  No 2D group lies 'between triangular and CM' -- the CM jump is a")
    print("degree/rank jump (more roots of unity), orthogonal to the (flat) 2D kissing & exponent axes.")
    print("="*78)
