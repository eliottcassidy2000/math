"""
monad-explorer-2026-06-06-S702
==============================
THE KISSING CAP IS NOT THE UNIT-DISTANCE CAP: the norm-layer spectrum.

Builds on opus-S699 (trienerment '=1' layer = norm-1 layer; 2D lattices kissing-capped
kappa<=6 => "<= 3n"; CM field escapes via field-norm). opus's OPEN PROBE:
  "is there a 2D-realizable GROUP between the triangular lattice (kappa=6) and the CM
   field whose norm-1 layer beats 3n at moderate n?"

CLAIM OF THIS SCRIPT (the correction): opus conflates TWO different "norm-1 layers".
  (A) MULTIPLICATIVE norm-1 layer = roots of unity / units of finite order in the order
      = the MINIMAL (kissing) vectors. Genuinely capped at 6 in 2D (ZZ[i]:4, ZZ[w]:6,
      everything else: 2). This is opus's kappa.
  (B) ADDITIVE norm-R layer = lattice vectors of a FIXED, FREELY-CHOSEN length sqrt(R).
      This is what the unit-distance graph actually uses (x-y has length 1). r_Q(R) is
      UNBOUNDED over R for EVERY 2D lattice with a multiplicative norm form.

The Erdos sqrt(n) x sqrt(n) grid (n^{1+c/loglog n} unit distances) IS a subset of the
square lattice ZZ^2 and ALREADY beats 3n by choosing a NON-minimal radius. So no group
"between" triangular and CM is needed: ZZ^2 itself beats 3n at moderate n.

We verify rigorously:
  1. r_Q(R) unbounded for square (r2(5^k)=4(k+1)) and triangular (r(7^k)=6(k+1)).
  2. Smallest R with r_Q(R) > 6, 12, 24 for each class-number-1 form.
  3. EXACT finite-patch unit-distance count: smallest disk patch N where U(N) > 3N,
     at the best non-minimal radius. (the explicit "beats 3n at moderate n" threshold.)
  4. The full spectrum over all 9 class-number-1 imaginary quadratic forms:
     kissing # (additive minimal layer), root-of-unity # (multiplicative), max additive layer.
"""
import math
from collections import defaultdict

# ----------------------------------------------------------------------------
# Binary quadratic forms.  A 2D lattice with Gram matrix [[2a,b],[b,2c]] has
# norm form Q(x,y) = a x^2 + b x y + c y^2.  We use the principal reduced form
# of each class-number-1 imaginary quadratic field, disc = b^2-4ac < 0.
#   disc -3:  x^2+xy+y^2       (Eisenstein ZZ[w], triangular lattice) kappa=6
#   disc -4:  x^2+y^2          (Gaussian   ZZ[i],  square lattice)    kappa=4
#   disc -7:  x^2+xy+2y^2
#   disc -8:  x^2+2y^2
#   disc -11: x^2+xy+3y^2
#   disc -19: x^2+xy+5y^2
#   disc -43: x^2+xy+11y^2
#   disc -67: x^2+xy+17y^2
#   disc -163:x^2+xy+41y^2
FORMS = {
    -3:  (1,1,1),
    -4:  (1,0,1),
    -7:  (1,1,2),
    -8:  (1,0,2),
    -11: (1,1,3),
    -19: (1,1,5),
    -43: (1,1,11),
    -67: (1,1,17),
    -163:(1,1,41),
}
# number of roots of unity (units of finite order) in O_{disc}: 6 for -3, 4 for -4, else 2
WK = {-3:6, -4:4}

def rep_counts(a,b,c, Rmax):
    """r_Q(R) for 0<R<=Rmax, by brute force over a box."""
    cnt = defaultdict(int)
    # bound the box: Q(x,y) <= Rmax.  Use crude radius bound.
    B = int(math.isqrt(Rmax)) + 3
    Bx = B*4
    for x in range(-Bx, Bx+1):
        for y in range(-Bx, Bx+1):
            R = a*x*x + b*x*y + c*y*y
            if 0 < R <= Rmax:
                cnt[R]+=1
    return cnt

def kissing(a,b,c):
    """minimal nonzero value of Q and its multiplicity = additive minimal layer = kappa."""
    cnt = rep_counts(a,b,c, 50)
    Rmin = min(cnt)
    return Rmin, cnt[Rmin]

def smallest_R_with(cnt, thresh):
    cand = [R for R,m in cnt.items() if m >= thresh]
    return (min(cand), cnt[min(cand)]) if cand else (None,None)

print("="*78)
print("1.  r_Q(R) IS UNBOUNDED (additive norm-R layer) -- closed-form check")
print("="*78)
# square form r2(R) = 4 * sum_{d|R} chi_{-4}(d); r2(5^k)=4(k+1)
def r2(R):
    s=0; d=1
    while d*d<=R:
        if R%d==0:
            for dd in ({d, R//d}):
                if dd%4==1: s+= 1
                elif dd%4==3: s-= 1
        d+=1
    return 4*s
print("Square ZZ^2 (disc -4): r2(5^k) = 4(k+1) -> infinity")
for k in range(0,7):
    R=5**k
    print(f"   R=5^{k}={R:<8d}  r2={r2(R):3d}   (formula 4(k+1)={4*(k+1)})")
# triangular form r(R) = 6 * sum_{d|R} chi_{-3}(d); r(7^k)=6(k+1)
def r3(R):
    s=0; d=1
    while d*d<=R:
        if R%d==0:
            for dd in ({d, R//d}):
                m=dd%3
                if m==1: s+=1
                elif m==2: s-=1
        d+=1
    return 6*s
print("Triangular ZZ[w] (disc -3): r(7^k) = 6(k+1) -> infinity")
for k in range(0,5):
    R=7**k
    print(f"   R=7^{k}={R:<8d}  r ={r3(R):3d}   (formula 6(k+1)={6*(k+1)})")
print("=> BOTH 2D lattices have UNBOUNDED additive norm-R layers; kappa<=6 caps ONLY R=Rmin.")

print()
print("="*78)
print("2.  SPECTRUM over the 9 class-number-1 forms")
print("    minimal(kissing) layer vs roots-of-unity vs MAX additive layer (R<=20000)")
print("="*78)
print(f"{'disc':>5} {'form':>16} {'kappa(min)':>10} {'#roots1':>8} {'maxLayer':>9} {'@R':>7} "
      f"{'R>6':>10} {'R>12':>11} {'R>24':>11}")
spectrum={}
for disc,(a,b,c) in FORMS.items():
    cnt = rep_counts(a,b,c, 20000)
    Rmin = min(cnt); kap = cnt[Rmin]
    maxlayer = max(cnt.values()); argmax = min(R for R,m in cnt.items() if m==maxlayer)
    wk = WK.get(disc,2)
    r6  = smallest_R_with(cnt,8)   # strictly >6 => >=8 (layers are even multiples here)
    r12 = smallest_R_with(cnt,14)  # strictly >12
    r24 = smallest_R_with(cnt,26)  # strictly >24
    formstr=f"{a}x2+{b}xy+{c}y2".replace("+0xy","")
    spectrum[disc]=dict(Rmin=Rmin,kap=kap,maxlayer=maxlayer,argmax=argmax,wk=wk)
    print(f"{disc:>5} {formstr:>16} {kap:>10} {wk:>8} {maxlayer:>9} {argmax:>7} "
          f"{str(r6):>10} {str(r12):>11} {str(r24):>11}")
print("kappa(min) = additive MINIMAL layer (kissing).  #roots1 = MULTIPLICATIVE norm-1")
print("(roots of unity).  maxLayer = largest additive norm-R layer found (still growing).")

print()
print("="*78)
print("3.  EXACT finite-patch threshold: smallest disk patch N with U(N) > 3N")
print("    (the explicit 'beats 3n at MODERATE n', at a NON-minimal radius)")
print("="*78)

# EXACT integer-index computation.  A point is index (i,j); squared distance between
# (i1,j1),(i2,j2) is Q(i1-i2, j1-j2) = a du^2 + b du dv + c dv^2 (exact integer).
# A "disk patch of radius rho" = { (i,j) : Q(i,j) <= rho^2 } (Q itself is the squared norm).
def patch_indices(a,b,c, R2):
    pts=[]
    B=int(math.isqrt(R2))+ max(3, c)  # generous index bound
    for i in range(-B,B+1):
        for j in range(-B,B+1):
            if a*i*i+b*i*j+c*j*j <= R2:
                pts.append((i,j))
    return pts

def count_unit_distances_exact(a,b,c, pts, d2):
    """# unordered pairs (p,q) with Q(p-q)=d2, EXACT integer arithmetic, via offset shells."""
    S=set(pts)
    # precompute the offset shell: all (du,dv) with Q=d2
    B=int(math.isqrt(d2))+max(3,c)
    offs=[(du,dv) for du in range(-B,B+1) for dv in range(-B,B+1)
          if a*du*du+b*du*dv+c*dv*dv==d2]
    U=0
    for (i,j) in pts:
        for (du,dv) in offs:
            if (i+du,j+dv) in S: U+=1
    return U//2  # each pair counted twice

# square ZZ^2 (disc -4): best small non-minimal radius with layer 12 is R=25 (dist 5)
# triangular ZZ[w] (disc -3): layer 12 at R=7 (dist sqrt7) -- even smaller R
configs = [
    ("square ZZ^2 (disc -4)", (1,0,1),  25, 12),  # Q=x^2+y^2=25, layer 12
    ("triangular  (disc -3)", (1,1,1),   7, 12),  # Q=x^2+xy+y^2=7, layer 12
]
for name,(a,b,c),d2,layer in configs:
    print(f"\n  {name}: non-minimal radius sqrt({d2}), infinite-lattice degree = {layer} "
          f"(=> infinite density {layer//2}n >> 3n)")
    found=False
    for R2 in range(d2, 600):
        pts=patch_indices(a,b,c,R2)
        N=len(pts)
        if N<8: continue
        U=count_unit_distances_exact(a,b,c,pts,d2)
        if U>3*N and not found:
            print(f"     FIRST N with U>3N:  patch radius^2={R2:3d}  N={N:4d}  U={U:5d}  "
                  f"U/N={U/N:.3f}  (3N={3*N})")
            found=True
        if R2==500:
            print(f"     at radius^2={R2}: N={N:4d}  U={U:6d}  U/N={U/N:.3f}  (ratio -> {layer//2})")
    if not found:
        print("     (no crossing in range)")

print()
print("="*78)
print("4.  VERDICT on opus-S699's probe")
print("="*78)
print("""  No 2D group BETWEEN triangular and CM is required to beat 3n.  The square lattice
  ZZ^2 ITSELF beats 3n at moderate N, using a NON-minimal radius (additive norm-R layer).
  The kappa<=6 cap binds ONLY the minimal-vector / root-of-unity (MULTIPLICATIVE) layer.

  The genuine distinction of the CM field is the GROWTH RATE, not the 3n line:
    - 2D lattice, radius FIXED (best non-min): degree = const K>6  => Theta(K/2 * n)  [linear, > 3n]
    - 2D lattice, radius GROWS with n (Erdos grid): n^{1 + Theta(1/loglog n)}          [super-linear]
    - CM field (Sawin/OpenAI):                       n^{1.0+ }  (> n^{1.014})           [larger exponent]
  All three exceed 3n.  The spectrum is a spectrum of EXPONENTS (1, 1+o(1), 1.014+),
  NOT a threshold of 'crossing 3n'.  opus's kissing-cap statement is about the MINIMAL
  layer; the unit-distance graph lives on the (unbounded) ADDITIVE layer.""")
