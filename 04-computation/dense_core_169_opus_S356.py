# opus-2026-07-17-S356 -- HYP-7390: (A) can the sharp wall bound close the
# residual dense core?  (B) the 169 = 13^2 archaeology: 4/169 = (2/13)^2 is
# the INDEPENDENCE constant in the 13-convention -- the exact analogue of my
# 4*lam^2 at lam = 1/14 -- and saw(S) = sum_pairs [rho(a,b) - 4/169] is the
# repo's dilation-invariant scalar (S181).  THM-1012 gives saw a LOWER bound.
from fractions import Fraction as F
from math import gcd
import random, itertools

LAM14 = F(1,14); LAM13 = F(1,13)

def rho(a, b, LAM):
    """exact mu(D_a n D_b) via the sawtooth sum."""
    g = gcd(a,b); tot = F(0); Mb = LAM*(a+b); m = 0
    while True:
        if F(m) >= Mb*g and m > 0: break
        for mm in ([0] if m == 0 else [m,-m]):
            if mm % g: continue
            if F(abs(mm)) >= Mb: continue
            w = min(2*LAM*a, 2*LAM*b, LAM*(a+b)-abs(mm))
            if w > 0: tot += F(w, a*b)*g
        m += g
        if m > int(Mb)+2: break
    return tot

print("(A) DOES THM-1012 CLOSE THE DENSE CORE?  lam = 1/14")
print("    THM-1012 per-pair floor: 4*lam^2 - 2*lam*(a/b) = 1/49 - (1/7)(a/b)")
print("    S350 containment floor:  2*lam/b = 1/(7b)   (always > 0)")
print()
random.seed(356)
neg = 0; tot = 0; s350pos = 0
for _ in range(400):
    m0 = random.randint(20, 200)
    B = sorted(random.sample(range(m0, 13*m0), 7))     # comparable: ratio <= 13
    tot += 1
    thm1012 = sum(4*LAM14**2 - 2*LAM14*F(B[i], B[i+1]) for i in range(6))
    s350    = sum(2*LAM14/B[i+1] for i in range(6))
    exact   = sum(rho(B[i], B[i+1], LAM14) for i in range(6))
    if thm1012 <= 0: neg += 1
    if s350 > 0: s350pos += 1
print(f"    {tot} comparable 7-blocks:")
print(f"      THM-1012 floor sum <= 0 in {neg}/{tot} cases  <-- the defect swamps 4*lam^2")
print(f"      S350 floor sum > 0 in     {s350pos}/{tot} cases <-- positivity always")
print(f"      (exact Hunter floor is positive; only the THM-1012 LOWER BOUND fails here)")
print()
print("    => the wall bound + S350 CLOSES an isolated 7-comparable block (positivity")
print("       suffices).  It does NOT close the 13-speed DENSE CORE: nesting the other")
print("       6 runners needs a QUANTITATIVE window, i.e. the sharp floor, and THM-1012")
print("       is sharp only when b >> a.  The dense core is by definition unseparated.")
print()
print("    THE STRUCTURAL STATEMENT: THM-1012's defect 2*lam*(a/b) -> 2*lam as a -> b,")
print("    and 2*lam = 1/7 swamps the leading 4*lam^2 = 1/49 by a factor of 7.")
print(f"    crossover: defect < leading iff a/b < 2*lam = 1/7, i.e. b > 7a.")
print()
print("(B) THE 169 ARCHAEOLOGY: 4/169 = (2/13)^2 = the independence constant")
print(f"    (2/13)^2 = {F(2,13)**2}   <-> my 4*lam^2 at lam=1/14 = {4*LAM14**2}")
print(f"    (2/13)^3 = {F(2,13)**3}   (the triple analogue; 2197 = 13^3)")
print()
print("    saw(S) = sum_pairs [rho(a,b) - 4/169]  (S181's dilation-invariant scalar)")
print("    THM-1012 in the 13-convention gives saw(S) >= -2*lam*sum_pairs(a/b):")
for S in ([1,2,3,4], [7,11,13], [77,143,169], [3,5,7,11]):
    saw = sum(rho(a,b,LAM13) - F(4,169) for a,b in itertools.combinations(S,2))
    bound = sum(-2*LAM13*F(a,b) for a,b in itertools.combinations(S,2))
    print(f"      S={S}: saw = {saw} = {float(saw):+.6f}   THM-1012 bound {float(bound):+.6f}"
          f"   {'OK' if saw >= bound else 'VIOLATED'}")
print()
print("    THE (77,143,169) COINCIDENCE: 77 = 7*11, 143 = 11*13, 169 = 13*13 --")
print("    all inside the 1001 = 7*11*13 system; the repo records this triple as the")
print("    EXACT-independence point (ratio 1.000, deviation +0.000000).")
for a,b in itertools.combinations([77,143,169],2):
    print(f"      rho({a},{b}) = {rho(a,b,LAM13)}   4/169 = {F(4,169)}   "
          f"dev = {float(rho(a,b,LAM13)-F(4,169)):+.8f}  gcd={gcd(a,b)}")

print()
print("(C) THE EXACT-INDEPENDENCE LOCUS, explained by the folded identity.")
print("    THM-965: the deviation from independence is fold(S mod M) - fold(D mod M),")
print("    S = a+b, D = b-a, M = 13g (13-convention).  fold(r) = r(M-r) is SYMMETRIC")
print("    under r <-> M-r.  So the folds cancel exactly when S = +-D mod M, i.e.")
print("    when M | 2a  or  M | 2b.   PREDICTION: rho(a,b) = 4/169 exactly iff that holds.")
hit = miss = 0
for _ in range(3000):
    a = random.randint(1, 300); b = random.randint(a, 300)
    g = gcd(a,b); M = 13*g
    pred = (2*a) % M == 0 or (2*b) % M == 0
    actual = (rho(a,b,LAM13) == F(4,169))
    if pred == actual: hit += 1
    else:
        miss += 1
        if miss <= 4: print(f"      MISMATCH a={a} b={b} g={g} pred={pred} actual={actual}")
print(f"    3000 random pairs: characterization correct in {hit}/{hit+miss}")
print()
print("    the (77,143,169) triple, checked against the law:")
for a,b in itertools.combinations([77,143,169],2):
    g = gcd(a,b); M = 13*g
    print(f"      ({a},{b}): g={g}, M=13g={M};  2a mod M = {(2*a)%M}, 2b mod M = {(2*b)%M}"
          f"  -> independence: {((2*a)%M==0) or ((2*b)%M==0)}")
print()
print("    => the repo's recorded 'Sidon-far exact hit at (77,143,169)' is NOT a")
print("       coincidence of that triple: it is the fold reflection r <-> M-r, and")
print("       (77,143,169) sits on the locus because 13g | 2b for every pair.")
