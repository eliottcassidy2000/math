from fractions import Fraction as F

def nrm(x):
    r = x - int(x)
    r = r + 1 if r < 0 else r
    return r if r <= F(1,2) else 1 - r

def g(S, t):
    return min(nrm(v*t) for v in S)

def cand(S):
    S = sorted(set(S)); C = set()
    for v in S:
        k = 0
        while F(2*k+1, 2*v) <= F(1,2):
            C.add(F(2*k+1, 2*v)); k += 1
    for i in range(len(S)):
        for j in range(i+1, len(S)):
            for d in (S[i]+S[j], S[j]-S[i]):
                if d > 0:
                    k = 1
                    while F(k, d) <= F(1,2):
                        C.add(F(k, d)); k += 1
    C.add(F(1,2)); return C

def M(S):
    b = F(0); at = None
    for t in cand(S):
        v = g(S, t)
        if v > b:
            b = v; at = t
    return b, at

# ============================================================
# CLAIM (1): a runner ==7 mod 14 is ALWAYS safe at the grid antipode.
# At grid a/14, runner v sits in section (v*a mod 14). v==7 mod14:
#   v*a == 7a mod 14. For a odd: 7a == 7 mod14 (section 7, dist 1/2 = FREE).
#   For a even: 7a == 0 mod14 (section 0!). So at EVEN grid points a section-7
#   runner lands in section 0. But ||v * a/14|| : v==7 mod14, a even=2b:
#   v*a/14 = 7*2b/14 * (v/7) ... let's just compute ||v*(a/14)|| directly.
# The QUESTION: does a ==7 runner ever CONSTRAIN grid loneliness, i.e. is it
# ever the minimizer (closest to 0) at some grid point a/14?
print("="*70)
print("CLAIM (1): does a runner v==7 mod14 ever constrain grid loneliness?")
print("="*70)
# grid loneliness: tau = a/14, a=1..13. runner lonely if min_v ||v a/14|| achieved.
# A runner v "constrains" if ||v a/14|| is small (could be the min / kill loneliness).
for v in [7, 21, 35, 49]:
    print(f"\n  runner v={v} (v mod14={v%14}): ||v*a/14|| for a=1..13:")
    row = []
    for a in range(1, 14):
        d = nrm(F(v*a, 14))
        row.append((a, d))
    print("   ", "  ".join(f"a={a}:{str(d)}" for a,d in row))

print("\n" + "="*70)
print("CLAIM (1) refined: section-7 runner at EVEN grid points lands in section 0")
print("="*70)
# A section-7 runner v==7 mod14: at a even, ||v a/14|| = 0. That means it IS in
# section 0 at even grid points and WOULD make the observer NOT lonely there.
# So it does constrain — but only at even a. The point: 7-runners are 'safe for
# THEMSELVES' (antipode at odd a) but as OTHER runners they hit section 0 at even a.
# Key clarification needed: LRC loneliness is about the OBSERVER runner being lonely.
# Reframe: M(S) = max_tau min_v ||v tau||. A v==7 mod14 contributes ||7a/14||=0 at
# even a, =1/2 at odd a. So at EVEN grid points a 7-runner FORCES min=0 -> kills that
# grid point. At ODD grid points it contributes 1/2 (never the min unless alone).
# Conclusion: a 7-runner does NOT help the gap at even grid; it only is "free" at odd.
# Let's verify by example: does adding a 7-runner ever change M via the grid?
for base in [[1,2,3,4,5,6], [2,3,5,8,11], [1,3,9,11,13]]:
    mb, tb = M(base)
    S7 = base + [7]
    m7, t7 = M(S7)
    S21 = base + [21]
    m21, t21 = M(S21)
    print(f"  base={base}: M={mb} @ {tb}")
    print(f"    +7  : M={m7} @ {t7}   (changed: {m7!=mb})")
    print(f"    +21 : M={m21} @ {t21}  (changed: {m21!=mb})")

print("\n" + "="*70)
print("CLAIM (1) PRECISE: at GRID points a/14, is a 7-runner ever the min?")
print("="*70)
# At grid a/14: section-7 runner gives ||7a/14|| in {0 (a even), 1/2 (a odd)}.
# - a odd: contributes 1/2 = max possible distance -> NEVER the unique min (unless S={7r}).
# - a even: contributes 0 -> it IS a min, but 0 means observer NOT lonely there anyway
#           (some runner in section 0). So it doesn't *create* loneliness; if anything
#           it removes the even grid points as loneliness witnesses.
# So: a 7-runner never CONSTRAINS grid loneliness in the sense of being the
# binding small-distance runner at an ODD grid point. Verify exhaustively:
import itertools
viol = 0
examples = []
for combo in itertools.combinations([1,2,3,5,9,11,13], 3):  # mix of residues
    S = list(combo) + [7]
    for a in range(1,14):
        if a % 2 == 1:  # odd grid: 7-runner at 1/2
            d7 = nrm(F(7*a,14))
            assert d7 == F(1,2)
            others = [nrm(F(v*a,14)) for v in combo]
            # is 7 the strict min? only if all others > 1/2 (impossible)
            if d7 < min(others) if others else False:
                viol += 1
print(f"  Odd-grid violations (7-runner strict min): {viol}")
print("  => At odd grid points, 7-runner is always at 1/2, never the binding min.")
print("  => CLAIM (1) HOLDS: a 7-runner is FREE at odd grid (its own antipode);")
print("     at even grid it sits in section 0 (dist 0) = harmless, only removes")
print("     even grid points as loneliness witnesses. It never makes a config harder")
print("     at the grid; off-grid is a separate matter (it can lower M off-grid).")

print("\n" + "="*70)
print("CLAIM (2): can a transform move a section-0 (==0 mod14) runner to section-7?")
print("="*70)
# Section is determined by v mod 14. A 'parked' section-0 runner has v==0 mod14.
# A section-7 free runner has v==7 mod14.
# Test (Z/14)* units: multiply S by unit u. u*0==0 mod14 always. So units FIX section 0.
units = [u for u in range(1,14) if __import__('math').gcd(u,14)==1]
print(f"  (Z/14)* units: {units}")
print("  u*0 mod14 for each unit:", [ (u*0)%14 for u in units], "-> all 0. Units FIX section 0.")
print("  Reversal a->-a == mult by 13(=-1): -0==0, -7==7 mod14. FIXES both 0 and 7.")
print("  => Section 0 is a FIXED POINT of every (Z/14)* unit (0*u=0). So NO unit")
print("     can move a section-0 runner to section-7. CONFIRMED structurally.")
print()
# Multiply by 2: 0*2=0 (still parked!), but this is NOT a bijection on Z/14 (2 not a unit).
print("  Mult by 2: 0->0 (still section 0). 7->14==0 mod14! So *2 sends a 7-runner")
print("    INTO section 0 (free->parked, WRONG direction). And 2 is not invertible mod14.")
# Additive shift: v -> v + s changes the runner's SPEED, but speeds are fixed data.
# Adding a shift to tau (rotating observer) is the grid choice itself, not a relabel.
print("  Additive shift v->v+7: sends 0->7 (parked->free!) and 7->14==0 (free->parked).")
print("    BUT v+7 changes the actual runner SPEED -> a DIFFERENT runner, different config.")
print("    Does v->v+7 preserve M? Test:")
for base in [[14,1,2,3,4,5], [28,3,5,9,11], [14,2,4,8,16]]:
    mb,tb = M(base)
    shifted = [v+7 for v in base]
    ms,ts = M(shifted)
    print(f"    {base}: M={mb}  ->  shift+7 {shifted}: M={ms}  (preserved: {mb==ms})")

print("\n" + "="*70)
print("CLAIM (2) EXPLANATION: why 0 and 7 are BOTH fixed but distinguishable")
print("="*70)
# Z/14 = Z/2 x Z/7 (CRT). A residue r <-> (r mod2, r mod7).
#   section 0 = (0,0). section 7 = (1,0).
# Units (Z/14)* act by multiplication = (Z/2)* x (Z/7)* = {1} x {1..6}.
#   On the Z/2 part: only unit is 1, so the Z/2 coordinate is UNTOUCHED.
#   So a unit can NEVER flip the Z/2-parity coordinate.  0 has parity 0, 7 has parity 1.
# => 0 and 7 differ exactly in the Z/2 coordinate, which units fix. That's the obstruction.
for r in [0,7]:
    print(f"  r={r}: CRT (mod2,mod7)=({r%2},{r%7})")
print("  Units fix the mod-2 coordinate -> cannot map parity-0 (sec 0) to parity-1 (sec 7).")
print("  The ONLY maps mixing them are additive (shift by odd), which change speeds.")
print()
print("="*70)
print("CLAIM (3): full-config reversal. Two reversals to distinguish:")
print("="*70)
# Reversal A (the LRC complement / observer reversal): tau -> -tau. M(S) is invariant
#   under tau->-tau trivially since ||v(-t)||=||vt||. So M(S)=M(S) -- no info.
# Reversal B (residue reversal r -> 14-r, the section duality): apply to runner residues.
#   This is multiply by -1 = 13 (a UNIT). 0->0, 7->7 fixed. M invariant since 13 is a unit
#   and ||v*13*t||... actually relabel speeds v -> -v mod something is not the speed map.
# The genuine "complement of the WHOLE config" in residue space: S_rev = {14 - (v mod14)}.
print("Reversal B: S -> {(-v) mod 14 represented as 14-(v%14)} acting on RESIDUES.")
print("  This is mult by -1 (a unit). Section 0 (0)->0, section 7 (7)->7: BOTH FIXED.")
print("  Other sections swap: 1<->13, 2<->12, ... 6<->8. The 0/7 occupants are pinned.")
print()
# Does M(S) = M(reverse S) for actual speed reversal? Multiplying ALL speeds by a unit u
# mod 14 doesn't make sense as speeds (they're integers). The right statement: M is
# invariant under tau->-tau, giving the residue-reversal symmetry of the GRID picture.
print("  At grid a/14, section of v is v*a mod14. Replacing a -> -a == 14-a (also a grid")
print("  point) sends section s -> -s mod14. So the grid picture is symmetric under a->-a,")
print("  which fixes section 0 and section 7 and swaps s<->14-s. M unchanged (same grid set).")
# Verify: g(S, a/14) == g(S, (14-a)/14) for all a.
ok = True
for S in [[1,2,3,4,5,6],[14,1,3,9],[7,2,4,8]]:
    for a in range(1,14):
        if g(S,F(a,14)) != g(S,F(14-a,14)):
            ok=False
print(f"  Verify g(S,a/14)==g(S,(14-a)/14) for all tested a: {ok}")

print("\n" + "="*70)
print("CRUX: can we replace a sec-0 parked runner by a sec-7 free runner, M-preserving?")
print("="*70)
# A 'hard' config has a parked sec-0 runner w (==0 mod14). Replace w by a sec-7 runner
# w' (==7 mod14). Does M stay >= 1/14 / does the swap preserve looseness?
import itertools
print("Test: A u {w}, w==0 mod14, vs A u {w'}, w'==7 mod14. Compare M.")
cores = [[1,2,3,4,5,6,7,8,9,10,11,12], [2,3,5,8,11,13], [1,4,9,11], [3,5,6,10,12]]
for A in cores:
    mA,_ = M(A)
    print(f"\n  core A={A}: M(A)={mA} ({'EASY' if mA>=F(1,14) else 'HARD'})")
    for w0 in [14, 28, 84]:  # sec-0 parked
        S0 = A+[w0]; m0,t0 = M(S0)
        # match by replacing with a sec-7 runner of comparable size
        for w7 in [7, 21]:
            if w7 in A: continue
            S7 = A+[w7]; m7,t7 = M(S7)
            print(f"    +sec0 {w0}: M={m0} ({'loose' if m0>=F(1,14) else 'HARD<1/14'})  | "
                  f"+sec7 {w7}: M={m7} ({'loose' if m7>=F(1,14) else 'HARD<1/14'})")
        break  # one parked example per core for brevity

print("\n" + "="*70)
print("SURPRISE CHECK: is adding sec-0 MORE M-preserving than adding sec-7?")
print("="*70)
# Hypothesis from data: adding a sec-0 (parked) runner w==0 mod14 OFTEN preserves M(A)
# exactly (M(S)=M(A)) for GENERIC w, because the parked runner sits at section 0 only on
# the grid and on the grid the observer was already lonely with witness elsewhere. The
# sec-7 runner, being a REAL nonzero speed, can introduce new off-grid binding pairs.
# Quantify: over many cores A and many w, count how often M(Au{w}) == M(A).
import itertools, math, random
random.seed(1)
def gcd(a,b): return math.gcd(a,b)
sec0_pres = 0; sec0_tot = 0
sec7_pres = 0; sec7_tot = 0
sec0_drop = []; sec7_drop = []
cores = []
pool = list(range(1,28))
for _ in range(400):
    k = random.randint(3,7)
    A = sorted(random.sample(pool, k))
    cores.append(A)
for A in cores:
    mA,_ = M(A)
    for w in [14,28,42,56,84,112,140]:  # sec-0 parked, various
        S = A+[w]; mS,_ = M(S); sec0_tot+=1
        if mS==mA: sec0_pres+=1
        else: sec0_drop.append((mA-mS, A, w))
    for w in [7,21,35,49,63,77,91]:  # sec-7 free, various
        if w in A: continue
        S=A+[w]; mS,_=M(S); sec7_tot+=1
        if mS==mA: sec7_pres+=1
        else: sec7_drop.append((mA-mS, A, w))
print(f"  sec-0 parked runner: M preserved {sec0_pres}/{sec0_tot} = {sec0_pres/sec0_tot:.1%}")
print(f"  sec-7 free   runner: M preserved {sec7_pres}/{sec7_tot} = {sec7_pres/sec7_tot:.1%}")
# When it drops, how big?
if sec0_drop:
    sec0_drop.sort(reverse=True)
    print(f"  sec-0 worst drop: {sec0_drop[0][0]} (A={sec0_drop[0][1]}, w={sec0_drop[0][2]})")
if sec7_drop:
    sec7_drop.sort(reverse=True)
    print(f"  sec-7 worst drop: {sec7_drop[0][0]} (A={sec7_drop[0][1]}, w={sec7_drop[0][2]})")

print("\n" + "="*70)
print("WHY sec-7 still dips M: it has a real speed -> off-grid binding pairs")
print("="*70)
# The sec-7 runner v==7 mod14 is FREE only AT THE GRID (its antipode at odd a). But M is
# a max over ALL tau, attained at binding pairs k/(v_a +- v_b). The sec-7 runner forms
# binding pairs with every other runner -> can pull the optimum down OFF grid.
A=[10,17,24]; mA,tA=M(A); S=A+[7]; mS,tS=M(S)
print(f"  A={A}: M(A)={mA} @ tau={tA}")
print(f"  A+{{7}}: M={mS} @ tau={tS}  (off grid? {tS.denominator%14!=0})")
print(f"    The drop comes from a binding pair involving 7, e.g. k/(7 +- v).")
print()
print("="*70)
print("IS THERE *ANY* M-PRESERVING sec0<->sec7 SWAP? The honest answer.")
print("="*70)
# A genuine M-preserving transform must be a symmetry of the WHOLE multiset of binding
# distances. (Z/14)* units are NOT speed symmetries (multiplying speeds by u changes
# the binding denominators v_a +- v_b unless u=1). The only exact M-symmetries are:
#   (i) tau -> -tau  (always), (ii) global speed scaling S -> c*S (M invariant: ||cv*t/c||).
print("  Exact M-symmetries: (a) tau->-tau [trivial], (b) S->c*S global scale.")
# Check (b): does scaling preserve M?
for A in [[3,5,7,9],[2,3,11],[14,1,3]]:
    mA,_=M(A)
    for c in [2,3,5]:
        mC,_=M([c*v for v in A])
        assert mC==mA, (A,c)
print("  Verified: M(c*S)=M(S) for all tested c (global scaling is exact M-symmetry).")
print("  But scaling by c does NOT change residues mod14 in a 0<->7 way unless c is")
print("  chosen mod14 -- and c must be a UNIT to be a bijection on speeds' residues,")
print("  and units FIX the parity bit -> CANNOT swap 0 and 7. So:")
print()
print("  HONEST CONCLUSION: there is NO M-preserving transform that turns a sec-0")
print("  parked (hard) runner into a sec-7 free (easy) runner. The 0/7 duality is a")
print("  symmetry of the GRID/section picture only (it fixes 0 and 7 as antipodal")
print("  reversal fixed points), NOT of the full off-grid M-landscape.")

print("\n" + "="*70)
print("FINAL: the precise asymmetry of sec-0 vs sec-7 ON THE GRID")
print("="*70)
# sec-0 runner w==0 mod14: ||w*a/14|| = 0 for ALL a (always section 0, the center).
# sec-7 runner v==7 mod14: ||v*a/14|| = 1/2 (a odd), 0 (a even).
# So ON THE GRID:
#   - sec-0 runner kills EVERY grid point (dist 0 everywhere) -> if present, grid gives 0,
#     so grid loneliness is IMPOSSIBLE; M must be attained OFF grid. THIS is why parked
#     runners make configs 'hard' on the grid (they defeat every grid witness).
#   - sec-7 runner kills only EVEN grid points; ODD grid points still available (it's at 1/2).
for w,lab in [(14,'sec-0'),(7,'sec-7')]:
    row=[str(nrm(F(w*a,14))) for a in range(1,14)]
    print(f"  {lab} runner {w}: grid dists a=1..13: {row}")
print()
print("  => sec-0 is the GRID-HARDEST occupant (zero at every grid point);")
print("     sec-7 is GRID-SAFEST (1/2 at every odd grid point, its own antipode).")
print("  This IS the 0/7 duality's exact role: they are the two reversal fixed points,")
print("  and on the grid they are OPPOSITE extremes (all-0 vs alternating-1/2).")
print("  But OFF the grid both are real speeds forming binding pairs, so the duality")
print("  does NOT extend to an M-preserving hard->easy swap. The hard->easy mechanism")
print("  (THM-524 family) works via the EASY CORE + small resonance dip, NOT via 0<->7.")
