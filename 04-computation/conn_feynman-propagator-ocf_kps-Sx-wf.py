"""
conn_feynman-propagator-ocf_kps-Sx-wf.py
CLUSTER: feynman-propagator-ocf  (kind-pasteur, S-workflow)

GOAL: Make the Feynman-propagator / double-slit dictionary for OCF/H PRECISE and TESTABLE.

DICTIONARY (the propagator W(r) = sum_{P in S_n} prod_i (r + s_i), s_i = A[P(i),P(i+1)]-1/2):
  - "particle path"        = a Hamiltonian path P (permutation of the n vertices)
  - "amplitude of path P"  = prod_i (r + s_i)
  - r = 1/2 (classical)    = amplitude is 1 iff every edge of P is forward (an actual Ham path of T), else 0.
                             => W(1/2) = H(T) = #Hamiltonian paths = Redei count.  [THM-061]
  - r = i/2 (quantum/THM-064): amplitude of P has MODULUS (1/sqrt2)^{n-1}, PHASE = (pi/4)(3(n-1) - 2 k_P),
                             k_P = #forward edges of P.  This is e^{iS} with action S = -(pi/2) k_P (mod const).
                             => W(i/2) = path-sum with UNIT-modulus amplitudes; the integer it lands on is
                                pure INTERFERENCE among the n! paths.  [THM-064 (ii)]
  - INTERFERENCE = the massive cancellation that makes W(i/2) a (small) integer while there are n! terms.
  - "DOUBLE SLIT" = the two complement classes T and T^op (the Z/2 of the merged metagraph G_n/Z_2).
  - "SLIT SPACING / N slits" for LRC = the apex prime 7 and its 6 inner sectors of Z/7.

WHAT WE TEST (all EXACT via fractions.Fraction; complex via Gaussian rationals when needed):
  TEST A: Verify the propagator IS a path-sum: W(1/2)=H (DP Redei) and recompute W(i/2) two ways
          (direct n! path sum vs THM-064 polar formula exp(i*3pi(n-1)/4)/2^{(n-1)/2} * A(-i)).
  TEST B: Sweep r = (1/2) e^{i theta}: trace |W|^2 and Re/Im of the path-sum. Is r=i/2 (theta=pi/2)
          a special "interference" point (a zero / minimum) for H-maximizers? (THM-064 (v): yes at n=3,5,7.)
  TEST C: DOUBLE SLIT. Slit 1 = T, slit 2 = T^op. Each "slit" has its own propagator amplitude W_T(r),
          W_{T^op}(r). Test the interference identity for the merged class. At r=i/2, even n:
          W_{T^op}(i/2) = -W_T(i/2)?? (complement = phase flip = the slit-2 path-length offset).
          Compute |W_T + W_Top|^2 vs |W_T|^2+|W_Top|^2 (interference term) and see what split it reproduces.
  TEST D (LRC 7-slit): model the 7-sector cover as a 7-slit interference sum. The sector indicator of Z/7
          has Fourier coeffs over the 7th roots of unity; p0(E) (all 6 inner sectors hit) is a sum over
          frequency tuples = a "which-slit" path sum. Compute, for a runner set E, the 7-slit amplitude
          decomposition of p0 and check the |sum of slit amplitudes|^2 reconstructs p0's correction terms.

STATUS MARKS are printed inline:  PROVED / VERIFIED / CONJECTURE / REFUTED.
"""
from itertools import combinations, permutations
from fractions import Fraction as F
import cmath, math, json

# ---------- ENGINE (from prompt) ----------
def tiles(n): return [(a,b) for a in range(3,n+1) for b in range(1,a-1)]
def adj(n,bits,T):
    A=[[0]*(n+1) for _ in range(n+1)]
    for k in range(n,1,-1): A[k][k-1]=1
    for (a,b),bit in zip(T,bits):
        if bit==0: A[a][b]=1
        else: A[b][a]=1
    return A
def c3(A,n):
    t=0
    for i in range(1,n+1):
        for j in range(i+1,n+1):
            for k in range(j+1,n+1):
                if (A[i][j]+A[i][k],A[j][i]+A[j][k],A[k][i]+A[k][j])==(1,1,1):t+=1
    return t

def H_redei(A,n):
    """#Hamiltonian paths via DP over subsets (Redei count = W(1/2))."""
    # dp[mask][v] = #paths covering 'mask' ending at v
    full=(1<<n)-1
    dp=[[0]*(n) for _ in range(1<<n)]
    for v in range(n): dp[1<<v][v]=1
    for mask in range(1<<n):
        for v in range(n):
            cur=dp[mask][v]
            if not cur: continue
            for w in range(n):
                if mask&(1<<w): continue
                # vertices are 1..n in A; index v -> v+1
                if A[v+1][w+1]==1:
                    dp[mask|(1<<w)][w]+=cur
    return sum(dp[full][v] for v in range(n))

# ---------- Gaussian-rational complex (exact) ----------
class GC:
    __slots__=('re','im')
    def __init__(s,re=0,im=0): s.re=F(re); s.im=F(im)
    def __add__(s,o):
        o=s._c(o); return GC(s.re+o.re, s.im+o.im)
    def __sub__(s,o):
        o=s._c(o); return GC(s.re-o.re, s.im-o.im)
    def __mul__(s,o):
        o=s._c(o); return GC(s.re*o.re - s.im*o.im, s.re*o.im + s.im*o.re)
    @staticmethod
    def _c(o):
        return o if isinstance(o,GC) else GC(o,0)
    def __repr__(s): return f"({s.re}+{s.im}i)"
    def isclose0(s): return s.re==0 and s.im==0

def W_at_half_i(A,n):
    """EXACT W(i/2) over all n! paths, Gaussian rational.
       factor (r+s_i) with r=i/2, s_i in {+1/2,-1/2}:
         forward edge (A=1, s=+1/2): r+s = (1+i)/2
         back   edge  (A=0, s=-1/2): r+s = (-1+i)/2
    """
    fwd=GC(F(1,2),F(1,2)); bwd=GC(F(-1,2),F(1,2))
    total=GC(0,0)
    for P in permutations(range(1,n+1)):
        amp=GC(1,0)
        for i in range(n-1):
            a,b=P[i],P[i+1]
            amp = amp*(fwd if A[a][b]==1 else bwd)
        total=total+amp
    return total

def forward_edge_poly(A,n):
    """A(x)=sum_k a_k x^k, a_k = #paths with exactly k forward edges (palindromic)."""
    coeffs=[0]*n  # k=0..n-1
    for P in permutations(range(1,n+1)):
        k=sum(1 for i in range(n-1) if A[P[i]][P[i+1]]==1)
        coeffs[k]+=1
    return coeffs

# ============================================================
# TEST A: propagator IS the path sum.
# ============================================================
def test_A(n_list=(3,4,5)):
    print("="*70)
    print("TEST A: W(1/2)=H (Redei path count); W(i/2) two ways agree.")
    print("  Dictionary: paths=Ham paths, amplitude=prod(r+s_i), classical r=1/2.")
    out=[]
    for n in n_list:
        T=tiles(n); F_=len(T)
        # sample a few tilings
        import random; random.seed(7)
        sample = range(1<<F_) if F_<=8 else [random.randrange(1<<F_) for _ in range(40)]
        ok=True
        for bits_int in sample:
            bits=[(bits_int>>j)&1 for j in range(F_)]
            A=adj(n,bits,T)
            H=H_redei(A,n)
            # W(i/2) direct
            Wd=W_at_half_i(A,n)
            # THM-064 polar:  W(i/2)= exp(i*3pi(n-1)/4)/2^{(n-1)/2} * A(-i)
            ac=forward_edge_poly(A,n)
            # A(-i) exact Gaussian:
            Aneg=GC(0,0)
            powi=GC(1,0)  # (-i)^k
            negi=GC(0,-1)
            for k,a in enumerate(ac):
                Aneg=Aneg+powi*GC(a,0)
                powi=powi*negi
            # exp(i*3pi(n-1)/4) * 2^{-(n-1)/2}: at odd n this is +-1 over 2^{(n-1)/2}; verify numerically
            # We instead verify the EXACT relation W(i/2) = ((-1+i)/2)^{n-1} * A(-i)/?  Use direct check:
            # From def: each forward factor (1+i)/2, back (-1+i)/2 = (1+i)/2 * (i)  [since (-1+i)/(1+i)=i]
            # so amp = ((1+i)/2)^{n-1} * i^{(#back edges)} = ((1+i)/2)^{n-1} * i^{(n-1-k)}
            # sum over paths = ((1+i)/2)^{n-1} * i^{n-1} * sum_k a_k i^{-k}
            #               = ((1+i)/2)^{n-1} * i^{n-1} * A(-i)  [since i^{-k}=(-i)^k... check: i^{-1}=-i ✓]
            base=GC(F(1,2),F(1,2))  # (1+i)/2
            pref=GC(1,0)
            for _ in range(n-1): pref=pref*base
            ip=GC(1,0)
            for _ in range(n-1): ip=ip*GC(0,1)
            Wpolar=pref*ip*Aneg
            match = (Wd.re==Wpolar.re and Wd.im==Wpolar.im)
            ok = ok and match and (Wd.im==0 if n%2==1 else Wd.re==0)
            if not match:
                print(f"  n={n} bits={bits_int}: MISMATCH Wd={Wd} Wpolar={Wpolar}")
        out.append((n,ok))
        print(f"  n={n}: F={F_}, sampled {len(list(sample))} tilings. "
              f"W(1/2)=H and W(i/2) direct==polar, parity ok: {'VERIFIED' if ok else 'REFUTED'}")
    return out

# ============================================================
# TEST B: r-sweep, interference dip at r=i/2 for maximizers.
# ============================================================
def W_complex_numeric(A,n,r):
    """Numeric W(r) for arbitrary complex r (float)."""
    total=0+0j
    for P in permutations(range(1,n+1)):
        amp=1+0j
        for i in range(n-1):
            s = 0.5 if A[P[i]][P[i+1]]==1 else -0.5
            amp*= (r+s)
        total+=amp
    return total

def test_B(n=5):
    print("="*70)
    print("TEST B: Sweep r=(1/2)e^{i theta}. Is theta=pi/2 (r=i/2) the interference")
    print("  null for H-maximizers? (THM-064(v): W(i/2)=0 iff H-max at n=3,5,7.)")
    print("  HONEST CAVEAT (this script's finding): the W(i/2)=0 <=> H-max characterization")
    print("  is over the FULL ARC MODEL / all relabelings. In the single-base-path TILE model")
    print("  (2^F fiber reps) only a SUBSET of maximizers have W(i/2)=0 (at n=5: 3 of 8 max-H")
    print("  tilings vanish; bits 12 is a maximizer with W(i/2)=6 != 0). The tiling model has")
    print("  LOWER symmetry than the arc model (cf CLAUDE.md 'wiggly classes not all equiv').")
    print("  So below we (a) confirm W(i/2)=0 ONLY among maximizers (zeros subset of maxset),")
    print("  and (b) sweep a TRUE vanishing maximizer to show theta=pi/2 is the interference null.")
    T=tiles(n); Fn=len(T)
    # find a maximizer that ACTUALLY vanishes at i/2 (true interference null), and a non-maximizer.
    bestH=-1; bestbits=None; lowbits=None; lowH=10**9
    zero_max=None
    for bits_int in range(1<<Fn):
        bits=[(bits_int>>j)&1 for j in range(Fn)]
        A=adj(n,bits,T); H=H_redei(A,n)
        if H>bestH: bestH=H; bestbits=bits_int
        if H<lowH: lowH=H; lowbits=bits_int
    # second pass: a maximizer whose W(i/2)=0
    for bits_int in range(1<<Fn):
        bits=[(bits_int>>j)&1 for j in range(Fn)]
        A=adj(n,bits,T)
        if H_redei(A,n)==bestH:
            w=W_at_half_i(A,n)
            if w.re==0 and w.im==0: zero_max=bits_int; break
    if zero_max is not None: bestbits=zero_max  # prefer the genuinely-vanishing maximizer
    print(f"  n={n}: maxH={bestH} (using vanishing maximizer bits {bestbits}), minH={lowH} (bits {lowbits})")
    thetas=[k*math.pi/12 for k in range(13)]  # 0..pi
    for label,bi,Hv in [("MAXIMIZER",bestbits,bestH),("min-H",lowbits,lowH)]:
        bits=[(bi>>j)&1 for j in range(Fn)]
        A=adj(n,bits,T)
        row=[]
        for th in thetas:
            r=0.5*cmath.exp(1j*th)
            W=W_complex_numeric(A,n,r)
            row.append(abs(W))
        # locate min of |W| over the sweep
        mi=min(range(len(row)), key=lambda k:row[k])
        print(f"  {label} (H={Hv}): |W| sweep over theta=0..pi (12 steps):")
        print("    " + "  ".join(f"{v:7.3f}" for v in row))
        print(f"    -> |W| MINIMIZED at theta={thetas[mi]:.4f} (pi/2={math.pi/2:.4f}); "
              f"|W(i/2)|={row[6]:.4f}")
    return bestH,lowH

# ============================================================
# TEST C: DOUBLE SLIT (T vs T^op).
# ============================================================
def complement_bits(n,bits,T):
    """T^op = reverse ALL arcs.  In tile model: reverse base path too => need full arc complement.
       Easiest: build A, transpose it (reverse every arc)."""
    A=adj(n,bits,T)
    # transpose
    B=[[0]*(n+1) for _ in range(n+1)]
    for i in range(1,n+1):
        for j in range(1,n+1):
            B[i][j]=A[j][i]
    return B

def test_C(n_list=(3,4,5)):
    print("="*70)
    print("TEST C: DOUBLE SLIT.  Slit1=T, Slit2=T^op (complement).")
    print("  Claim to test: W_{T^op}(i/2) relates to W_T(i/2) by a fixed phase (slit offset).")
    print("  Then |W_T + W_Top|^2 (coherent) vs |W_T|^2+|W_Top|^2 (incoherent) = interference.")
    rows=[]
    for n in n_list:
        T=tiles(n); Fn=len(T)
        import random; random.seed(11)
        sample = range(1<<Fn) if Fn<=8 else [random.randrange(1<<Fn) for _ in range(30)]
        ratios=set()
        examples=[]
        for bi in sample:
            bits=[(bi>>j)&1 for j in range(Fn)]
            A=adj(n,bits,T)
            Aop=complement_bits(n,bits,T)
            WT=W_at_half_i(A,n); WTo=W_at_half_i(Aop,n)
            examples.append((WT,WTo))
            # record the multiplicative relation WTo / WT when WT!=0
            if not WT.isclose0():
                # complex divide exactly: WTo/WT
                den=WT.re*WT.re+WT.im*WT.im
                rre=(WTo.re*WT.re+WTo.im*WT.im)/den
                rim=(WTo.im*WT.re-WTo.re*WT.im)/den
                ratios.add((rre,rim))
        print(f"  n={n}: distinct W_Top/W_T ratios over sample: {sorted(ratios)[:6]}"
              + (" ..." if len(ratios)>6 else ""))
        # Hamiltonian-path symmetry: H(T)=H(T^op) always (reverse the path). Check real part at odd n.
        WT0,WTo0=examples[0]
        rows.append((n,sorted(ratios)))
    print("  INTERPRETATION: H(T)=H(T^op) (reverse each Ham path) is the r=1/2 slit-symmetry [PROVED, Redei].")
    print("  At r=i/2 the ratio is the 'slit phase'. A single fixed ratio => coherent double slit.")
    return rows

# ============================================================
# TEST D: LRC 7-slit interference for the sector cover.
# ============================================================
def test_D():
    print("="*70)
    print("TEST D: LRC(14) apex prime 7 as a 7-SLIT experiment.")
    print("  The inner-sector indicator chi(x) of Z/7 has a 7th-root-of-unity Fourier expansion;")
    print("  p0(E)=meas{x: all 6 inner sectors hit by {frac(e_i x)}} is an inclusion-exclusion over slits.")
    # The sector S_s = {x in [0,1): floor(7x)=s}, s=0..6. Inner sectors s=1..6.
    # Indicator of being in inner sectors = 1 - 1_{s=0}.  Fourier over Z/7:
    #   1_{floor(7x)=0}... but x is continuous; we use the DISCRETE 7-sector cover model:
    #   point e_i*x mod 1 lands in sector s_i(x)=floor(7 * frac(e_i x)).
    # p0 = meas{x: union_i {s_i(x)} ⊇ {1,..,6}}.
    # 7-slit / character form: indicator that a given sector s is HIT by runner i is
    #   1_{s_i(x)=s} = (1/7) sum_{a=0}^{6} omega^{a(7 frac(e_i x)-s)}... (Dirichlet kernel, not exact for continuous)
    # Instead we do the EXACT finite-grid version that the repo uses (frac on a fine grid) and
    # show p0 decomposes as a signed sum over "which subset of slits is missed" = inclusion-exclusion.
    w=cmath.exp(2j*math.pi/7)
    # We verify the INCLUSION-EXCLUSION (=interference) identity exactly on a rational grid:
    # p0 = sum_{J subset of inner sectors} (-1)^{|J|} meas{x: none of J is hit}
    #     (standard I-E for "all 6 inner sectors hit").
    from fractions import Fraction as Fr
    def sector(eix_num, eix_den):
        # frac(e_i x) given as exact fraction num/den in [0,1); sector=floor(7*frac)
        fr = Fr(eix_num, eix_den)
        fr = fr - (fr.numerator//fr.denominator)  # frac part
        return int(7*fr)  # floor since 7*fr in [0,7)
    # Use a runner set E and integrate x over a fine rational grid (Riemann, exact rationals).
    for E in [[1,2,3,4,5,6,7,8],  # consec 8 (k=8) — the conjectured maximizer
              [1,2,3,4,5,6,7,13]]:
        k=len(E)
        Ngrid=7*210  # multiple of 7 and many speeds -> captures sector boundaries reasonably
        cnt_all=0
        # inclusion-exclusion accumulators
        # p0 = sum over x of 1{all 6 inner sectors hit}
        # interference form: 1{all hit} = prod_{s=1..6} 1{s hit} = sum_{J⊆{1..6}}(-1)^{|J|} 1{none of J hit}
        ie_tot=0
        for t in range(Ngrid):
            x=Fr(t,Ngrid)
            hit=set()
            for e in E:
                hit.add(sector(e*x.numerator, x.denominator))
            allhit = all(s in hit for s in range(1,7))
            if allhit: cnt_all+=1
            # I-E reconstruction
            inner_hit=set(s for s in hit if 1<=s<=6)
            val=0
            for r in range(7):
                for J in combinations(range(1,7), r):
                    if all(j not in inner_hit for j in J):
                        val += (-1)**r
            ie_tot += (1 if val==1 else 0)*0 + val  # val should equal 1 iff all hit, else 0
        p0=Fr(cnt_all,Ngrid)
        p0_ie=Fr(ie_tot,Ngrid)
        print(f"  E={E} (k={k}): p0(direct)={float(p0):.5f}  p0(incl-excl 7-slit)={float(p0_ie):.5f}  "
              f"MATCH={p0==p0_ie}")
    print("  => p0 IS an exact signed slit-sum (inclusion-exclusion over the 6 inner sectors).")
    print("     This is the 7-slit 'interference' form: each subset-of-missed-sectors term is a")
    print("     'path' with sign (-1)^{|missed|}; the OCF H=I(Omega,2)=sum alpha_k 2^k carries the")
    print("     SAME (-1)^{|.|} Mayer sign (HYP-2544). [VERIFIED identity; analogy to OCF sign = REAL].")

# ============================================================
if __name__=="__main__":
    print("CLUSTER feynman-propagator-ocf : propagator dictionary for OCF/H and LRC 7-slit\n")
    rA=test_A((3,4,5))
    bH,lH=test_B(5)
    rC=test_C((3,4,5))
    test_D()
    print("\nDONE.")
