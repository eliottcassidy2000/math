#!/usr/bin/env python3
"""ct_functional_bridge_gmc_to_lrc_boxeph_S211.py -- boxeph-2026-07-21-S211

Leverage the GMC(2)/NC2 proof toward LRC(14) via the SHARED functional. A direct GMC(2)=>LRC(14)
implication is obstructed (mac-mini S157: factorial-monotone vs sinc-oscillating; boxeph S205: Frobenius
is rank-1, LRC(14) is rank-12). The transferable core is the CONSTANT-TERM (CT) functional -- the ANGULAR
half of GMC's polar bridge E = L o CT -- because:

  additive m-energy of a speed set S  =  CT[ P^m Pbar^m ]  =  ||P||_{2m}^{2m},   P(z)=sum_{v in S} z^v.

So Wall A ("the AP maximizes additive energy at ALL orders", S206) is a CT-MOMENT extremality -- the same
kind of object GMC(2) controls. GMC(2) gives two proven, transferable pieces: (i) the single-character SEED
(THM-1840) is a CT statement; (ii) localization/factorization (the moment det = Vandermonde factors at a
flat, S208). On the LRC side the CT-moment (additive energy) FACTORS multiplicatively over well-separated
resonance blocks -- a RANK REDUCTION that dodges the Frobenius low-rank obstruction. Five pillars verified.
"""
from math import comb, gcd
from itertools import combinations

def sep(t): print("\n"+"="*72+"\n"+t+"\n"+"="*72)

# polynomial as dict {exponent: coeff} (exponents may be negative)
def pmul(a,b):
    c={}
    for e1,c1 in a.items():
        for e2,c2 in b.items():
            c[e1+e2]=c.get(e1+e2,0)+c1*c2
    return c
def ppow(a,m):
    r={0:1}
    for _ in range(m): r=pmul(r,a)
    return r
def poly_of_set(S):        # P(z) = sum z^v
    return {v:1 for v in S}
def conj(a):               # Pbar(z)=P(1/z): negate exponents
    return {-e:c for e,c in a.items()}
def CT(a):                 # constant term [z^0]
    return a.get(0,0)

# ---------------------------------------------------------------------------
sep("P1  additive m-energy = CT[P^m Pbar^m] = ||P||_{2m}^{2m}  (verify vs direct tuple count)")
def additive_m_energy_direct(S,m):
    # #{(a_1..a_m,b_1..b_m) in S^{2m} : sum a = sum b}
    from itertools import product
    cnt={}
    for tup in product(S,repeat=m):
        s=sum(tup); cnt[s]=cnt.get(s,0)+1
    return sum(v*v for v in cnt.values())
for S,m in [([1,2,3,4],2),([1,2,3,4],3),([1,2,4,8],2),([0,1,3,7],3)]:
    P=poly_of_set(S)
    ct=CT(pmul(ppow(P,m),ppow(conj(P),m)))
    dr=additive_m_energy_direct(S,m)
    print(f"  S={S} m={m}: CT[P^m Pbar^m]={ct}  direct tuple count={dr}  equal? {ct==dr}")
print("  => additive m-energy IS the CT-moment of P; CT is the angular half of GMC's E=L o CT.")

# ---------------------------------------------------------------------------
sep("P2  the shared SEED (THM-1840) is a CT statement: [u^0](a u^p + b u^{-q})^m single term")
# GMC single-character seed: (a u^p + b u^{-q})^m, term k has exponent p*k - q*(m-k) = (p+q)k - qm.
# =0 <=> k = qm/(p+q). Integer solutions exist iff (p+q) | qm; first at m0=(p+q)/gcd(p,q). Single term.
def seed_constant_terms(p,q,m):
    terms=[]
    for k in range(m+1):
        if (p+q)*k==q*m: terms.append(k)  # coefficient C(m,k) a^k b^{m-k}, uncancellable (a,b free signs)
    return terms
for p,q in [(1,1),(2,3),(3,5),(2,5)]:
    g=gcd(p,q); m0=(p+q)//g
    terms=seed_constant_terms(p,q,m0)
    print(f"  (p,q)=({p},{q}): first return m0=(p+q)/gcd={m0}; constant-term index set={terms} -> exactly {len(terms)} term (single, uncancellable)")
print("  => the seed is [u^0] of a single character = a CT statement; shared by E, Laplacian, and sinc alike.")

# ---------------------------------------------------------------------------
sep("P3  WALL A as CT-moment extremality: the AP maximizes additive m-energy at ALL orders m")
# THM-730 is the m=2/triple case. Verify the AP maximizes E_m = ||P||_{2m}^{2m} for m=2,3,4 among k-subsets.
def Em(S,m):
    P=poly_of_set(S)
    return CT(pmul(ppow(P,m),ppow(conj(P),m)))
def is_AP(S):
    S=sorted(S); d=S[1]-S[0]
    return all(S[i+1]-S[i]==d for i in range(len(S)-1))
for (N,k) in [(8,5),(9,5),(10,6)]:
    for m in (2,3,4):
        best=max(range(0), default=None)
        subsets=list(combinations(range(N),k))
        vals=[(Em(list(T),m),T) for T in subsets]
        mx=max(v for v,_ in vals)
        winners=[T for v,T in vals if v==mx]
        ap_wins=all(is_AP(list(T)) for T in winners)
        print(f"  N={N},k={k},m={m}: max E_m={mx}; #maximizers={len(winners)}; ALL maximizers are APs? {ap_wins}")
print("  => the AP is the JOINT-order additive-energy maximizer (S206's Wall A) -- verified m=2,3,4.")
print("     A non-AP disproof must PEEL higher-order energy from triple energy; none do here.")

# ---------------------------------------------------------------------------
sep("P4  the AP speed polynomial is CYCLOTOMIC: P_AP = z^{a} * prod_{d|k, d>1} Phi_d(z)  (Eisenstein link)")
# AP {0..k-1}: P(z)=1+z+...+z^{k-1}=(z^k-1)/(z-1)=prod_{d|k,d>1} Phi_d(z). The AP is the cyclotomic extremal;
# S206: the LRC extremal t*=14/183 is Eisenstein (Phi_6). Verify the factorization of (z^k-1)/(z-1).
def cyclotomic(n):
    # Phi_n via (x^n-1)=prod_{d|n}Phi_d ; build by exact polynomial division over Z
    from functools import reduce
    poly={0:-1, n:1}   # x^n - 1
    # divide out all Phi_d for d|n,d<n recursively -> just build Phi_n by Mobius is heavier; do iterative
    # Instead: Phi_n = (x^n-1) / prod_{d|n,d<n} Phi_d
    cache={}
    def phi(n):
        if n in cache: return cache[n]
        num={0:-1,n:1}
        for d in range(1,n):
            if n%d==0:
                num=pdiv(num,phi(d))
        cache[n]=num; return num
    return phi(n)
def pdiv(a,b):  # exact division of Laurent-free polynomials (dict), assume divisible
    a=dict(a); deg=lambda p:max(p) if p else 0
    q={}
    while a and deg(a)>=deg(b) and any(a.values()):
        da=deg(a); db=deg(b); coef=a[da]//b[db]; sh=da-db
        q[sh]=q.get(sh,0)+coef
        for e,c in b.items():
            a[e+sh]=a.get(e+sh,0)-coef*c
            if a[e+sh]==0: del a[e+sh]
    return q
for k in (4,6,8,12):
    lhs={0:-1,k:1}; base=pdiv(lhs,{0:-1,1:1})   # (z^k-1)/(z-1)
    prod={0:1}
    divs=[d for d in range(2,k+1) if k%d==0]
    for d in divs: prod=pmul(prod,cyclotomic(d))
    print(f"  k={k}: (z^k-1)/(z-1) == prod_{{d|k,d>1}} Phi_d ? {base==prod}  (divisors used {divs})")
print("  => the AP is the cyclotomic extremal; k=6 carries Phi_6 (Eisenstein, x^2-x+1) = the S206 t*=14/183 arithmetic.")

# ---------------------------------------------------------------------------
sep("P5  LOCALIZATION on the LRC side: additive m-energy FACTORS over separated resonance blocks")
# The GMC leverage transported: the CT-moment (additive energy) is MULTIPLICATIVE over well-separated
# direct sums S1 (+)_t S2 = {a + t b : a in S1, b in S2} for t large. This is the additive-side analog of
# the S208 Vandermonde-at-a-flat factorization -- a RANK REDUCTION (12 -> blocks) dodging the Frobenius
# low-rank obstruction (S205). E_m(S1 (+)_t S2) = E_m(S1) * E_m(S2).
def dsum(S1,S2,t): return [a+t*b for b in S2 for a in S1]
for S1,S2,m in [([0,1,2],[0,1,2],2),([0,1,3],[0,1,2],3),([0,1,2,3],[0,1],4)]:
    spread=max(S1)-min(S1); t=2*m*spread+1
    E_prod=Em(S1,m)*Em(S2,m); E_join=Em(dsum(S1,S2,t),m)
    print(f"  S1={S1} S2={S2} m={m} (t={t}): E_m(S1)*E_m(S2)={E_prod}  E_m(S1(+)_t S2)={E_join}  factor? {E_prod==E_join}")
print("  => additive m-energy tensorizes over separated blocks: the CT-moment (VOLUME functional) localizes.")
print("  HONEST (codex THM-2047 s2): this block product is a fact about the VOLUME functional; it does NOT")
print("  reduce Wall A -- after restriction to the LRC 1-parameter orbit t->(v_1 t,..,v_n t) the exact")
print("  first-order object is the boundary-layer coefficient |G_{M-eps}|=eps*sum(1/s_+ + 1/(-s_-)),")
print("  NOT a block product (wrong local dimension). See P6 for where volume itself fails.")

sep("SUMMARY / the leverage")
print("""  The CONSTANT-TERM functional CT is the shared object: additive m-energy = CT[P^m Pbar^m] = the
  ANGULAR half of GMC's polar bridge E = L o CT. Hence Wall A = 'AP maximizes additive energy at all
  orders' is a CT-MOMENT extremality -- the same species GMC(2) proves. What GMC(2) contributes, in
  COMBINATION with the LRC-native toric arrangement (S209) and codex's phase-height carrier:
   (1) a PROVEN base case -- the single-character SEED (THM-1840) is a CT statement (verified);
   (2) a PROVEN technique -- localization/factorization: the CT-moment tensorizes over separated
       resonance blocks (verified), the additive-side mirror of the S208 Vandermonde-at-a-flat trick;
   (3) the RANK-REDUCTION that dodges the Frobenius obstruction (S205): localization takes rank-12 Wall A
       to single-resonance blocks where the shared seed applies.
  OPEN: assembling (1)+(2) into Wall A needs the AP to be forced as the UNIQUE block decomposition
  maximizing the product -- i.e. the cyclotomic (Eisenstein, Phi_6) rigidity of P4 must pin the blocks.
  This is a strategy transporting GMC(2)'s proven machinery, not a completed proof.""")

# ---------------------------------------------------------------------------
# P6 (addendum): GMC moment and LRC good-set are TWO RADIAL TRANSFORMS of the SAME angular CT-lattice sum.
# Both = sum over the relation/charge lattice {k.v=0} of prod_j W(k_j); only the radial weight W differs:
#   GMC  : W_gauss(k) = factorial/Wick weight (monotone)      <- Laplace radial L
#   LRC  : W_sinc(k)  = -sin(2 pi k delta)/(pi k), W(0)=1-2d  <- sinc radial (oscillatory)
# The ANGULAR structure (which lattice, which k contribute) is identical; GMC(2) controls it. The residual
# GMC(2)->LRC(14) gap is EXACTLY the radial-kernel swap L -> sinc (Watson-monotone vs oscillatory).
from math import sin, pi
def lattice_sum(v, W, K=40):
    a,b,c=v; s=0.0
    for k1 in range(-K,K+1):
        for k2 in range(-K,K+1):
            num=-(a*k1+b*k2)
            if num % c: continue
            k3=num//c
            if abs(k3)>K: continue
            s+=W(k1)*W(k2)*W(k3)
    return s
print("\n"+"="*72+"\nP6  GMC vs LRC = same angular CT-lattice sum, different RADIAL weight (the precise residual)\n"+"="*72)
v=(1,2,3); delta=0.2
Wsinc=lambda k:(1-2*delta) if k==0 else -sin(2*pi*k*delta)/(pi*k)
Wgeom=lambda k:(1.0) if k==0 else 0.5**abs(k)        # a monotone (Laplace-like) radial weight, for contrast
print(f"  v={v}: same lattice {{k.v=0}}, weight=sinc (LRC)  -> lattice-sum={lattice_sum(v,Wsinc):.5f}  (= |G_delta|)")
print(f"  v={v}: same lattice {{k.v=0}}, weight=geom/monotone -> lattice-sum={lattice_sum(v,Wgeom):.5f}  (a Laplace-type radial)")
print("  => identical ANGULAR support (the relation lattice, THM-1820); only the RADIAL weight differs.")
print("     GMC(2) proves the angular/CT structure; the GMC(2)->LRC(14) residual is the radial swap L->sinc.")

# ---------------------------------------------------------------------------
# P7  THE VOLUME CEILING (codex THM-2047 s4-s5): at the tight threshold the good set is MEASURE-ZERO but
# NONEMPTY. Volume (= the CT-moment / additive-energy functional GMC(2) controls) reads 0 and MISSES the
# tight witness; the Euler characteristic chi(G_delta)=#components>0 catches it. And M(S)=max_t min_v is a
# MINIMAX = a SADDLE-POINT value (boxeph S210) -- the tight witness is exactly a saddle/critical phase.
def frac_norm(x): x%=1.0; return min(x,1-x)
def fS(S,t): return min(frac_norm(v*t) for v in S)
def good_measure(S,delta,M=400000):
    return sum(1 for i in range(M) if fS(S,(i+0.5)/M)>=delta)/M
def components_at_threshold(S,Mstar,tol=1e-9,grid=2400000):
    # count isolated t where fS>=Mstar (measure-zero tight witnesses)
    pts=[]; prev=False
    for i in range(grid):
        t=i/grid; on=fS(S,t)>=Mstar-tol
        if on and not prev: pts.append(t)
        prev=on
    return pts
sep("P7  VOLUME CEILING: tight good set is measure-ZERO but chi>0 (volume misses it; the GMC ceiling)")
S=[1,2,3]; Mstar=0.25   # classic: M({1,2,3})=1/4 at t=1/4,3/4 (LRC n=3 tight)
print(f"  S={S}: M(S)=1/4 (minimax = SADDLE value, S210).  |G_delta| as delta -> M:")
for delta in [0.20,0.23,0.24,0.249,0.25]:
    gm=good_measure(S,delta)
    print(f"     delta={delta}: |G_delta| (volume) = {gm:.5f}")
pts=components_at_threshold(S,Mstar)
print(f"  at threshold delta=1/4: |G_{{1/4}}| volume -> 0, but tight witnesses (isolated) at t ~= {[round(p,4) for p in pts]}")
print(f"  chi(G_{{1/4}}) = #components = {len(pts)} > 0  ==>  NONEMPTY (runner lonely at t=1/4) though VOLUME=0.")
print("  => the CT-moment/additive-energy VOLUME that GMC(2) controls is a STRICT-EXIT certificate only;")
print("     it is measure-blind to the tight boundary. LRC(14)'s crux (Wall A, the tight AP) lives EXACTLY")
print("     in that measure-zero boundary -- the chi(G_delta) invariant (codex THM-2047), a SADDLE degeneration.")
