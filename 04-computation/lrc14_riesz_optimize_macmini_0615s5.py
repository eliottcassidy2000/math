#!/usr/bin/env python3
"""
lrc14_riesz_optimize_macmini_0615s5.py  (mac-mini-2026-06-15-S5; OPEN-Q-104 + OPEN-Q-097)

ATTACK THE LRC(14) ENDGAME inf_S L(S) > 0 ON TWO ROUTES.

ROUTE 104 (Bedert Riesz-product certificate). S is LOOSE if a probability density
R(τ)=∏_{j∈D}(1 + a_j cos 2π m_j τ) ≥ 0 has  ∫ M·R / ∫R < 1,  M(τ)=Σ_v 1_{||vτ||≤1/14}.
In Fourier:  ∫M·R = Σ_v Σ_k s(k) R̂(vk),  ∫R = R̂(0),  s(k)=sin(πk/7)/(πk).
R̂(n) = Σ_{ε∈{-1,0,1}^D, Σε_j m_j = n} ∏_{j:ε_j≠0}(a_j/2)·sign.  We OPTIMIZE D, a_j to push
the ratio below 1 for the extremizer cores.

ROUTE 097 (cross-level decomposition). L = (6/7)^13 + Σ_{|T|≥1}(-1)^|T|(6/7)^{13-|T|}Σ_T,
Σ_T = Σ_{exact relations supp=T}∏s. Compute the per-LEVEL mass Λ_k=Σ_{|T|=k}(6/7)^{13-k}Σ_T
to test whether |Σ_{|T|≥3}| < (6/7)^13 (the OPEN-Q-097 literal target) or whether the
cross-level alternation is essential.
"""
import sys, itertools, math
from math import gcd, sin, pi
from collections import defaultdict

sys.stdout.reconfigure(line_buffering=True)

def s(k):
    return sin(pi*k/7)/(pi*k) if k != 0 else 1/7

def lonely_measure(S, Q):
    rad = Q//14; c=0
    for a in range(Q):
        if all(not ((v*a)%Q<=rad or (v*a)%Q>=Q-rad) for v in S): c+=1
    return c/Q

# ---- ROUTE 097: level decomposition ----
def level_masses(S, B, Kmax_support=4):
    """Λ_k = (6/7)^{13-k} Σ_{|T|=k, exact 7-primitive relations, |coeff|<=B} ∏ s(t).
       Returns dict k -> signed level mass (the (-1)^k already folded in via term sign in L)."""
    n=len(S); base=6/7
    raw=defaultdict(float)  # k -> Σ_{|T|=k} (6/7)^{13-k} ∏ s(t)  (unsigned product of s; sign (-1)^k applied to total)
    seen=set()
    for sz in range(1, Kmax_support+1):
        for T in itertools.combinations(range(n), sz):
            for co in itertools.product([c for c in range(-B,B+1) if c!=0], repeat=sz):
                if sum(c*S[i] for c,i in zip(co,T))==0:
                    key=tuple(sorted(zip(T,co)))
                    if key in seen: continue
                    seen.add(key)
                    prod=1.0
                    for c in co: prod*=s(c)
                    raw[sz]+= (base**(n-sz))*prod
    # the contribution to L of level k is (-1)^k * raw[k]
    return {k:((-1)**k)*raw[k] for k in raw}, base**n

# ---- ROUTE 104: Riesz product certificate ----
def riesz_hat(D, amps):
    """R̂(n) for R=∏_j(1+a_j cos 2π m_j τ): enumerate ε∈{-1,0,1}^|D|, n=Σ ε m, weight ∏(a/2)^|ε| sign."""
    Rhat=defaultdict(float)
    halves=[a/2 for a in amps]
    for eps in itertools.product([-1,0,1], repeat=len(D)):
        n=sum(e*m for e,m in zip(eps,D))
        w=1.0
        for e,hf in zip(eps,halves):
            if e!=0: w*=hf  # |e|=1; sign carried by e in n, magnitude a/2
        # the COEFFICIENT sign: cos contributes (a/2)(e(m)+e(-m)); the coeff of e(n) is +(a/2) for each ±,
        # so the weight is +∏(a/2) (no extra sign) -- BUT with (1 - a cos) we'd flip. We use (1 + a cos);
        # to get NEGATIVE level-1 coeffs use a_j < 0.
        Rhat[n]+=w
    return Rhat

def MR_ratio(S, D, amps, Kmax=14):
    Rhat=riesz_hat(D, amps)
    MR=0.0
    for v in S:
        for k in range(-Kmax, Kmax+1):
            sk = s(k) if k!=0 else 1/7
            MR += sk*Rhat.get(v*k, 0.0)
    R0=Rhat.get(0,0.0)
    return MR, R0, (MR/R0 if R0 else float('inf'))

cores = {
  "extremizer {1..13}\\{6}∪{56}": sorted(set([x for x in range(1,14) if x!=6]+[56])),
  "dilated d=1 {1..12}∪{14}":     sorted(set(list(range(1,13))+[14])),
  "tight 14*{1..13} (L=0 ref)":   [14*i for i in range(1,14)],
}

print("="*74)
print("ROUTE 097: LEVEL DECOMPOSITION  L = (6/7)^13 + Σ_k (contribution of |T|=k)")
print("="*74)
for name,S in cores.items():
    if len(S)!=13: continue
    lev,main = level_masses(S, 6, 4)
    Lm = lonely_measure(S, 11000)
    pos2 = lev.get(2,0); m3 = lev.get(3,0); m4 = lev.get(4,0)
    ge3 = sum(v for k,v in lev.items() if k>=3)
    print(f"  {name}:")
    print(f"     main (6/7)^13 = {main:.4f}; level-2 = {pos2:+.4f}; level-3 = {m3:+.4f}; level-4 = {m4:+.4f}")
    print(f"     Σ_{{|T|≥3}} contribution = {ge3:+.4f}   |Σ≥3| {'>' if abs(ge3)>main else '<'} main "
          f"=> OPEN-Q-097 literal target {'FALSE' if abs(ge3)>main else 'OK'}")
    print(f"     partial L (main+Σlevels, B=6) = {main+sum(lev.values()):.4f}  vs true L = {Lm:.4f}")

print("\n" + "="*74)
print("ROUTE 104: RIESZ-PRODUCT CERTIFICATE — optimize ∫M·R / ∫R  (target < 1 ⟹ loose)")
print("="*74)
Sx = sorted(set([x for x in range(1,14) if x!=6]+[56]))
print(f"  target core {Sx} (L≈0.0056)")
# generators = the speeds; scan uniform amplitude a (negative => negative level-1 coeffs)
print("  D = speeds, uniform amplitude a (a<0 makes level-1 R̂ negative, the wanted sign):")
best=(2.0,None)
for a in [-0.2,-0.4,-0.6,-0.8,-1.0]:
    MR,R0,ratio=MR_ratio(Sx, Sx, [a]*len(Sx), Kmax=14)
    print(f"     a={a:+.1f}: ∫M·R={MR:.4f}, ∫R={R0:.4f}, ratio={ratio:.4f}  {'< 1  CERTIFICATE!' if ratio<1 else ''}")
    if ratio<best[0]: best=(ratio,a)
print(f"  best uniform-speed ratio = {best[0]:.4f} at a={best[1]}")
# try DISSOCIATED subset of speeds (a Sidon-like subset) to avoid level-0 inflation
print("\n  D = a dissociated/Sidon-like subset (reduce R̂(0) inflation from generator relations):")
for sub in [[1,2,4,8,13,56], [1,3,5,11,12,56], [3,5,9,11,13,56]]:
    sub=[x for x in sub if x in Sx]
    for a in [-0.7,-0.9]:
        MR,R0,ratio=MR_ratio(Sx, sub, [a]*len(sub), Kmax=14)
        print(f"     D={sub}, a={a}: ∫M·R={MR:.4f}, ∫R={R0:.4f}, ratio={ratio:.4f}  {'< 1 CERT!' if ratio<1 else ''}")
print("\n  NOTE: a < 1 certificate at the extremizer would prove looseness of THAT core; a UNIFORM")
print("  family certificate (one construction, all primitive mult-of-14 S) is OPEN-Q-104.")
print("\nDONE.")
