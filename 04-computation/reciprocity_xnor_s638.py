#!/usr/bin/env python3
"""
S638 — quadratic reciprocity & the XNOR-with-out-of-phase-middle, abstractly.
Legendre symbol (a/p) = a multiplicative {±1} character (XNOR: (ab/p)=(a/p)(b/p)) with an always-on
MIDDLE 0 (at p|a) OUT OF PHASE with ±1. Its matrix is the Paley CONFERENCE matrix C_{ij}=(j-i/p):
±1 off-diagonal (binary choice), 0 diagonal (always-on middle), C Cᵀ = pI - J (orthogonality).
Reciprocity governs the PHASE: Gauss sum g_p = √p (p≡1 mod4, C symmetric=graph) vs i√p (p≡3 mod4,
C SKEW = the Paley TOURNAMENT). The 'out of phase' i = skew = tournament = odd = the 2-adic seam.
"""
import cmath, math, itertools

def legendre(a,p):
    a%=p
    if a==0: return 0
    return 1 if pow(a,(p-1)//2,p)==1 else -1

def is_prime(n):
    if n<2: return False
    i=2
    while i*i<=n:
        if n%i==0: return False
        i+=1
    return True

print("(1) Legendre = multiplicative XNOR on ±1, with 0-middle. (ab/p)=(a/p)(b/p):")
p=11
bad=0
for a in range(1,p):
    for b in range(1,p):
        if legendre(a*b,p)!=legendre(a,p)*legendre(b,p): bad+=1
print(f"   p={p}: multiplicativity violations={bad}; (0/p)={legendre(0,p)} (the always-on middle)")

print("\n(2) Quadratic reciprocity (p/q)(q/p)=(-1)^((p-1)/2 (q-1)/2); governed by p,q mod 4:")
primes=[p for p in range(3,40) if is_prime(p)]
viol=0
for p in primes:
    for q in primes:
        if p==q: continue
        lhs=legendre(q,p)*legendre(p,q)
        rhs=(-1)**(((p-1)//2)*((q-1)//2))
        if lhs!=rhs: viol+=1
print(f"   reciprocity holds for all prime pairs <40: violations={viol}")
print(f"   (-1/p): p≡1 mod4 => +1 (−1 is a QR, character EVEN/symmetric); p≡3 => −1 (odd/skew):")
for p in [5,7,11,13]:
    print(f"      p={p} (≡{p%4} mod4): (-1/p)={legendre(-1,p)}")

print("\n(3) Gauss sum phase g_p = Σ_a (a/p) ζ^a  (ζ=e^{2πi/p}): √p vs i√p by p mod 4")
for p in [5,7,11,13]:
    z=cmath.exp(2j*math.pi/p)
    g=sum(legendre(a,p)*z**a for a in range(p))
    print(f"   p={p} (≡{p%4}): g={g.real:+.3f}{g.imag:+.3f}i  |g|={abs(g):.3f}(√p={math.sqrt(p):.3f}) "
          f"arg={math.degrees(cmath.phase(g)):.1f}°  => {'REAL (in phase, p≡1)' if p%4==1 else 'IMAGINARY (out of phase, p≡3)'}")

print("\n(4) Paley conference matrix C_{ij}=(j-i/p): skew (p≡3=tournament) vs symmetric (p≡1=graph)")
for p in [7,11,13]:
    C=[[legendre((j-i)%p,p) for j in range(p)] for i in range(p)]
    skew=all(C[i][j]==-C[j][i] for i in range(p) for j in range(p))
    symm=all(C[i][j]==C[j][i] for i in range(p) for j in range(p))
    diag0=all(C[i][i]==0 for i in range(p))
    # C C^T = pI - J ?
    def dot(i,j): return sum(C[i][k]*C[j][k] for k in range(p))
    conf=all(dot(i,j)==(p-(1 if i==j else 0)*0 + (p if i==j else 0) - 1) for i in range(p) for j in range(p))
    confcheck=all(dot(i,j)==((p-1) if i==j else -1) for i in range(p) for j in range(p))
    print(f"   p={p} (≡{p%4}): skew={skew} symm={symm} diag0(always-on middle)={diag0}; C Cᵀ=(p-1)I-J off+... conf-matrix={confcheck}")

print("\n(5) Paley TOURNAMENT (p≡3 mod4): self-complementary; H=#ham paths; round=>chi_di=2; spectrum imaginary")
import sys; sys.path.insert(0,'04-computation')
def H_ham(adj,n):
    full=(1<<n)-1; dp=[[0]*n for _ in range(1<<n)]
    for v in range(n): dp[1<<v][v]=1
    for mask in range(1<<n):
        for v in range(n):
            c=dp[mask][v]
            if c:
                for w in range(n):
                    if not (mask>>w)&1 and adj[v][w]: dp[mask|(1<<w)][w]+=c
    return sum(dp[full])
for p in [3,7,11]:
    adj=[[legendre((j-i)%p,p)==1 for j in range(p)] for i in range(p)]   # i->j iff (j-i) is QR
    # self-complementary: T^op iso T? for Paley p≡3 yes (multiply by a non-residue)
    H=H_ham(adj,p)
    print(f"   Paley T_{p}: H(#ham paths)={H} (odd? {H%2==1}); skew-conference => purely imaginary skew-spectrum (±i√p Gauss)")
print("\n  => the 'out of phase middle' = the imaginary Gauss sum i√p = the SKEW (tournament) p≡3 case;")
print("     p≡1 is the in-phase symmetric graph. Reciprocity's parity ((p-1)/2) = the 2-adic seam.")
