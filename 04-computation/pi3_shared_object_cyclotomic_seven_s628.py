import math, itertools, random
# ===== 1) the forbidden values are CYCLOTOMIC at 2: 7 = Phi_3(2), 21 = 3*Phi_3(2) =====
def Phi3(x): return x*x+x+1   # 3rd cyclotomic polynomial X^2+X+1 (roots = primitive cube roots of unity)
print("=== forbidden H-values are cyclotomic at the cube-root (pi/3) angle ===")
print(f"  Phi_3(2) = 2^2+2+1 = {Phi3(2)}  = the forbidden atom 7")
print(f"  3*Phi_3(2) = {3*Phi3(2)} = the forbidden 21 = Phi_2(2)*Phi_3(2) (Phi_2(2)=3)")
print(f"  Phi_3 roots = e^{{+-2pi i/3}} (primitive cube roots) -> angle 2pi/3; 7 = norm-2 evaluation")
print(f"  cf Jacobsthal J(4)=21=3*7=3*Phi_3(2) (the path-UDG hits the forbidden value at the cyclotomic point)")
print(f"  Mersenne-like: Phi_3(2)=(2^3-1)/(2-1)=7; generally (q^3-1)/(q-1)=q^2+q+1")

# ===== 2) the 0.014 excess = Cl_2(pi/3) - 1 (Clausen = ideal regular tetrahedron volume) =====
def Cl2(theta, N=200000):
    # Clausen Cl_2(theta) = -int_0^theta log|2 sin(t/2)| dt  (Lobachevsky), max at pi/3
    s=0.0; dt=theta/N
    for k in range(1,N):
        t=k*dt; s+=-math.log(abs(2*math.sin(t/2)))
    return s*dt
cl=Cl2(math.pi/3)
print(f"\n=== the 0.014 exponent excess = Cl_2(pi/3) - 1 (Clausen max = ideal regular tetrahedron volume) ===")
print(f"  Cl_2(pi/3) = {cl:.6f}   (Lobachevsky/Catalan-type const; ideal-tetrahedron vol = 3*Lobachevsky(pi/3))")
print(f"  excess Cl_2(pi/3) - 1 = {cl-1:.6f}  ~ the 0.014 unit-distance / SC-shape growth excess (HYP-2184)")

# ===== 3) the pi/3 shared-object web =====
print("\n=== the pi/3 (Eisenstein zeta_6 / 60 deg) shared object across all threads ===")
for thread,manifest in [
 ("unit distance","triangular/EISENSTEIN lattice Z[omega]=Q(sqrt-3); 60deg=pi/3; chord=1 <=> dZ=1/6 (S623)"),
 ("tournament forbidden H","7=Phi_3(2), 21=3Phi_3(2); Lee-Yang/3-cycle skew-eigenvalues at +-2pi/3 (cube roots)"),
 ("LRC","gap delta=1/6 = the hexagonal chord; n=14 doubling, the 2-3 seam"),
 ("Collatz","x3 vs 2-adic; the cube-root / Phi_3 angle in 2^K=3^L resonance"),
 ("0.014 exponent","Cl_2(pi/3)-1 = ideal-tetrahedron volume excess (Dehn/scissors-congruence signature)"),
 ("Jacobsthal chain","J(4)=21=3*Phi_3(2): the path-UDG hits forbidden at the cyclotomic value"),
]:
    print(f"  {thread:24}: {manifest}")

# ===== 4) SC (self-complementary) tournaments: H-values, the alpha_2=1 / norm-1 family =====
def make_beats(n,bits):
    b=[[False]*n for _ in range(n)]; idx=0
    for i in range(n):
        for j in range(i+1,n):
            if (bits>>idx)&1: b[i][j]=True
            else: b[j][i]=True
            idx+=1
    return b
def ham_paths(n,beats):
    full=(1<<n)-1; dp=[[0]*n for _ in range(1<<n)]
    for v in range(n): dp[1<<v][v]=1
    for mask in range(1<<n):
        for last in range(n):
            c=dp[mask][last]
            if not c: continue
            for nx in range(n):
                if not (mask>>nx)&1 and beats[last][nx]: dp[mask|(1<<nx)][nx]+=c
    return sum(dp[full][v] for v in range(n))
def is_SC(n,beats):
    # self-complementary: exists permutation p with beats[i][j] == reverse[p_i][p_j]=beats[p_j][p_i]
    for perm in itertools.permutations(range(n)):
        ok=True
        for i in range(n):
            for j in range(n):
                if i!=j and beats[i][j]!=(not beats[perm[i]][perm[j]] if perm[i]!=perm[j] else False):
                    # complement: arc reversed
                    if beats[i][j]==beats[perm[i]][perm[j]]: ok=False;break
            if not ok: break
        if ok: return True
    return False
print("\n=== SC (self-complementary, T iso T^op = sigma/CM-conjugation fixed) tournament H-values ===")
for n in [4,5]:
    m=n*(n-1)//2; Hs=set()
    for bits in range(1<<m):
        beats=make_beats(n,bits)
        if is_SC(n,beats): Hs.add(ham_paths(n,beats))
    print(f"  n={n} (n mod 4 = {n%4}): SC tournament H-values = {sorted(Hs)}")
print("  SC = fixed points of the complement involution tau (=CM conjugation c = sigma); the alpha_2=1 'norm-1' family (roots rho_1 rho_2 = 1 = |beta|^2).")
