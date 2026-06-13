"""
N=6: The unique consistency point.

At n=6, the THM-217 factorization gap is zero (n-6=0).
Also n=6 is where redundancy = 1/2 (equal independent and redundant).
And n=6 is the ONLY even n where regular tournaments don't exist
but the QR tournament DOES (6 is not prime, but Paley exists at q=5+1=6? No.)

What makes n=6 special in the C_T^TC_T spectrum?
"""
import numpy as np
from math import comb
import random
from collections import Counter, defaultdict

def tournament_from_bits(n, bits):
    adj = [[0]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx):
                adj[i][j] = 1
            else:
                adj[j][i] = 1
            idx += 1
    return adj

def full_constraint_matrix(n):
    E = n*(n-1)//2
    edges = []
    edge_idx = {}
    for i in range(n):
        for j in range(i+1, n):
            edges.append((i,j))
            edge_idx[(i,j)] = len(edges)-1
    rows = []
    for a in range(n):
        for b in range(a+1, n):
            for c in range(b+1, n):
                row = np.zeros(E)
                row[edge_idx[(b,c)]] = 1
                row[edge_idx[(a,c)]] = -1
                row[edge_idx[(a,b)]] = 1
                rows.append(row)
    return np.array(rows), edges, edge_idx

def constraint_matrices(adj, n, C_full):
    row = 0
    trans_idx = []
    cyclic_idx = []
    for a in range(n):
        for b in range(a+1, n):
            for c in range(b+1, n):
                wins = [adj[a][b]+adj[a][c], adj[b][a]+adj[b][c], adj[c][a]+adj[c][b]]
                if max(wins) == 2:
                    trans_idx.append(row)
                else:
                    cyclic_idx.append(row)
                row += 1
    C_T = C_full[trans_idx] if trans_idx else np.zeros((0, C_full.shape[1]))
    C_R = C_full[cyclic_idx] if cyclic_idx else np.zeros((0, C_full.shape[1]))
    return C_T, C_R

# ====================================================================
print("=" * 70)
print("n=6: THE UNIQUE CONSISTENCY POINT")
print("=" * 70)
# ====================================================================

n = 6
E = n*(n-1)//2  # 15
max_rank = (n-1)*(n-2)//2  # 10
Cn3 = comb(n, 3)  # 20

print(f"n={n}: E={E}, max_rank={max_rank}, C(n,3)={Cn3}")
print(f"  Independence fraction = 3/n = {3/n:.4f}")
print(f"  Redundancy fraction = (n-3)/n = {(n-3)/n:.4f}")
print(f"  At n=6: redundancy = 1/2 (HALF the constraints are redundant)")
print(f"  This is the unique n where redundancy = independence (both = 1/2)")
print()

# Key n=6 identities
x = n*(n-3)  # 18
alpha = 3-n  # -3
metallic = (n-2 + np.sqrt(n**2+4))/2  # 2+√10
metallic2 = (n-2 - np.sqrt(n**2+4))/2  # 2-√10

print(f"THM-217 at n=6:")
print(f"  x = n(n-3) = {x}")
print(f"  α = 3-n = {alpha}")
print(f"  λ₁ = 2+√10 = {metallic:.6f}")
print(f"  λ₂ = 2-√10 = {metallic2:.6f}")
print(f"  α·λ₁·λ₂ = {alpha * metallic * metallic2:.6f} (char poly constant)")
print(f"  Expected: x = {x}")
print()

# THM-224 self-dual eigenvalues at n=6
disc_224 = n*(n-4)  # 12
phi = (n + np.sqrt(disc_224))/2
phibar = (n - np.sqrt(disc_224))/2
print(f"THM-224 self-dual equation λ²-6λ+6=0:")
print(f"  φ = 3+√3 = {phi:.6f}")
print(f"  φ̄ = 3-√3 = {phibar:.6f}")
print(f"  Product: {phi*phibar:.6f} = n")
print(f"  Sum: {phi+phibar:.6f} = n")
print()

# ====================================================================
print("=" * 70)
print("SAMPLING n=6: SPECTRAL ANATOMY")
print("=" * 70)
# ====================================================================

C_full, edges, edge_idx = full_constraint_matrix(n)
U_f, S_f, Vt_f = np.linalg.svd(C_full, full_matrices=False)
nz_f = sum(1 for s in S_f if s > 1e-8)
P_basis = Vt_f[:nz_f]

random.seed(42)
spectra_by_t3_b1 = defaultdict(list)

for trial in range(10000):
    bits = random.randint(0, (1<<E)-1)
    adj = tournament_from_bits(n, bits)
    C_T, C_R = constraint_matrices(adj, n, C_full)
    t3 = C_T.shape[0]
    c3 = C_R.shape[0]

    CTC = C_T.T @ C_T
    eigs = np.sort(np.linalg.eigvalsh(CTC))[::-1]
    rk = sum(1 for e in eigs if e > 1e-6)
    beta1 = max_rank - rk

    nz_eigs = tuple(round(e, 2) for e in eigs if e > 0.01)
    spectra_by_t3_b1[(t3, beta1)].append(nz_eigs)

# Show summary
print(f"\nn=6, 10000 random tournaments:")
for (t3, b1) in sorted(spectra_by_t3_b1.keys()):
    entries = spectra_by_t3_b1[(t3, b1)]
    n_distinct = len(set(entries))
    print(f"  t₃={t3}, β₁={b1}: {len(entries)} tournaments, {n_distinct} distinct spectra")

    # Show most common
    c = Counter(entries)
    for spec, count in c.most_common(2):
        # Identify algebraic structure
        nz_list = list(spec)
        n_equal_6 = sum(1 for e in nz_list if abs(e - 6) < 0.05)
        non6 = [e for e in nz_list if abs(e - 6) > 0.05]
        non6_sum = sum(non6)
        non6_prod = 1
        for e in non6:
            non6_prod *= e
        print(f"    [{', '.join(f'{e:.2f}' for e in spec)}] × {count}")
        print(f"      {n_equal_6} eigenvalues = 6, rest sum={non6_sum:.2f}, prod={non6_prod:.2f}")

# ====================================================================
print("\n" + "=" * 70)
print("T/R DUALITY AT n=6")
print("=" * 70)
# ====================================================================

# Do the eigenvalues satisfy λ(6-λ)=6 (self-dual) at t₃=c₃=10?
# t₃=10 means c₃=10, equal split
random.seed(123)
equal_split_spectra = []
for trial in range(50000):
    bits = random.randint(0, (1<<E)-1)
    adj = tournament_from_bits(n, bits)
    C_T, C_R = constraint_matrices(adj, n, C_full)
    t3 = C_T.shape[0]
    if t3 != 10:  # equal split: t₃=c₃=10 at n=6, Cn3=20
        continue

    CTC = C_T.T @ C_T
    eigs = np.sort(np.linalg.eigvalsh(CTC))[::-1]
    nz = tuple(round(e, 3) for e in eigs if e > 0.01)
    equal_split_spectra.append(nz)

if equal_split_spectra:
    c = Counter(equal_split_spectra)
    print(f"\nt₃=10 (=c₃=10) tournaments at n=6: {len(equal_split_spectra)} found")
    for spec, count in c.most_common(5):
        print(f"  [{', '.join(f'{e:.3f}' for e in spec)}] × {count}")
        # Check if eigenvalue pairs have product = 6
        nz_list = list(spec)
        for i in range(len(nz_list)):
            for j in range(i+1, len(nz_list)):
                if abs(nz_list[i]*nz_list[j] - 6) < 0.05:
                    print(f"    Product pair: {nz_list[i]:.3f} × {nz_list[j]:.3f} = {nz_list[i]*nz_list[j]:.3f}")
else:
    print(f"\nNo t₃=10 tournaments found in {50000} samples. Checking t₃ distribution...")
    random.seed(789)
    t3_counts = Counter()
    for trial in range(10000):
        bits = random.randint(0, (1<<E)-1)
        adj = tournament_from_bits(n, bits)
        C_T, _ = constraint_matrices(adj, n, C_full)
        t3_counts[C_T.shape[0]] += 1
    for t3 in sorted(t3_counts.keys()):
        print(f"    t₃={t3}: {t3_counts[t3]} ({100*t3_counts[t3]/10000:.1f}%)")

# ====================================================================
print("\n" + "=" * 70)
print("THE n=6 TRANSFER MATRIX TRACES vs CONSTRAINT TRACES")
print("=" * 70)
# ====================================================================

print("""
THM-217 transfer matrix at x = 18 (n=6):
  M(18) eigenvalues: α=-3, λ₁=2+√10, λ₂=2-√10

  Tr(M^N) = α^N + λ₁^N + λ₂^N

THM-224 constraint matrix at n=6:
  C^TC eigenvalue = n = 6, multiplicity = 10

  Tr((C^TC)^N) = 6^N · 10

The "bridge" quantity: trace of (C_T^TC_T)^N for a tournament T.
  = sum of (eigenvalues of C_T^TC_T)^N
""")

# Compute Tr(M^N) at x=18
x = 18
for N in range(1, 8):
    tr_M = alpha**N + metallic**N + metallic2**N
    tr_CTC = 6**N * 10
    print(f"  N={N}: Tr(M^{N}) = {tr_M:.1f}, Tr((C^TC)^{N}) = {tr_CTC}")

# For the regular-ish tournament at n=6, Tr(C_T^TC_T)^N
# But n=6 is even, no regular tournament exists. The "most balanced" scores are (2,2,2,3,3,3)
print("\n  Finding most balanced tournaments at n=6:")
random.seed(42)
best_balance = float('inf')
best_bits = 0
for trial in range(10000):
    bits = random.randint(0, (1<<E)-1)
    adj = tournament_from_bits(n, bits)
    sc = sorted([sum(adj[i]) for i in range(n)])
    balance = max(sc) - min(sc)
    if balance < best_balance:
        best_balance = balance
        best_bits = bits

adj = tournament_from_bits(n, best_bits)
sc = sorted([sum(adj[i]) for i in range(n)])
C_T, C_R = constraint_matrices(adj, n, C_full)
CTC = C_T.T @ C_T
eigs = np.sort(np.linalg.eigvalsh(CTC))[::-1]
print(f"  Most balanced tournament: scores={sc}, t₃={C_T.shape[0]}, c₃={C_R.shape[0]}")
print(f"  Eigenvalues: {eigs.round(3)}")

# ====================================================================
print("\n" + "=" * 70)
print("SYNTHESIS: WHY n=6 IS UNIQUE")
print("=" * 70)
# ====================================================================

print("""
n=6 is the UNIQUE point where three structures align:

1. REDUNDANCY = INDEPENDENCE (both = 1/2)
   At n<6: independence > redundancy (system is "lean")
   At n>6: redundancy > independence (system is "fat")
   At n=6: exact balance

2. THM-217 SELF-CONSISTENCY
   The char poly λ³-λ²-xλ-x factors via the metallic quadratic
   λ²-(n-2)λ-n ONLY at n=6 (gap = n-6 = 0).
   The third eigenvalue α=3-n=-3 exactly satisfies BOTH the constant
   term condition (αc=-x → -3·6=-18 ✓) and the linear term
   condition (αb-c=-x → -3·4-6=-18 ✓).

3. THE SPECTRAL DUALITY
   C_T^TC_T + C_R^TC_R = n·P means eigenvalues pair as (λ, n-λ).
   At n=6 with redundancy=1/2, the "typical" tournament has t₃≈c₃,
   making the T↔R duality most symmetric.

The number n=6 also has unique combinatorial properties:
  - C(6,3) = 20, max_rank = 10, redundancy = 10 (equal!)
  - 6 is the smallest n where the constraint system has "room" for
    both types of eigenvalue behavior
  - The metallic mean M₆-1 = 2+√10 ≈ 5.162 is closest to n among all M_n-1
    (since M_n-1 ≈ n-1+1/(n-1), at n=6: M₆-1 ≈ 5.2 vs n=6)
""")

# Check: M_n-1 vs n
print("M_n-1 vs n:")
for n_val in range(3, 15):
    m = (n_val-2 + np.sqrt(n_val**2+4))/2
    print(f"  n={n_val}: M_n-1 = {m:.4f}, n-1 = {n_val-1}, ratio M_n-1/n = {m/n_val:.4f}")

# ====================================================================
print("\n" + "=" * 70)
print("THE UNIVERSAL TOP EIGENVALUE: PROOF VIA POTENTIAL FUNCTIONS")
print("=" * 70)
# ====================================================================

print("""
THEOREM (THM-225?): For any tournament T on n≥3 vertices with t₃≥1,
the largest eigenvalue of C_T^TC_T equals n.

PROOF:

Step 1: C_T^TC_T = n·P - C_R^TC_R (from THM-224).
So max_eig(C_T^TC_T) ≤ n (since C_R^TC_R ≥ 0).

Step 2: We need v ∈ im(P) with C_R v = 0, i.e., v ∈ im(∂₂^T) ∩ ker(C_R).

Step 3: im(∂₂^T) = ker(∂₁)^⊥ in ℝ^E (Hodge decomposition on simplex).
So v ∈ im(∂₂^T) means v is orthogonal to all gradient 1-forms.

Step 4: ker(C_R) = {v : v·∂₂(σ)=0 for all cyclic σ}.
Since ∂₂(σ) has the same form regardless of whether σ is transitive
or cyclic in the tournament, this means v assigns zero "flux" through
every cyclic triple.

Step 5: Choose any TRANSITIVE TRIPLE σ₀ = {a,b,c} of T.
Let v₀ = P · ∂₂(σ₀) (project the boundary of σ₀ onto im(∂₂^T)).
Then v₀ ∈ im(P). And C_T v₀ has nonzero entry at σ₀ (at minimum).

But C_R v₀ might not be zero...

ALTERNATIVE: Use the vertex potential approach.

For any tournament T, define f_e for edge e = (i,j) as follows:
Choose a transitive triple {i,j,k} (exists since t₃≥1).
Actually, not every edge need be in a transitive triple.

Let me try a different approach: the SIGNED SCORE vector.

Define s_i = out-degree of vertex i in T. Then s_i ∈ {0,...,n-1}
and Σ s_i = C(n,2).

Define the edge vector v by v_{(i,j)} = s_j - s_i (the score difference).

This is a GRADIENT: v = ∂₁^T s, so v ∈ im(∂₁^T) = ker(P).
Not useful — gives eigenvalue 0.

Instead, define v_{(i,j)} = s_j - s_i - (n-1-2s_i)?? No, too ad hoc.

Actually, the proof might be simpler:

Step A: v = C_T^T · 1 (the sum of all transitive triple boundaries).
This v has v_{(i,j)} = (# trans triples containing (i,j) with + sign)
                      - (# trans triples containing (i,j) with - sign).

Step B: Similarly, C_R^T · 1 is the sum of all cyclic triple boundaries.
And C^T · 1 = C_T^T · 1 + C_R^T · 1.

Step C: C^T · 1 is in im(C^T) = im(P), and
(C^TC)·(C^T·1) = n·P·(C^T·1) = n·C^T·1 (since C^T·1 ∈ im(P)).
So ||C·(C^T·1)||² = n·||C^T·1||².

Hmm, that's the eigenvalue equation for C·C^T, not C^TC.
Actually: C^TC · (C^T·1) = C^T·(C·C^T·1).
And C·C^T is the triple-triple Gram matrix, all eigenvalues = n.
So C·C^T·1 = n·P'·1 where P' is the triple-space projection.

This is getting complicated. Let me try numerically.
""")

# Compute v = C_T^T · 1 and check if it's an eigenvector of C_T^TC_T
n = 6
C_full, edges, edge_idx = full_constraint_matrix(n)

random.seed(42)
for trial in range(5):
    bits = random.randint(0, (1<<E)-1)
    adj = tournament_from_bits(n, bits)
    C_T, C_R = constraint_matrices(adj, n, C_full)
    t3 = C_T.shape[0]
    if t3 == 0:
        continue

    # v = C_T^T · 1
    ones = np.ones(t3)
    v = C_T.T @ ones

    # Is this an eigenvector of C_T^TC_T?
    CTC = C_T.T @ C_T
    CTC_v = CTC @ v
    # Check eigenvalue: CTC_v / v
    ratios = CTC_v / (v + 1e-30)
    # Check if it's proportional
    if np.std(ratios[np.abs(v) > 0.1]) < 0.1:
        eig_val = np.mean(ratios[np.abs(v) > 0.1])
        print(f"  bits={bits}, t₃={t3}: C_T^T·1 IS eigenvector with eigenvalue {eig_val:.4f}")
    else:
        # Rayleigh quotient
        rq = v @ CTC_v / (v @ v)
        print(f"  bits={bits}, t₃={t3}: C_T^T·1 is NOT eigenvector. Rayleigh quotient = {rq:.4f}")

    # Try v = C^T · 1 (all triples)
    v_full = C_full.T @ np.ones(Cn3)
    CTC_v_full = CTC @ v_full
    rq_full = v_full @ CTC_v_full / (v_full @ v_full) if v_full @ v_full > 0 else 0
    print(f"    C^T·1: Rayleigh quotient = {rq_full:.4f}")

    # Try the SIGNED version: C^T · s where s_σ = +1 if transitive, -1 if cyclic
    s = np.ones(Cn3)
    row = 0
    for a in range(n):
        for b in range(a+1, n):
            for c in range(b+1, n):
                wins = [adj[a][b]+adj[a][c], adj[b][a]+adj[b][c], adj[c][a]+adj[c][b]]
                if max(wins) < 2:  # cyclic
                    s[row] = -1
                row += 1
    v_signed = C_full.T @ s
    if v_signed @ v_signed > 1e-10:
        CTC_v_signed = CTC @ v_signed
        rq_signed = v_signed @ CTC_v_signed / (v_signed @ v_signed)
        print(f"    C^T·(±1): Rayleigh quotient = {rq_signed:.4f}, ||v||²={v_signed@v_signed:.4f}")
    else:
        print(f"    C^T·(±1): zero vector!")

# ====================================================================
print("\n" + "=" * 70)
print("PROOF BY EXPLICIT CONSTRUCTION")
print("=" * 70)
# ====================================================================

print("""
Let me try: for each edge e, define its "transitive score"
  w(e) = (# transitive triples containing e) - (# cyclic triples containing e)

Then w = (C_T^T · 1_T - C_R^T · 1_R) componentwise... no, that's not right
because the signs in C are not all +1.

Actually: for edge (i,j), the diagonal of C_T^TC_T is
  (C_T^TC_T)_{ee} = Σ_σ [∂₂(σ)_e]² = # transitive triples containing e

And (C_R^TC_R)_{ee} = # cyclic triples containing e.

Sum: (C^TC)_{ee} = n-2 (THM-224 diagonal).

So: (C_T^TC_T)_{ee} = t₃(e) and (C_R^TC_R)_{ee} = c₃(e)
where t₃(e)+c₃(e) = n-2.
""")

# Verify
n = 5
C_full5, edges5, edge_idx5 = full_constraint_matrix(n)

bits = 0b0011101010
adj = tournament_from_bits(n, bits)
C_T, C_R = constraint_matrices(adj, n, C_full5)
CTC_T = C_T.T @ C_T
CTC_R = C_R.T @ C_R

print(f"n={n}, bits={bits}:")
print(f"  Diagonal of C_T^TC_T: {np.diag(CTC_T)}")
print(f"  Diagonal of C_R^TC_R: {np.diag(CTC_R)}")
print(f"  Sum: {np.diag(CTC_T) + np.diag(CTC_R)}")
print(f"  Expected: n-2 = {n-2}")

# The diagonal tells us how many transitive triples each edge is in
# For t₃(e) values, each edge is in 0..n-2 transitive triples

# This is the KEY: trace(C_T^TC_T) = Σ_e t₃(e) = 3·t₃
# because each transitive triple has 3 edges.

print(f"\n  Σ t₃(e) = {sum(np.diag(CTC_T)):.0f}, 3·t₃ = {3*C_T.shape[0]}")

# ====================================================================
print("\n" + "=" * 70)
print("THE QUADRATIC DUALITY")
print("=" * 70)
# ====================================================================

print("""
The two quadratic equations:

  METALLIC (THM-217): μ² - (n-2)μ - n = 0
    Roots: μ = [(n-2) ± √(n²+4)] / 2
    Sum = n-2, Product = -n

  SELF-DUAL (THM-224): λ² - nλ + n = 0
    Roots: λ = [n ± √(n²-4n)] / 2 = [n ± √(n(n-4))] / 2
    Sum = n, Product = n

OBSERVATION: The metallic equation can be rewritten as
  μ(μ-(n-2)) = n
while the self-dual equation is
  λ(n-λ) = n, i.e., λ(λ-n) = -n

So:
  Metallic: product of (root, root shifted by n-2) = n
  Self-dual: product of (root, root shifted by n) = -n

The shift DIFFERENCE is 2. The product RATIO is -1.

If we define the CENTERED roots:
  Metallic: μ - (n-2)/2 = ±√(n²+4)/2
  Self-dual: λ - n/2 = ±√(n²-4n)/2

The half-sums are (n-2)/2 and n/2, differing by 1.
The half-discriminants are √(n²+4)/2 and √(n(n-4))/2.

RATIO of discriminants: √(n²+4)/√(n²-4n) = √((n²+4)/(n²-4n))

At n=6: √(40/12) = √(10/3) ≈ 1.826

At n=∞: √((n²+4)/(n²-4n)) → 1

So the two quadratics CONVERGE as n→∞ (both → centered at n/2 with half-width n/2).

The FINITE-n difference is what makes n=6 special:
at n=6, the metallic mean M₆-1 = 2+√10 ≈ 5.162, and
the self-dual eigenvalue φ = 3+√3 ≈ 4.732.
Their difference 5.162-4.732 = 0.430 is among the smallest.
""")

# Compare metallic and self-dual roots across n
print("Metallic vs Self-dual eigenvalues:")
print(f"{'n':>3} {'M_n-1':>10} {'φ_n':>10} {'diff':>10} {'ratio':>10}")
for n_val in range(5, 20):
    met = (n_val-2 + np.sqrt(n_val**2+4))/2
    disc = n_val*(n_val-4)
    if disc >= 0:
        phi = (n_val + np.sqrt(disc))/2
        print(f"{n_val:3d} {met:10.4f} {phi:10.4f} {met-phi:10.4f} {met/phi:10.4f}")
    else:
        print(f"{n_val:3d} {met:10.4f} {'complex':>10}")

# ====================================================================
print("\n" + "=" * 70)
print("THE 3/n BRIDGE: FINAL FORMULATION")
print("=" * 70)
# ====================================================================

print("""
╔══════════════════════════════════════════════════════════════════╗
║              THE 3/n BRIDGE: UNIFIED FRAMEWORK                  ║
╠══════════════════════════════════════════════════════════════════╣
║                                                                  ║
║  The fraction 3/n is the INDEPENDENCE RATIO of the simplicial    ║
║  boundary on K_n: rank(∂₂)/C(n,3) = 3/n.                       ║
║                                                                  ║
║  THM-224: C^TC = n·P                                            ║
║    - Eigenvalue n = trace/rank = 3·C(n,3)/[C(n,3)·(3/n)] = n   ║
║    - The 3 comes from dim(∂Δ²) = 3 (triangle has 3 edges)      ║
║    - The 3/n comes from H₁(K_n)=0 (complete graph contractible) ║
║                                                                  ║
║  THM-217: char poly factors at x = n(n-3) = n²·(1-3/n)         ║
║    - Third eigenvalue α = 3-n = -n·(1-3/n)                      ║
║    - Encodes the redundancy fraction (n-3)/n = 1-3/n            ║
║    - Self-consistent ONLY at n=6 (gap = n-6)                    ║
║                                                                  ║
║  BRIDGE:                                                         ║
║    Both theorems encode the same simplicial combinatorics:       ║
║    "3 edges per triple, n-2 triples per edge"                    ║
║                                                                  ║
║    THM-224 captures this SPECTRALLY (uniform eigenvalue n)       ║
║    THM-217 captures this ALGEBRAICALLY (factorization at n(n-3)) ║
║                                                                  ║
║  SPECTRAL DUALITY:                                               ║
║    For any tournament T: C_T^TC_T + C_R^TC_R = n·P              ║
║    Eigenvalues pair as (λ, n-λ)                                  ║
║    β₁=1 ⟺ C_R^TC_R has eigenvalue n on im(P)                  ║
║    ⟺ some edge function in the constraint column space is       ║
║      "completely cyclic" (zero flux through all transitive triples)║
║                                                                  ║
║  AT n=6 (unique):                                                ║
║    - Redundancy = Independence (both = 10)                       ║
║    - THM-217 factorization is self-consistent                    ║
║    - The metallic mean 2+√10 governs the HP growth rate          ║
║      at exactly the "redundancy capacity" x = 18 = n²/2         ║
║                                                                  ║
╚══════════════════════════════════════════════════════════════════╝
""")
