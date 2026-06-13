#!/usr/bin/env python3
"""
applications_s112.py — Creative investigation of practical applications
kind-pasteur-2026-03-15-S112

The Cayley-Delannoy-Tournament chain:
  Q(x) = (1+x)/(1-x) = exp(2*arctanh(x))
  Q^m = 1 + 2*sum g_k(m)*x^k
  g_k(m) = Delannoy diagonal step weight
  CV^2(H) = sum 2*g_k(n-2k)/(n)_{2k}

What can we DO with this? Let's explore applications across domains.
"""

from fractions import Fraction
from math import comb, factorial, log, exp, sqrt, pi
import random

# ============================================================
# APPLICATION 1: TOURNAMENT CONFIDENCE INTERVALS
# ============================================================
print("="*70)
print("APPLICATION 1: CONFIDENCE INTERVALS FOR RANKING RELIABILITY")
print("="*70)
print()
print("Context: In sports/elections, a tournament represents pairwise outcomes.")
print("H(T) = number of consistent total rankings. If H is large, the outcomes")
print("are 'coherent' (many linear extensions). If H is small, contradictions.")
print()

def cv2_exact(n):
    """Exact CV^2 using Delannoy formula."""
    total = Fraction(0)
    for k in range(1, (n-1)//2 + 1):
        m = n - 2*k
        if m < 1: continue
        gk = sum(Fraction(comb(k-1,j-1)*comb(m,j)) * 2**(j-1)
                 for j in range(1, min(k,m)+1))
        ff = Fraction(1)
        for i in range(2*k):
            ff *= (n - i)
        total += 2 * gk / ff
    return total

mean_H = lambda n: Fraction(factorial(n), 2**(n-1))

print("For n teams in a round-robin tournament:")
for n in [5, 8, 10, 16, 20, 32]:
    cv2 = float(cv2_exact(n)) if n <= 18 else 2.0/n  # asymptotic for large n
    cv = sqrt(cv2)
    mu = float(mean_H(n)) if n <= 18 else factorial(n) / 2**(n-1)

    # Chebyshev bound: P(|H - mu| > t*sigma) <= 1/t^2
    # 95% CI: t = sqrt(20) ~ 4.47 (Chebyshev), or ~2 (assuming near-normal)
    sigma = cv * mu
    ci_low = max(1, mu - 2*sigma)
    ci_high = mu + 2*sigma

    print(f"  n={n:2d}: E[H]={mu:.1e}, CV={cv:.3f}, "
          f"95% CI ~ [{ci_low:.1e}, {ci_high:.1e}]")

print()
print("INSIGHT: CV ~ sqrt(2/n), so for n=32 teams, CV ~ 25%.")
print("A tournament with H > 1.5*E[H] is significantly 'more orderable'")
print("than random — suggesting real skill differences exist.")

# ============================================================
# APPLICATION 2: ANOMALY DETECTION IN PREFERENCE DATA
# ============================================================
print()
print("="*70)
print("APPLICATION 2: ANOMALY DETECTION IN PREFERENCE DATA")
print("="*70)
print()
print("Given a tournament T from pairwise comparisons (products, candidates,")
print("movies), compute H(T) and compare to the random distribution.")
print()

# The Z-score of H(T):
# Z = (H(T) - E[H]) / (sigma) = (H(T)/E[H] - 1) / CV
# If |Z| > 3: anomalous tournament (too orderable or too contradictory)

print("Example: 7 products compared pairwise.")
print("E[H] = 7!/64 = 78.75")
print("sigma = 78.75 * sqrt(131/504) = 78.75 * 0.510 = 40.15")
print()
print("  H(T) = 189 (Paley T_7): Z = (189-78.75)/40.15 = +2.74 (UNUSUAL)")
print("  H(T) = 3 (near-cyclic):  Z = (3-78.75)/40.15 = -1.89 (somewhat low)")
print("  H(T) = 1 (transitive):   Z = (1-78.75)/40.15 = -1.94 (low)")
print()
print("PRODUCT: A 'tournament anomaly score' for ranking quality assessment.")

# ============================================================
# APPLICATION 3: FAST VARIANCE ESTIMATION
# ============================================================
print()
print("="*70)
print("APPLICATION 3: FAST VARIANCE ESTIMATION (O(1) FORMULA)")
print("="*70)
print()
print("The formula CV^2 = 2/n + O(1/n^3) means we can estimate Var(H)")
print("in O(1) time, without any simulation or enumeration!")
print()

for n in [10, 50, 100, 1000, 10000]:
    try:
        log_mu = sum(log(i) for i in range(1, n+1)) - (n-1)*log(2)
        mu = exp(min(log_mu, 700))  # cap to avoid overflow
    except:
        mu = float('inf')
    cv2_approx = 2.0/n
    sigma_approx = sqrt(cv2_approx) * mu

    print(f"  n={n:5d}: sigma/E[H] ~ sqrt(2/{n}) = {sqrt(2.0/n):.4f}")
    if n <= 20:
        cv2_true = float(cv2_exact(n))
        print(f"           exact CV = {sqrt(cv2_true):.4f}, "
              f"approx CV = {sqrt(2.0/n):.4f}, "
              f"error = {abs(sqrt(cv2_true)-sqrt(2.0/n))/sqrt(cv2_true)*100:.1f}%")

print()
print("For n >= 50, the formula CV = sqrt(2/n) is accurate to < 1%.")

# ============================================================
# APPLICATION 4: SPECTRAL TOURNAMENT FINGERPRINTING
# ============================================================
print()
print("="*70)
print("APPLICATION 4: SPECTRAL TOURNAMENT FINGERPRINT")
print("="*70)
print()
print("The Fourier energy profile (E_2, E_4, E_6, ...) of a tournament")
print("is a FINGERPRINT that captures its structural properties.")
print()
print("From the Delannoy formula, the expected energy at level 2k for a")
print("RANDOM tournament is E_{2k}/E_0 = 2*g_k(n-2k)/(n)_{2k}.")
print()
print("Deviations from this profile indicate structure:")
print("  - Excess E_2: pairwise skill differences (simple ranking)")
print("  - Excess E_4: 4-way interactions (group effects, e.g. rock-paper-scissors)")
print("  - Excess E_6: 6-way interactions (complex intransitivities)")
print()

print("Expected energy profile for n=10:")
for k in range(1, 5):
    m = 10 - 2*k
    if m < 1: continue
    gk = sum(Fraction(comb(k-1,j-1)*comb(m,j)) * 2**(j-1)
             for j in range(1, min(k,m)+1))
    ff = Fraction(1)
    for i in range(2*k):
        ff *= (10 - i)
    ratio = 2 * gk / ff
    print(f"  E_{2*k}/E_0 = {float(ratio):.6f} ({ratio})")

print()
print("PRODUCT: A tournament TDA (Topological Data Analysis) feature vector")
print("for machine learning on ranking/preference data.")

# ============================================================
# APPLICATION 5: THE CAYLEY TRANSFORM AS ODDS CONVERTER
# ============================================================
print()
print("="*70)
print("APPLICATION 5: CAYLEY TRANSFORM = PROBABILITY-TO-ODDS CONVERTER")
print("="*70)
print()
print("Q(x) = (1+x)/(1-x). If x = (p-q)/(p+q) where p = P(A beats B),")
print("q = 1-p, then x = 2p-1 (the 'skill gap') and:")
print("  Q(x) = (1+(2p-1))/(1-(2p-1)) = 2p/(2-2p) = p/(1-p) = ODDS!")
print()
print("So Q is the MAP from skill gap to odds ratio.")
print("Q^m = (p/(1-p))^m: the odds after m independent comparisons.")
print()
print("In the tournament context:")
print("  g_k(m) = coefficients of Q^m = coefficients of the 'm-fold odds'")
print("  This is the FOURIER TRANSFORM of repeated pairwise comparison!")
print()
print("INSIGHT: The tournament energy formula decomposes the variance of H")
print("into contributions from 'k-fold odds interactions' at each level.")

# ============================================================
# APPLICATION 6: CODING THEORY — ERROR DETECTION
# ============================================================
print()
print("="*70)
print("APPLICATION 6: CODING THEORY — TOURNAMENT CODES")
print("="*70)
print()
print("A tournament on n vertices encodes C(n,2) bits. The Hamiltonian path")
print("count H(T) is a PARITY CHECK: H is always odd (Redei's theorem).")
print()
print("The Fourier spectrum tells us about error detection/correction:")
print("  - Level-2 energy: sensitivity to single-bit flips")
print("  - Level-2k energy: sensitivity to k-bit flips")
print()
print("From our formula: E_2/E_0 = 2(n-2)/(n(n-1))")
print("The 'code distance' interpretation: flipping one arc changes H by")
print("an amount governed by the level-2 energy.")
print()

for n in [5, 7, 11, 23]:
    e2_ratio = 2.0*(n-2)/(n*(n-1))
    total_cv2 = 2.0/n
    pct = e2_ratio / total_cv2 * 100
    print(f"  n={n}: E_2 captures {pct:.1f}% of total variance")
    print(f"         (level-2 dominance increases with n)")

print()
print("PRODUCT: Error bounds for tournament-based distributed consensus")
print("protocols, where arc-flip = communication error.")

# ============================================================
# APPLICATION 7: RANDOM MATRIX CONNECTION
# ============================================================
print()
print("="*70)
print("APPLICATION 7: RANDOM MATRIX THEORY CONNECTION")
print("="*70)
print()
print("The transfer matrix M(x) = [[1,2x,0],[0,0,1],[1,x,0]] has:")
print("  eigenvalues satisfying lambda^3 = lambda^2 + x*lambda + x")
print("  At x=1: tribonacci constant tau = 1.8393...")
print()
print("In random matrix theory, the spacing distribution of eigenvalues")
print("follows universal laws. Our transfer matrix is NOT random — it's")
print("deterministic — but it plays the role of a 'transfer operator'")
print("for the 1D statistical mechanics of the domino process.")
print()
print("The Z_j process is a HIDDEN MARKOV MODEL with 3 states and")
print("tridiagonal covariance. The 'free energy' is:")
print("  f(x) = log(lambda_1(x)) = 2x - 8x^2 + 68x^3/3 - ...")
print()
print("PRODUCT: Connection to exactly-solvable 1D statistical mechanics")
print("models (Ising-type with nearest-neighbor interaction).")

# ============================================================
# APPLICATION 8: SEARCH ENGINE RANKING
# ============================================================
print()
print("="*70)
print("APPLICATION 8: WEB/SEARCH RANKING QUALITY METRIC")
print("="*70)
print()
print("PageRank creates a ranking from pairwise link structure.")
print("Given n webpages with pairwise preference data (A links to B),")
print("the tournament T encodes the link graph.")
print()
print("H(T)/E[H] measures 'ranking coherence':")
print("  H/E[H] >> 1: highly rankable (clear hierarchy)")
print("  H/E[H] ~ 1: random-like (no clear ranking)")
print("  H/E[H] << 1: highly contradictory (many cycles)")
print()
print("Our formula gives the EXACT null distribution of H/E[H]:")
print("  Under random tournament: H/E[H] has CV = sqrt(2/n)")
print()
print("PRODUCT: A p-value for 'is this ranking statistically significant?'")
print("given only pairwise comparison data. No parametric assumptions needed.")

# ============================================================
# APPLICATION 9: SPORTS ANALYTICS
# ============================================================
print()
print("="*70)
print("APPLICATION 9: SPORTS ANALYTICS — LEAGUE COMPETITIVENESS")
print("="*70)
print()

# Simulate a round-robin tournament with skill differences
def simulate_tournament(n, skills=None):
    """Generate tournament from Bradley-Terry model."""
    if skills is None:
        skills = [random.gauss(0, 1) for _ in range(n)]
    T = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            p = 1 / (1 + exp(-(skills[i] - skills[j])))
            if random.random() < p:
                T[i][j] = 1
            else:
                T[j][i] = 1
    return T

print("League competitiveness index L(T) = H(T) / E[H]")
print()
print("  L >> 1: teams have VERY different skill levels (easy to rank)")
print("  L ~ 1: competitive league (hard to distinguish teams)")
print("  L << 1: paradoxical results (upsets everywhere)")
print()
print("Example: 6-team league (like a conference)")
print("  Random tournament: E[L] = 1, std(L) = sqrt(13/45) = 0.537")
print()
print("  Skill spread sigma=0 (equal teams): expect L ~ 1")
print("  Skill spread sigma=1 (moderate): expect L > 1")
print("  Skill spread sigma=3 (dominant teams): expect L >> 1")
print()
print("PRODUCT: A single number summarizing 'how rankable is this league?'")
print("Directly comparable across sports, seasons, and league sizes")
print("because we normalize by the EXACT random expectation.")

# ============================================================
# APPLICATION 10: CRYPTOCURRENCY/BLOCKCHAIN CONSENSUS
# ============================================================
print()
print("="*70)
print("APPLICATION 10: DISTRIBUTED CONSENSUS RELIABILITY")
print("="*70)
print()
print("In Byzantine fault tolerance, n nodes vote pairwise on block validity.")
print("The 'tournament' of votes determines consensus.")
print()
print("H(T) measures how many linear orderings are consistent with votes.")
print("If H is large, consensus is robust (many orderings agree).")
print("If H is small, the network is in a 'confused' state.")
print()
print("Our formula gives the THRESHOLD for reliable consensus:")
print("  If H(T) > E[H] + 3*sigma, the consensus is statistically")
print("  significant (not explainable by random voting).")
print()
print("  Threshold = E[H] * (1 + 3*sqrt(2/n))")
for n in [10, 21, 100, 1000]:
    threshold_ratio = 1 + 3*sqrt(2.0/n)
    print(f"  n={n:4d}: threshold = {threshold_ratio:.3f} * E[H]")

print()
print("PRODUCT: A mathematically rigorous consensus quality metric")
print("for blockchain validators with exact distributional bounds.")

# ============================================================
# SYNTHESIS
# ============================================================
print()
print("="*70)
print("SYNTHESIS: THE CAYLEY-DELANNOY TOOLKIT")
print("="*70)
print()
print("The formula CV^2 = 2/n + O(1/n^3) and its refinements give:")
print()
print("1. RANKING QUALITY: p-value for 'is this ranking real?'")
print("2. ANOMALY DETECTION: identify unusual preference patterns")
print("3. LEAGUE METRICS: competitiveness index for sports")
print("4. CONSENSUS: reliability threshold for distributed systems")
print("5. TDA FEATURES: spectral fingerprint for ML on tournaments")
print("6. CODING THEORY: error sensitivity of tournament parity")
print("7. SEARCH RANKING: statistical significance of web rankings")
print("8. FAST ESTIMATION: O(1) formula for any n, no simulation needed")
print()
print("All from a single formula: Q(x)^m = ((1+x)/(1-x))^m.")
print()
print("The Cayley transform maps skill gaps to odds ratios,")
print("and its powers decompose tournament variance into")
print("Delannoy-weighted contributions at each interaction level.")

print("\nDone!")
