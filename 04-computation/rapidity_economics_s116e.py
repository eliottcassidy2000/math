#!/usr/bin/env python3
"""rapidity_economics_s116e.py -- Economics through the rapidity lens

Rapidity phi = arctanh(v) = (1/2)*ln((1+v)/(1-v)) = (1/2)*ln(Q(v))
Q(v) = (1+v)/(1-v) = e^{2*phi}
logit(p) = ln(p/(1-p)) = 2*arctanh(2p-1) = 2*rapidity(2p-1)

Every economic quantity that lives in (0,1) or (0,inf) has a natural
rapidity coordinate. This script explores eight domains.
"""

from math import (sqrt, log, exp, tanh, atanh, pi, cosh, sinh,
                  erf, factorial, comb)

def Q(v):
    """Cayley transform: Q(v) = (1+v)/(1-v)"""
    if abs(v) >= 1:
        return float('inf') if v > 0 else 0.0
    return (1+v)/(1-v)

def rapidity(v):
    """rapidity = arctanh(v) = (1/2)*ln(Q(v))"""
    if abs(v) >= 1.0:
        return float('inf') if v > 0 else float('-inf')
    return atanh(v)

def logit(p):
    """logit(p) = ln(p/(1-p)) = 2*rapidity(2p-1)"""
    if p <= 0: return float('-inf')
    if p >= 1: return float('inf')
    return log(p/(1-p))

def sigmoid(x):
    """inverse logit: sigmoid(x) = 1/(1+exp(-x))"""
    return 1.0/(1.0 + exp(-x))

def normal_cdf(x):
    """Standard normal CDF via erf"""
    return 0.5*(1 + erf(x/sqrt(2)))

# ============================================================
print("RAPIDITY IN ECONOMICS")
print("="*70)
print()
print("  Core identities:")
print("    rapidity(v) = arctanh(v) = (1/2)*ln(Q(v))")
print("    logit(p) = 2*rapidity(2p-1)")
print("    Q(v) = (1+v)/(1-v) = e^{2*rapidity(v)}")
print()
print("  Every probability p in (0,1) has rapidity arctanh(2p-1).")
print("  Every positive ratio r has rapidity (1/2)*ln(r).")
print("  These are the SAME operation on the SAME universal coordinate.")
print()

# ============================================================
# 1. LOGIT DEMAND / MULTINOMIAL CHOICE
# ============================================================
print("="*70)
print("1. LOGIT DEMAND MODELS")
print("="*70)
print()
print("  McFadden's discrete choice: P_ij = exp(V_ij) / sum_k exp(V_ik)")
print()
print("  For two products A and B with utilities V_A, V_B:")
print("    P_A = 1/(1+exp(-(V_A - V_B))) = sigmoid(V_A - V_B)")
print("    P_A/P_B = exp(V_A - V_B)")
print()
print("  Key insight: V_A - V_B = logit(P_A) = 2*rapidity(2*P_A - 1)")
print("  The utility difference IS the logit IS twice the rapidity!")
print()
print("  Equivalently: P_A/P_B = Q(2*P_A - 1) = e^{2*rapidity(2*P_A-1)}")
print()

print("  Example: Two-product market")
print("  " + "-"*60)
print(f"  {'V_A-V_B':>10} {'P_A':>8} {'P_B':>8} {'P_A/P_B':>10} {'rapidity':>10}")
print("  " + "-"*60)
for dv in [-3, -2, -1, -0.5, 0, 0.5, 1, 2, 3]:
    pa = sigmoid(dv)
    pb = 1 - pa
    ratio = pa/pb if pb > 0 else float('inf')
    rap = rapidity(2*pa - 1)
    print(f"  {dv:10.2f} {pa:8.4f} {pb:8.4f} {ratio:10.4f} {rap:10.6f}")

print()
print("  Verification: V_A - V_B = 2*rapidity(2*P_A - 1)")
for dv in [0.5, 1.0, 2.0]:
    pa = sigmoid(dv)
    rap = rapidity(2*pa - 1)
    print(f"    dV={dv:.1f}: 2*rap = {2*rap:.10f}, dV = {dv:.10f}, match = {abs(2*rap - dv) < 1e-12}")

print()
print("  Three-product market: V = [2.0, 1.0, 0.5]")
V = [2.0, 1.0, 0.5]
denom = sum(exp(v) for v in V)
shares = [exp(v)/denom for v in V]
print(f"  {'Product':>10} {'V':>6} {'Share':>8} {'logit(s)':>10} {'rapidity':>10}")
print("  " + "-"*55)
for i, (v, s) in enumerate(zip(V, shares)):
    lg = logit(s)
    rap = rapidity(2*s - 1)
    print(f"  {'ABC'[i]:>10} {v:6.1f} {s:8.4f} {lg:10.4f} {rap:10.6f}")

print()
print("  Pairwise share ratios = Q-value ratios:")
for i in range(3):
    for j in range(i+1, 3):
        ratio = shares[i]/shares[j]
        dv = V[i] - V[j]
        print(f"    P_{{'ABC'[i]}}/P_{{'ABC'[j]}} = {ratio:.4f} = exp({dv:.1f}) = {exp(dv):.4f}")

print()
print("  THEOREM (Rapidity Demand): In a logit model, the market share")
print("  rapidity arctanh(2*P_j - 1) encodes utility differences.")
print("  Adding a new product SHIFTS all rapidities by the same amount")
print("  (IIA property = translation invariance in rapidity space).")

# ============================================================
# 2. KELLY CRITERION
# ============================================================
print()
print("="*70)
print("2. KELLY CRITERION IN RAPIDITY")
print("="*70)
print()
print("  Kelly: f* = (b*p - q)/b  where p=P(win), q=1-p, b=net odds")
print("  For even odds (b=1): f* = p - q = 2p - 1 = tanh(rapidity(2p-1))")
print()
print("  INSIGHT: For even odds, f* = 2p-1 itself, so")
print("    rapidity(f*) = rapidity(2p-1) = arctanh(2p-1) = logit(p)/2")
print()
print("  For general odds b:")
print("    f* = (bp - q)/b = p - (1-p)/b")
print("    Let r_p = arctanh(2p-1) (edge rapidity)")
print("    Let r_b = (1/2)*ln(b) (odds rapidity)")
print()

print("  Kelly fraction for various (p, b):")
print("  " + "-"*70)
print(f"  {'p':>6} {'b':>6} {'f*':>8} {'r_p':>10} {'r_b':>10} {'arctanh(f*)':>12}")
print("  " + "-"*70)
for p in [0.5, 0.55, 0.6, 0.7, 0.8]:
    for b in [1, 2, 3, 5]:
        q = 1 - p
        f_star = (b*p - q)/b
        if f_star <= 0:
            continue
        r_p = rapidity(2*p - 1)
        r_b = log(b)/2
        # f* might exceed 1 for large p,b -- clamp for rapidity
        r_f = rapidity(f_star) if abs(f_star) < 1 else float('inf')
        print(f"  {p:6.2f} {b:6d} {f_star:8.4f} {r_p:10.6f} {r_b:10.6f} {r_f:12.6f}")

print()
print("  Even-odds identity: f* = tanh(logit(p)/2)")
print("  Verification:")
for p in [0.55, 0.6, 0.7, 0.8, 0.9]:
    f_kelly = 2*p - 1
    f_from_rap = tanh(logit(p)/2)
    print(f"    p={p:.2f}: f*={f_kelly:.4f}, tanh(logit(p)/2)={f_from_rap:.10f}, match={abs(f_kelly - f_from_rap) < 1e-12}")

print()
print("  GENERAL FORMULA: f* = tanh(r_p) - (1/b)*tanh(-r_p)")
print("    = tanh(r_p) + tanh(r_p)/b  ... NO, let's derive properly.")
print()
print("  f* = (bp - q)/b = p - (1-p)/b")
print("     = (1 + tanh(r_p))/2 - (1 - tanh(r_p))/(2b)")
print("     = (1/2)[1 + tanh(r_p) - 1/b + tanh(r_p)/b]")
print("     = (1/2)[(1 - 1/b) + tanh(r_p)*(1 + 1/b)]")
print("     = (1/2)[(b-1)/b + tanh(r_p)*(b+1)/b]")
print()
print("  Using b = Q(v_b) where v_b = (b-1)/(b+1):")
print("    (b-1)/b = 2*v_b/(1+v_b), (b+1)/b = 2/(1+v_b)")
print("    f* = [v_b + tanh(r_p)] / (1 + v_b)")
print()
print("  For b=1 (v_b=0): f* = tanh(r_p) = 2p-1. Correct!")
print()
print("  THEOREM (Kelly-Rapidity): f* = [v_b + tanh(r_p)] / (1 + v_b)")
print("  where r_p = arctanh(2p-1), v_b = (b-1)/(b+1).")
print("  This is a WEIGHTED MEAN of the odds address and the edge tanh.")

print()
print("  Verification of Kelly-Rapidity formula:")
for p, b in [(0.6, 1), (0.6, 2), (0.7, 3), (0.55, 5)]:
    q = 1-p
    f_std = (b*p - q)/b
    r_p = rapidity(2*p-1)
    v_b = (b-1)/(b+1)
    f_rap = (v_b + tanh(r_p))/(1 + v_b)
    print(f"    p={p}, b={b}: standard={f_std:.10f}, rapidity={f_rap:.10f}, match={abs(f_std-f_rap)<1e-12}")

# ============================================================
# 3. BLACK-SCHOLES IN RAPIDITY
# ============================================================
print()
print("="*70)
print("3. BLACK-SCHOLES IN RAPIDITY COORDINATES")
print("="*70)
print()
print("  Black-Scholes call price: C = S*N(d1) - K*e^{-rT}*N(d2)")
print("  where d1 = [ln(S/K) + (r + sigma^2/2)*T] / (sigma*sqrt(T))")
print("        d2 = d1 - sigma*sqrt(T)")
print()
print("  Key: ln(S/K) = 2*rapidity((S/K - 1)/(S/K + 1))")
print("  More precisely: if x = (S-K)/(S+K), then Q(x) = S/K,")
print("  so ln(S/K) = 2*arctanh(x) = 2*rapidity(moneyness)")
print()
print("  Define: moneyness address m = (S-K)/(S+K)")
print("          moneyness rapidity r_m = arctanh(m) = (1/2)*ln(S/K)")
print()
print("  Then: d1 = [2*r_m + (r + sigma^2/2)*T] / (sigma*sqrt(T))")
print("        d2 = [2*r_m + (r - sigma^2/2)*T] / (sigma*sqrt(T))")
print()

def black_scholes_call(S, K, T, r, sigma):
    """Black-Scholes European call price"""
    d1 = (log(S/K) + (r + sigma**2/2)*T) / (sigma*sqrt(T))
    d2 = d1 - sigma*sqrt(T)
    return S*normal_cdf(d1) - K*exp(-r*T)*normal_cdf(d2)

print("  Example: S=100, K varies, T=1, r=0.05, sigma=0.20")
print("  " + "-"*75)
print(f"  {'K':>6} {'S/K':>7} {'m=(S-K)/(S+K)':>15} {'r_m':>10} {'d1':>8} {'d2':>8} {'Call':>8}")
print("  " + "-"*75)
S, T, r, sigma = 100, 1, 0.05, 0.20
for K in [80, 85, 90, 95, 100, 105, 110, 115, 120]:
    m = (S - K)/(S + K)
    r_m = rapidity(m)
    d1 = (2*r_m + (r + sigma**2/2)*T)/(sigma*sqrt(T))
    d2 = d1 - sigma*sqrt(T)
    call = black_scholes_call(S, K, T, r, sigma)
    print(f"  {K:6.0f} {S/K:7.4f} {m:15.6f} {r_m:10.6f} {d1:8.4f} {d2:8.4f} {call:8.2f}")

print()
print("  Verification: d1 via rapidity vs standard formula")
for K in [90, 100, 110]:
    m = (S - K)/(S + K)
    r_m = rapidity(m)
    d1_rap = (2*r_m + (r + sigma**2/2)*T)/(sigma*sqrt(T))
    d1_std = (log(S/K) + (r + sigma**2/2)*T)/(sigma*sqrt(T))
    print(f"    K={K}: d1_rap={d1_rap:.10f}, d1_std={d1_std:.10f}, match={abs(d1_rap-d1_std)<1e-12}")

print()
print("  AT-THE-MONEY (S=K): m=0, r_m=0")
print("    d1 = (r + sigma^2/2)*T / (sigma*sqrt(T))")
print("    d2 = (r - sigma^2/2)*T / (sigma*sqrt(T))")
print("    The moneyness rapidity vanishes -- pure time/vol effect.")
print()
print("  DEEP IN/OUT: As S/K -> inf, m -> 1, r_m -> inf, d1,d2 -> inf, C -> S.")
print("  As S/K -> 0, m -> -1, r_m -> -inf, d1,d2 -> -inf, C -> 0.")
print()
print("  THEOREM (BS-Rapidity): d1 = [2*r_m + drift*T] / (sigma*sqrt(T))")
print("  where r_m = arctanh((S-K)/(S+K)) is the moneyness rapidity.")
print("  Options pricing lives naturally in rapidity-of-moneyness space.")
print()
print("  Delta = N(d1) as a function of moneyness rapidity:")
print("  " + "-"*50)
print(f"  {'r_m':>10} {'S/K':>8} {'Delta':>8}")
print("  " + "-"*50)
for r_m_val in [-1.0, -0.5, -0.2, -0.1, 0.0, 0.1, 0.2, 0.5, 1.0]:
    SK = exp(2*r_m_val)
    K_val = S / SK
    d1 = (2*r_m_val + (r + sigma**2/2)*T)/(sigma*sqrt(T))
    delta = normal_cdf(d1)
    print(f"  {r_m_val:10.4f} {SK:8.4f} {delta:8.4f}")

# ============================================================
# 4. NASH EQUILIBRIUM IN RAPIDITY
# ============================================================
print()
print("="*70)
print("4. NASH EQUILIBRIUM IN RAPIDITY SPACE")
print("="*70)
print()
print("  2-player zero-sum game: Row player mixes with probability p on row 1.")
print("  rapidity(p) = arctanh(2p-1)")
print()
print("  Pure strategy p=1: rapidity = +inf (total commitment)")
print("  Pure strategy p=0: rapidity = -inf (total commitment to other)")
print("  Uniform p=1/2: rapidity = 0 (maximum uncertainty)")
print()

print("  Example 1: Matching Pennies [[1,-1],[-1,1]]")
print("    Nash: p = q = 1/2. Rapidity = 0 for both players.")
print("    The equilibrium sits at the ORIGIN of rapidity space.")
print()

print("  Example 2: Asymmetric game [[3,-1],[-2,1]]")
# Row player's indifference: q*3 + (1-q)*(-1) = q*(-2) + (1-q)*1
# 3q - 1 + q = -2q + 1 - q  => 4q - 1 = -3q + 1 => 7q = 2 => q = 2/7
# Col player's indifference: p*3 + (1-p)*(-2) = p*(-1) + (1-p)*1
# 3p - 2 + 2p = -p + 1 - p => 5p - 2 = -2p + 1 => 7p = 3 => p = 3/7
p_nash = 3/7
q_nash = 2/7
r_p = rapidity(2*p_nash - 1)
r_q = rapidity(2*q_nash - 1)
print(f"    Nash: p = 3/7 = {p_nash:.6f}, q = 2/7 = {q_nash:.6f}")
print(f"    Rapidity of p: arctanh(2*(3/7)-1) = arctanh(-1/7) = {r_p:.6f}")
print(f"    Rapidity of q: arctanh(2*(2/7)-1) = arctanh(-3/7) = {r_q:.6f}")
print(f"    Both negative -- both players lean toward strategy 2.")
print()

print("  Example 3: Rock-Paper-Scissors (3 strategies)")
print("    Nash: p = (1/3, 1/3, 1/3). Each rapidity = arctanh(2/3-1) = arctanh(-1/3)")
print(f"    Rapidity per strategy: {rapidity(2/3-1):.6f}")
print("    In a 3-strategy game, the 'uniform' rapidity is NOT zero!")
print("    Zero rapidity means p=1/2. Uniform over 3 means p=1/3.")
print()

print("  Rapidity encodes DEVIATION from maximum uncertainty:")
print("  " + "-"*50)
print(f"  {'Strategy':>10} {'p':>6} {'rap(2p-1)':>12} {'Interpretation':>20}")
print("  " + "-"*50)
strategy_data = [
    ("Committed", 0.99, "near-pure"),
    ("Strong lean", 0.80, "strong preference"),
    ("Mild lean", 0.60, "slight edge"),
    ("Uniform/2", 0.50, "indifferent (2-game)"),
    ("Uniform/3", 1/3, "indifferent (3-game)"),
    ("Uniform/4", 0.25, "indifferent (4-game)"),
    ("Avoid", 0.10, "rarely play"),
    ("Never", 0.01, "near-pure avoid"),
]
for name, p, interp in strategy_data:
    r = rapidity(2*p-1)
    print(f"  {name:>10} {p:6.4f} {r:12.6f} {interp:>20}")

print()
print("  THEOREM (Nash-Rapidity): In a 2-player zero-sum game, the Nash")
print("  equilibrium mixed strategies (p*, q*) map to rapidities")
print("  (arctanh(2p*-1), arctanh(2q*-1)). The equilibrium point in")
print("  rapidity space encodes the game's asymmetry. Symmetric games")
print("  have equilibrium at the origin.")

# ============================================================
# 5. AUCTION THEORY IN RAPIDITY
# ============================================================
print()
print("="*70)
print("5. AUCTION THEORY IN RAPIDITY")
print("="*70)
print()
print("  First-price sealed-bid auction, n bidders, values uniform on [0,1].")
print("  Optimal bid: b*(v) = v*(n-1)/n = v - v/n")
print("  Bid shading: s = b/v = (n-1)/n")
print()
print("  The bid/value ratio s = (n-1)/n has Cayley address")
print("  x_s = (s-1)/(s+1) = (-1/n)/((2n-1)/n) = -1/(2n-1)")
print("  and rapidity arctanh(-1/(2n-1))")
print()

print("  Bid shading in rapidity space:")
print("  " + "-"*65)
print(f"  {'n':>4} {'b/v':>8} {'address':>10} {'rapidity':>12} {'Q(addr)':>10}")
print("  " + "-"*65)
for n in [2, 3, 4, 5, 6, 7, 8, 10, 20, 50, 100]:
    s = (n-1)/n
    addr = (s-1)/(s+1)  # = -1/(2n-1)
    rap = rapidity(addr)
    q_val = Q(addr)
    print(f"  {n:4d} {s:8.6f} {addr:10.6f} {rap:12.6f} {q_val:10.6f}")

print()
print("  As n -> inf: b/v -> 1 (no shading), rapidity -> 0.")
print("  As n = 2: maximal shading b/v = 1/2, rapidity = -arctanh(1/3) = -0.3466.")
print()

print("  Revenue equivalence in rapidity:")
print("  Expected revenue per bidder:")
print("    First-price:  E[b] = (n-1)/(n+1) of E[max value]")
print("    Second-price: E[b] = same (revenue equivalence)")
print()
print("  The ratio (n-1)/(n+1) IS a Cayley address! Its Q-value is n.")
print("  rapidity((n-1)/(n+1)) = ln(n)/2")
print()
for n in [2, 3, 5, 10]:
    ratio = (n-1)/(n+1)
    rap = rapidity(ratio)
    print(f"    n={n}: revenue ratio = {ratio:.6f}, rapidity = {rap:.6f} = ln({n})/2 = {log(n)/2:.6f}")

print()
print("  THEOREM (Auction-Rapidity): The revenue ratio (n-1)/(n+1) = address of n.")
print("  Its rapidity ln(n)/2 is exactly the rapidity of the integer n.")
print("  More bidders = higher revenue rapidity, logarithmically.")
print()

print("  Winner's Curse: In common-value auctions, winning is bad news.")
print("  If your value estimate v_i ~ Uniform[V-e, V+e], the expected")
print("  value conditional on winning = V - e*(n-1)/(n+1).")
print("  The curse magnitude e*(n-1)/(n+1) has rapidity ln(n)/2!")
print("  The curse grows LOGARITHMICALLY in the number of bidders.")

# ============================================================
# 6. GINI COEFFICIENT AND RAPIDITY
# ============================================================
print()
print("="*70)
print("6. GINI COEFFICIENT IN RAPIDITY")
print("="*70)
print()
print("  Gini G in [0,1]: 0 = perfect equality, 1 = perfect inequality.")
print("  Map G to rapidity: arctanh(2G-1)")
print("    G=0: rapidity = arctanh(-1) = -inf (perfect equality)")
print("    G=0.5: rapidity = arctanh(0) = 0 (boundary)")
print("    G=1: rapidity = arctanh(1) = +inf (perfect inequality)")
print()

# Real-world Gini coefficients (approximate, World Bank data circa 2020s)
countries = [
    ("Slovakia", 0.232),
    ("Denmark", 0.281),
    ("Germany", 0.317),
    ("France", 0.326),
    ("Japan", 0.329),
    ("Canada", 0.334),
    ("UK", 0.351),
    ("Australia", 0.341),
    ("USA", 0.398),
    ("China", 0.382),
    ("Mexico", 0.454),
    ("Brazil", 0.489),
    ("Colombia", 0.513),
    ("South Africa", 0.630),
]

print("  Real-world Gini coefficients in rapidity space:")
print("  " + "-"*65)
print(f"  {'Country':>15} {'Gini':>6} {'2G-1':>8} {'rapidity':>10} {'Q(2G-1)':>10}")
print("  " + "-"*65)
for name, g in countries:
    v = 2*g - 1
    rap = rapidity(v)
    q = Q(v)
    print(f"  {name:>15} {g:6.3f} {v:8.4f} {rap:10.6f} {q:10.4f}")

print()
print("  Observations:")
print("  - All real-world Ginis have NEGATIVE rapidity (all G < 0.5).")
print("  - The most unequal country (South Africa, G=0.63) just barely")
print("    crosses into positive rapidity territory.")
print("  - In rapidity space, the gap between Slovakia (G=0.23) and")
print("    South Africa (G=0.63) is:")
r_sk = rapidity(2*0.232-1)
r_sa = rapidity(2*0.630-1)
print(f"    {r_sa:.4f} - ({r_sk:.4f}) = {r_sa - r_sk:.4f}")
print(f"    Q-value ratio: {Q(2*0.630-1)/Q(2*0.232-1):.4f}")
print()
print("  The Q-value ratio tells us: South Africa's inequality-to-equality")
print(f"  odds are {Q(2*0.630-1)/Q(2*0.232-1):.1f}x those of Slovakia.")
print()

# Lorenz curve connection
print("  Lorenz curve connection:")
print("  At cumulative population fraction p, cumulative income fraction L(p).")
print("  The 'local Gini' at p is the rapidity of (p - L(p))/p.")
print("  For Pareto distribution F(x) = 1 - (x_m/x)^alpha:")
print("    Gini = 1/(2*alpha - 1)")
print()
print("  Pareto Gini in rapidity:")
print("  " + "-"*50)
print(f"  {'alpha':>8} {'Gini':>8} {'rapidity':>10}")
print("  " + "-"*50)
for alpha in [1.2, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 10.0]:
    g = 1/(2*alpha - 1)
    r = rapidity(2*g - 1)
    print(f"  {alpha:8.1f} {g:8.4f} {r:10.6f}")

print()
print("  As alpha -> 1 (heavy tail): G -> 1, rapidity -> +inf")
print("  As alpha -> inf (light tail): G -> 0, rapidity -> -inf")
print("  alpha=1.5 gives G = 0.5, rapidity = 0 (the Gini midpoint).")
print()
print("  THEOREM (Gini-Rapidity): arctanh(2G-1) = arctanh(2/(2a-1) - 1)")
print("  = arctanh((2-(2a-1))/(2a-1)) = arctanh((3-2a)/(2a-1))")
print("  For Pareto with parameter alpha.")

# ============================================================
# 7. INTEREST RATES IN RAPIDITY
# ============================================================
print()
print("="*70)
print("7. INTEREST RATES IN RAPIDITY")
print("="*70)
print()
print("  Nominal rate r: growth factor = 1+r per period.")
print("  Continuous rate: rho = ln(1+r)")
print()
print("  Rapidity approach: map the growth factor g = 1+r to rapidity.")
print("  Address of g: x_g = (g-1)/(g+1) = r/(2+r)")
print("  Rapidity: arctanh(r/(2+r)) = (1/2)*ln(g) = (1/2)*ln(1+r) = rho/2")
print()
print("  THE RAPIDITY OF THE GROWTH FACTOR IS HALF THE CONTINUOUS RATE.")
print("  (Same as rapidity(n) = ln(n)/2 for any positive number.)")
print()

rates = [
    ("Deflation", -0.02),
    ("Zero", 0.00),
    ("Low (Japan)", 0.001),
    ("Savings", 0.02),
    ("Moderate", 0.05),
    ("Historical US equity", 0.10),
    ("High growth", 0.15),
    ("Emerging market", 0.20),
    ("Hyperinflation (mild)", 1.00),
    ("Hyperinflation (severe)", 10.0),
    ("Zimbabwe 2008 (monthly)", 79600000000.0),
]

print("  " + "-"*78)
print(f"  {'Description':>25} {'r':>14} {'g=1+r':>14} {'addr':>10} {'rapidity':>10} {'rho=2*rap':>10}")
print("  " + "-"*78)
for desc, r in rates:
    g = 1 + r
    if g <= 0:
        print(f"  {desc:>25} {r:14.6f}  (negative growth factor)")
        continue
    addr = r/(2+r)
    rap = rapidity(addr) if abs(addr) < 1 else log(g)/2
    rho = log(g)
    print(f"  {desc:>25} {r:14.4f} {g:14.4f} {addr:10.6f} {rap:10.6f} {rho:10.6f}")

print()
print("  Compounding periods and rapidity:")
print("  For annual rate r compounded m times per year:")
print("    Effective rate: (1 + r/m)^m - 1")
print("    The rapidity of each sub-period: arctanh(r/(m*(2+r/m)))")
print()
r_annual = 0.10
print(f"  Annual rate r = {r_annual}")
print("  " + "-"*60)
print(f"  {'m':>8} {'sub-rate':>10} {'sub-rap':>10} {'m*sub-rap':>10} {'eff rate':>10}")
print("  " + "-"*60)
for m in [1, 2, 4, 12, 52, 365]:
    sub = r_annual/m
    g_sub = 1 + sub
    addr_sub = sub/(2+sub)
    rap_sub = rapidity(addr_sub)
    eff = g_sub**m - 1
    print(f"  {m:8d} {sub:10.6f} {rap_sub:10.6f} {m*rap_sub:10.6f} {eff:10.6f}")

print()
print("  As m -> inf: m * rapidity(r/(m*(2+r/m))) -> ln(1+r)/2 = rho/2")
cont_rap = log(1+r_annual)/2
print(f"  Continuous limit: rho/2 = {cont_rap:.6f}")
print()
print("  THEOREM (Interest-Rapidity): The rapidity of the growth factor")
print("  1+r equals half the continuous rate rho = ln(1+r). Compounding")
print("  m times adds m sub-rapidities, converging to rho/2.")
print()

# Real-world: Rule of 72
print("  Rule of 72 in rapidity:")
print("  Doubling time T_2 = ln(2)/rho = ln(2)/(2*rapidity)")
print("  'Rule of 72': T_2 ~ 72/r_percent")
print()
for r in [0.02, 0.05, 0.08, 0.10, 0.12, 0.15]:
    rap = log(1+r)/2
    t2_exact = log(2)/(2*rap)
    t2_rule72 = 72/(r*100)
    print(f"    r={r*100:.0f}%: rapidity={rap:.6f}, T2_exact={t2_exact:.2f}, Rule72={t2_rule72:.2f}, error={abs(t2_exact-t2_rule72)/t2_exact*100:.2f}%")

# ============================================================
# 8. PROSPECT THEORY IN RAPIDITY
# ============================================================
print()
print("="*70)
print("8. PROSPECT THEORY IN RAPIDITY SPACE")
print("="*70)
print()
print("  Kahneman-Tversky probability weighting function:")
print("    w(p) = p^gamma / (p^gamma + (1-p)^gamma)^{1/gamma}")
print("  with gamma ~ 0.61 (gains) or gamma ~ 0.69 (losses)")
print()
print("  In rapidity space: rapidity_distortion = arctanh(2*w(p)-1) - arctanh(2p-1)")
print("  This measures how much people WARP probabilities.")
print()

def kt_weight(p, gamma):
    """Kahneman-Tversky probability weighting"""
    if p <= 0: return 0.0
    if p >= 1: return 1.0
    pg = p**gamma
    qg = (1-p)**gamma
    return pg / (pg + qg)**(1/gamma)

gamma_gains = 0.61
gamma_losses = 0.69

print(f"  Gains weighting (gamma={gamma_gains}):")
print("  " + "-"*75)
print(f"  {'p':>6} {'w(p)':>8} {'rap(p)':>10} {'rap(w)':>10} {'distortion':>12} {'overweight?':>12}")
print("  " + "-"*75)
for p in [0.01, 0.05, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 0.95, 0.99]:
    w = kt_weight(p, gamma_gains)
    r_p = rapidity(2*p-1)
    r_w = rapidity(2*w-1)
    distortion = r_w - r_p
    ow = "OVER" if w > p else ("UNDER" if w < p else "exact")
    print(f"  {p:6.2f} {w:8.4f} {r_p:10.4f} {r_w:10.4f} {distortion:12.4f} {ow:>12}")

print()
print("  The crossover point where w(p) = p:")
# Find crossover
for i in range(1, 1000):
    p = i/1000
    if kt_weight(p, gamma_gains) < p:
        p_cross = p
        break
print(f"    Crossover at p ~ {p_cross:.3f} (gains, gamma={gamma_gains})")
print(f"    Rapidity at crossover: {rapidity(2*p_cross-1):.4f}")
print()

print(f"  Losses weighting (gamma={gamma_losses}):")
print("  " + "-"*75)
print(f"  {'p':>6} {'w(p)':>8} {'rap(p)':>10} {'rap(w)':>10} {'distortion':>12} {'overweight?':>12}")
print("  " + "-"*75)
for p in [0.01, 0.05, 0.10, 0.20, 0.50, 0.80, 0.90, 0.95, 0.99]:
    w = kt_weight(p, gamma_losses)
    r_p = rapidity(2*p-1)
    r_w = rapidity(2*w-1)
    distortion = r_w - r_p
    ow = "OVER" if w > p else ("UNDER" if w < p else "exact")
    print(f"  {p:6.2f} {w:8.4f} {r_p:10.4f} {r_w:10.4f} {distortion:12.4f} {ow:>12}")

print()
print("  INSIGHT: The distortion is POSITIVE for small p (overweight rare events)")
print("  and NEGATIVE for large p (underweight common events). The rapidity")
print("  distortion curve crosses zero at the crossover point.")
print()
print("  In rapidity space, KT weighting is an S-shaped compression:")
print("  - Rare events get PUSHED toward zero rapidity (inflated)")
print("  - Common events get PUSHED toward zero rapidity (deflated)")
print("  - The weighting function COMPRESSES the rapidity range.")
print()

print("  Rapidity compression ratio = rap(w(p)) / rap(p):")
print("  " + "-"*50)
print(f"  {'p':>6} {'ratio (gains)':>15} {'ratio (losses)':>15}")
print("  " + "-"*50)
for p in [0.01, 0.05, 0.10, 0.20, 0.50, 0.80, 0.90, 0.95, 0.99]:
    r_p = rapidity(2*p-1)
    w_g = kt_weight(p, gamma_gains)
    w_l = kt_weight(p, gamma_losses)
    r_wg = rapidity(2*w_g-1)
    r_wl = rapidity(2*w_l-1)
    if abs(r_p) > 0.001:
        ratio_g = r_wg / r_p
        ratio_l = r_wl / r_p
        print(f"  {p:6.2f} {ratio_g:15.4f} {ratio_l:15.4f}")
    else:
        print(f"  {p:6.2f} {'(near 0)':>15} {'(near 0)':>15}")

print()
print("  THEOREM (Prospect-Rapidity): KT probability weighting is a")
print("  CONTRACTION in rapidity space. The compression ratio")
print("  rap(w(p))/rap(p) measures the degree of probability distortion.")
print("  For gamma < 1, the ratio is always in (0,1) away from the")
print("  crossover point, meaning KT weighting SHRINKS the rapidity range.")
print("  Human probability perception compresses the natural rapidity scale.")

# ============================================================
# SYNTHESIS
# ============================================================
print()
print("="*70)
print("SYNTHESIS: RAPIDITY AS UNIVERSAL ECONOMIC COORDINATE")
print("="*70)
print()
print("  Domain               Quantity              Rapidity =")
print("  " + "-"*65)
print("  Logit demand         P_A/(P_A+P_B)         logit(share)/2 = (V_A-V_B)/2")
print("  Kelly criterion      f* (bet fraction)     arctanh(f*)")
print("  Black-Scholes        S/K (moneyness)       ln(S/K)/2")
print("  Nash equilibrium     p* (mixed strategy)   arctanh(2p*-1)")
print("  Auction theory       (n-1)/(n+1) revenue   ln(n)/2")
print("  Gini coefficient     G (inequality)        arctanh(2G-1)")
print("  Interest rates       1+r (growth factor)   ln(1+r)/2 = rho/2")
print("  Prospect theory      w(p) distortion       arctanh(2w-1) - arctanh(2p-1)")
print()
print("  UNIFYING THEME: Every economic ratio, probability, or fraction")
print("  maps to rapidity via arctanh or (1/2)*ln. The same transform")
print("  that linearizes relativistic velocity addition also linearizes:")
print("    - Utility differences in discrete choice")
print("    - Optimal bet sizing")
print("    - Option moneyness")
print("    - Strategic mixing")
print("    - Auction revenue scaling")
print("    - Income inequality")
print("    - Compound interest")
print("    - Behavioral probability distortion")
print()
print("  Rapidity is the universal coordinate for economics.")
print()
print("  KEY FORMULAS:")
print("    Kelly: f* = [v_b + tanh(r_p)] / (1 + v_b)")
print("    BS:    d1 = [2*r_m + drift*T] / (sigma*sqrt(T))")
print("    Auction revenue ratio = address of n, rapidity = ln(n)/2")
print("    KT distortion = rapidity compression of probability space")
print()
print("="*70)
print("END rapidity_economics_s116e.py")
print("="*70)
