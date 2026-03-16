#!/usr/bin/env python3
"""logit_rapidity_s116d.py — The logit function IS twice the rapidity.

This is perhaps the most PRACTICAL connection in the entire project.
Every logistic regression, every neural network sigmoid, every Elo rating,
every Bradley-Terry model is computing in rapidity space.

logit(p) = ln(p/(1-p)) = 2 * arctanh(2p-1) = 2 * rapidity(2p-1)
"""
from math import sqrt, log, exp, atanh, tanh

print("LOGIT = 2 * RAPIDITY: THE MOST PRACTICAL CONNECTION")
print("="*70)
print()

# ============================================================
print("="*70)
print("THE IDENTITY")
print("="*70)
print()
print("  logit(p) = ln(p/(1-p))")
print("  rapidity(v) = arctanh(v) = (1/2)*ln((1+v)/(1-v))")
print()
print("  Set v = 2p - 1 (map [0,1] to [-1,1]):")
print("  rapidity(2p-1) = (1/2)*ln((1+(2p-1))/(1-(2p-1)))")
print("                 = (1/2)*ln(2p/(2-2p))")
print("                 = (1/2)*ln(p/(1-p))")
print("                 = logit(p)/2")
print()
print("  Therefore: logit(p) = 2 * rapidity(2p-1). QED.")
print()
print("  And inversely: sigmoid(x) = (1 + tanh(x/2))/2.")
print("  The SIGMOID is the shifted, scaled TANH.")
print("  (This is well-known but rarely stated in rapidity language.)")
print()

# Verify
print("  Verification:")
for p in [0.01, 0.1, 0.25, 0.5, 0.75, 0.9, 0.99]:
    logit = log(p/(1-p))
    rap = 2 * atanh(2*p - 1)
    print(f"  p={p:.2f}: logit={logit:+8.4f}, 2*rapidity={rap:+8.4f}, match={abs(logit-rap)<1e-10}")
print()

# ============================================================
print("="*70)
print("CONSEQUENCE 1: ELO RATINGS ARE RAPIDITY")
print("="*70)
print()
print("  Elo: P(A wins) = 1/(1 + 10^{(R_B - R_A)/400})")
print("  = sigmoid((R_A - R_B) * ln(10) / 400)")
print()
print("  Elo difference R_A - R_B = 400/ln(10) * logit(P(A wins))")
print("  = 400/ln(10) * 2 * rapidity(2P-1)")
print("  = 800/ln(10) * rapidity(2P-1)")
print(f"  = {800/log(10):.2f} * rapidity(2P-1)")
print()
print("  So 1 rapidity unit = {:.1f} Elo points.".format(800/log(10)))
print()
print("  Key Elo differences and their rapidities:")
print("  Elo diff   P(win)    rapidity   musical interval?")
print("  " + "-"*55)

elo_scale = 800/log(10)
for elo_diff in [0, 100, 200, 400, 800]:
    p = 1 / (1 + 10**(-elo_diff/400))
    rap = elo_diff / elo_scale
    # Check musical
    q_val = exp(2*rap)
    interval = ""
    if abs(q_val - 1) < 0.01: interval = "unison"
    elif abs(q_val - 2) < 0.1: interval = "~ octave"
    elif abs(q_val - 1.5) < 0.05: interval = "~ fifth"
    elif abs(q_val - 4/3) < 0.05: interval = "~ fourth"
    elif abs(q_val - 5/4) < 0.05: interval = "~ maj 3rd"
    print(f"  {elo_diff:5d}     {p:.4f}     {rap:.4f}     {interval}")

print()
print(f"  400 Elo points = rapidity {400/elo_scale:.4f}")
print(f"  = Q-value {exp(2*400/elo_scale):.4f}")
print(f"  = the odds ratio p/(1-p) = 10:1")
print()
print(f"  ONE OCTAVE of rapidity = {elo_scale * log(2)/2:.1f} Elo points")
print(f"  This corresponds to P(win) = {1/(1+exp(-log(2))):.4f}")
print()

# ============================================================
print("="*70)
print("CONSEQUENCE 2: NEURAL NETWORK WEIGHTS ARE RAPIDITIES")
print("="*70)
print()
print("  A single neuron: y = sigmoid(w*x + b)")
print("  = (1 + tanh((w*x + b)/2))/2")
print()
print("  The pre-activation w*x + b = logit(y) = 2*rapidity(2y-1).")
print("  So neural network weights operate in RAPIDITY SPACE.")
print()
print("  The weight w is a RAPIDITY PER INPUT UNIT:")
print("  changing x by 1 changes the rapidity by w/2.")
print()
print("  The bias b shifts the rapidity origin:")
print("  b = 0 => y = 0.5 (fair coin, rapidity = 0)")
print("  b > 0 => y > 0.5 (biased toward 1)")
print("  b < 0 => y < 0.5 (biased toward 0)")
print()
print("  PRACTICAL INSIGHT:")
print("  Weight initialization with std = 1/sqrt(fan_in) ensures")
print("  that initial rapidities are O(1), preventing saturation.")
print("  This is standard practice (Glorot init) but the RAPIDITY")
print("  interpretation explains WHY it works:")
print("  you want initial rapidities near 0 (fair coin),")
print("  not near infinity (saturated).")
print()

# ============================================================
print("="*70)
print("CONSEQUENCE 3: BRADLEY-TERRY MODEL IS RAPIDITY COMPARISON")
print("="*70)
print()
print("  Bradley-Terry: P(i > j) = exp(s_i) / (exp(s_i) + exp(s_j))")
print("  = sigmoid(s_i - s_j)")
print()
print("  The skill parameter s_i is a RAPIDITY.")
print("  The win probability depends on RAPIDITY DIFFERENCE s_i - s_j.")
print("  Rapidity differences ADD: if A is 1 rapidity ahead of B,")
print("  and B is 1 rapidity ahead of C, then A is 2 rapidity ahead of C.")
print()
print("  This means TOURNAMENTS BUILT FROM BRADLEY-TERRY are")
print("  tournaments where arc probabilities come from rapidity differences.")
print("  The expected H(T) depends on the rapidity spread of the skills.")
print()
print("  Small spread (all s_i similar): near-random tournament, H ~ E[H]")
print("  Large spread (s_i far apart): near-transitive, H ~ 1")
print()
print("  THE RAPIDITY SPREAD IS THE TOURNAMENT'S 'TEMPERATURE':")
print("  high temperature = random = many paths = high H")
print("  low temperature = ordered = few paths = low H")
print()
print("  This connects DIRECTLY to the Boltzmann distribution on tournaments")
print("  from the thermodynamics exploration!")
print()

# ============================================================
print("="*70)
print("CONSEQUENCE 4: INFORMATION GEOMETRY")
print("="*70)
print()
print("  The Fisher information metric on Bernoulli(p):")
print("  ds^2 = dp^2 / (p*(1-p)) = 4 * d(rapidity)^2")
print()
print("  The factor of 4 is because logit = 2*rapidity.")
print("  If we parametrize by rapidity phi = arctanh(2p-1),")
print("  then ds^2 = 4*d(phi)^2.")
print()
print("  The DISTANCE between two Bernoulli distributions p and q:")
print("  d_FR(p, q) = 2 * |arctanh(2p-1) - arctanh(2q-1)|")
print("  = |logit(p) - logit(q)|")
print("  = |rapidity(p) - rapidity(q)| * 2")
print()
print("  EVERY statistical divergence between Bernoulli distributions")
print("  is a RAPIDITY DISTANCE (up to factor 2).")
print()
print("  This means: KL divergence, Jensen-Shannon, Hellinger distance")
print("  all simplify in rapidity coordinates.")
print()

# ============================================================
print("="*70)
print("CONSEQUENCE 5: LOG-ODDS IN MEDICAL STATISTICS")
print("="*70)
print()
print("  In epidemiology:")
print("  - Odds ratio OR = (p1/(1-p1)) / (p0/(1-p0))")
print("  - ln(OR) = logit(p1) - logit(p0) = 2*(rapidity1 - rapidity0)")
print("  - The log-odds ratio IS twice the rapidity difference.")
print()
print("  In clinical trials:")
print("  - Number needed to treat NNT = 1/(p1-p0)")
print("  - In rapidity: NNT = 1/(tanh(rap1+...)/2 - tanh(rap0+...)/2)")
print()
print("  The RAPIDITY of a disease risk:")
risk_examples = [
    ("Baseline cancer", 0.001),
    ("Smoker cancer", 0.01),
    ("Coin flip", 0.5),
    ("Common cold/year", 0.3),
    ("Flu vaccine efficacy", 0.6),
]
for name, p in risk_examples:
    rap = atanh(2*p - 1)
    logit_val = log(p/(1-p))
    print(f"  {name:25s}: p={p:.3f}, rapidity={rap:+.4f}, logit={logit_val:+.4f}")

print()

# ============================================================
print("="*70)
print("THE GRAND PRACTICAL INSIGHT")
print("="*70)
print()
print("  Rapidity = arctanh = logit/2 = half-log-odds.")
print()
print("  Every field that uses log-odds is using rapidity:")
print("  - Machine learning (logistic regression, neural networks)")
print("  - Sports analytics (Elo, Glicko, TrueSkill)")
print("  - Medical statistics (odds ratios, risk analysis)")
print("  - Information theory (Fisher-Rao metric)")
print("  - Voting theory (Bradley-Terry model)")
print("  - Economics (logit demand models)")
print("  - Psychology (item response theory, Rasch models)")
print("  - Signal processing (bilinear transform)")
print()
print("  Our Cayley-Delannoy framework gives all these fields:")
print("  1. A UNIFIED LANGUAGE (rapidity instead of ad-hoc log-odds)")
print("  2. COMPOSITION RULES (rapidity adds under relativistic velocity)")
print("  3. A MUSICAL INTERPRETATION (consonance = small rapidity)")
print("  4. GEOMETRIC MEANING (hyperbolic distance on the Poincare line)")
print("  5. TOURNAMENT STRUCTURE (H(T) counts consistent orderings)")
print()
print("  The value proposition: you're already using rapidity everywhere.")
print("  We just showed you what it IS and what else it can do.")
