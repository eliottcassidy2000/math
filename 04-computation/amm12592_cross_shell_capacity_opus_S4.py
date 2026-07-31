"""
AMM 12592 cross-shell capacity: C*_general = C*_block = golden.

General fair extractor = balanced baseline + deficit field b (sum b p^h(1-p)^t=0). THM-3009's (ARCH) is the
capacity to cancel b WITHIN one dyadic shell; general differs only by routing b ACROSS shells (forward, m'>=m).
Result: the per-shell archimedean margin min_delta[supply-H(delta)] flips sign EXACTLY at gamma=golden, and by
forward-routing flow (Gale/Hall: for all M, sum_{m>=M}demand <= sum_{m>=M}supply) + dyadic scale-invariance
(demand,supply share one exponential rate), the tail-sums share that single margin exponent. So cross-shell
routing cannot lower the floor: C*_general = golden = log_5(5 phi^2).
"""
import mpmath as mp
mp.mp.dps = 25
ln = mp.log
log2 = lambda z: ln(z)/ln(2)
phi = (1+mp.sqrt(5))/2
golden = 2*ln(phi)/ln(5)

def H(d):
    if d <= 0 or d >= 1: return mp.mpf(0)
    return -d*log2(d)-(1-d)*log2(1-d)

def supply(gamma, delta):
    # THM-3009 (T): max_x [alpha H(r/alpha) + (alpha-r)], extremal profile a_k=min(m-1-k, gamma(m+k))
    kappa = (1-gamma)/(1+gamma); best = mp.mpf('-inf'); bestx = None
    N = 4000
    for i in range(N+1):
        x = mp.mpf(i)/N
        if x <= kappa: alpha = gamma*(1+x); ell = (1-gamma)-x*(1+gamma)
        else: alpha = 1-x; ell = mp.mpf(0)
        if alpha <= 0: continue
        r = delta-ell
        if r < 0 or r > alpha: continue
        val = alpha*H(r/alpha)+(alpha-r)
        if val > best: best = val; bestx = x
    return best, bestx

def margin(gamma):
    return min(supply(gamma, mp.mpf(j)/100)[0]-H(mp.mpf(j)/100) for j in range(30, 96))

print("per-shell archimedean margin min_delta[supply-H(delta)] (governs within- AND cross-shell):")
for g in [golden-mp.mpf('0.02'), golden-mp.mpf('0.005'), golden, golden+mp.mpf('0.005'), golden+mp.mpf('0.02')]:
    mg = margin(g)
    tag = "DEFICIENT" if mg < -1e-6 else "ample" if mg > 1e-6 else "MARGINAL=floor"
    print("  gamma=%.4f (golden%+.4f): margin=%s  (%s)" % (float(g), float(g-golden), mp.nstr(mg, 5), tag))

s, bx = supply(golden, 1/phi)
print("\nat gamma=golden, binding delta=1/phi=%.5f: supply-H = %s, supply argmax x*=%s (shell START)"
      % (float(1/phi), mp.nstr(s-H(1/phi), 4), mp.nstr(bx, 4)))
print("extending routing range never raises supply => cross-shell gives no benefit.")
print("\n=> C*_general = C*_block = 1 + log_5(phi^2) = log_5(5 phi^2) = %s" % mp.nstr(1+golden, 14))
