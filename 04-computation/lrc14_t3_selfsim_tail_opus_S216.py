"""
opus-2026-07-11-S216: the t>=3 self-similarity + tail-contraction investigation.

Tests death-star-S11's exact self-similarity  A(psi) = (1-2beta) conj(chat(psi))  and whether the
(5/7)-per-level contraction bounds the layer tail so LM/Q = (b/Q)^13 + sum_t layer_t > 0 a priori
for GENERIC (dissociated) families, with the relation-lattice part handled exactly.

Framework (klein THM-684, corrected S233):
  band B = { r in [1,q-1] : ceil(q/14) <= r <= floor(13q/14) },  b = |B|,  Q = q-1,  beta = b/Q
  A_t(U) = #{ c in (Z/q)^x : c*u mod q in B for every u in U }   (common-multiplier / partial-live count)
  A_1 = b, A_13(full v) = LM(q).
  Pure_t(U) via Mobius/inclusion-exclusion over sub-supports; dev_t = Pure_t / Q^{t-1}.
  layer_t = (b/Q)^{13-t} * (1/Q) * sum_{|U|=t} dev_t(U);   LM/Q = (b/Q)^13 + sum_{t>=1} layer_t  (exact).
"""
import sympy
from itertools import combinations
from fractions import Fraction

def band(q):
    lo = -(-q // 14)            # ceil(q/14)
    hi = (13 * q) // 14         # floor(13q/14)
    return lo, hi, set(range(lo, hi + 1))

def A_count(q, U, Bset):
    """#{c in [1,q-1] : (c*u) % q in B for all u in U}."""
    cnt = 0
    for c in range(1, q):
        ok = True
        for u in U:
            if (c * u) % q not in Bset:
                ok = False; break
        if ok:
            cnt += 1
    return cnt

def pure_layers(q, v):
    """Return dict t -> list of (U, dev_t) for t=1..len(v), plus LM and main-term deficit.
    Pure_t(U) = sum_{V subseteq U} (-b/Q)^{|U\V|} * A_{|V|}(V)   [A_0 := Q]; then dev_t = Pure_t / Q^{t-1}."""
    lo, hi, Bset = band(q)
    b = len(Bset); Q = q - 1
    beta = Fraction(b, Q)
    idx = list(range(len(v)))
    # precompute A over all subsets up to size cap (t<=4 for speed; A_13 = LM separately)
    Acache = {}
    def A(U):
        key = tuple(sorted(U))
        if key not in Acache:
            Acache[key] = A_count(q, [v[i] for i in U], Bset) if U else Q
        return Acache[key]
    layers = {}
    for t in range(1, 5):
        tot = Fraction(0)
        for U in combinations(idx, t):
            # inclusion-exclusion pure part
            pure = Fraction(0)
            for r in range(t + 1):
                for V in combinations(U, r):
                    pure += (-beta) ** (t - r) * A(V)
            dev = pure / Q ** (t - 1)
            tot += dev
        layers[t] = tot
    LM = A_count(q, v, Bset)
    return b, Q, beta, layers, LM

def selfsim_check(q):
    """Verify A(psi) = (1-2beta) conj(chat(psi)) via a direct multiplicative-character sum at one psi.
    Work on (Z/q)^x with generator g; chat(psi) = (1/Q) sum_{y in B} conj(psi(y)); psi(x)=zeta^{k*ind(x)}."""
    import cmath
    lo, hi, Bset = band(q); Q = q - 1
    beta = Fraction(len(Bset), Q)
    g = sympy.primitive_root(q)
    ind = {}
    x = 1
    for e in range(Q):
        ind[x] = e
        x = (x * g) % q
    zeta = cmath.exp(2j * cmath.pi / Q)
    def chat(k):  # k-th mult character coeff of 1_B fluctuation g=1_B-beta ; ghat(k)= (1/Q) sum_B zeta^{-k ind}
        s = sum(zeta ** ((-k * ind[y]) % Q) for y in Bset) / Q
        return s  # this is chat for k!=0 (beta only hits k=0)
    # pick psi = character index 1
    k = 3 % Q
    A = sum(chat(j) * (chat((j + k) % Q)).conjugate() for j in range(Q))  # sum_chi ghat(chi)conj(ghat(chi psi))
    rhs = float(1 - 2 * beta) * (chat(k)).conjugate()
    return abs(A - rhs)

print(f"{'q':>6} {'b/Q':>8} {'1-2beta':>9} {'LM/Q':>10} {'(b/Q)^13':>10} {'L1':>9} {'L2':>9} {'L3':>10} {'L4':>10} {'selfsim_err':>12}")
# GEN: spread, coprime, few relations; DIL: AP {1..13}, maximal relations
fams = {
  "GEN": [1, 5, 11, 17, 23, 28, 36, 41, 49, 55, 61, 67, 73],
  "DIL": [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13],
}
for name, v in fams.items():
    print(f"--- {name}  v={v}")
    for q in [149, 251, 401, 601, 1009, 1499, 2003]:
        if not sympy.isprime(q): continue
        b, Q, beta, layers, LM = pure_layers(q, v)
        main = beta ** 13
        LMQ = Fraction(LM, Q)
        err = selfsim_check(q) if q <= 601 else float('nan')
        print(f"{q:>6} {float(beta):>8.4f} {float(1-2*beta):>9.4f} {float(LMQ):>10.5f} {float(main):>10.5f} "
              f"{float(layers[1]):>9.4f} {float(layers[2]):>9.4f} {float(layers[3]):>10.5f} {float(layers[4]):>10.6f} {err:>12.2e}")
    # tail-contraction test: is |L_{t+1}| / |L_t| ~ 5/7?
print("\n=== tail contraction ratios |L_{t+1}/L_t| (target ~ 5/7 = 0.714 if geometric) ===")
for name, v in fams.items():
    for q in [401, 1009, 2003]:
        if not sympy.isprime(q): continue
        b, Q, beta, layers, LM = pure_layers(q, v)
        r23 = abs(float(layers[3] / layers[2])) if layers[2] != 0 else float('nan')
        r34 = abs(float(layers[4] / layers[3])) if layers[3] != 0 else float('nan')
        print(f"{name} q={q}: |L3/L2|={r23:.3f}  |L4/L3|={r34:.3f}   LM/Q-(b/Q)^13 = {float(Fraction(LM,Q)-beta**13):+.5f}")
