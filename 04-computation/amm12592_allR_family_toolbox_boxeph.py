"""AMM 12592 angle-4 toolbox: ballot coordinates, closed-form C, witness anatomy.

Exact integer arithmetic only. Polynomials = list of ints, index = power of x.
Bernstein cells at degree d: Delta = sum_k delta[k] * x^(d-k) * (1-x)^k.
"""
import json, os
from math import comb

HERE = os.path.dirname(os.path.abspath(__file__))

# ---------- exact poly ops ----------
def padd(a, b):
    n = max(len(a), len(b)); r = [0]*n
    for i,v in enumerate(a): r[i] += v
    for i,v in enumerate(b): r[i] += v
    while r and r[-1] == 0: r.pop()
    return r

def pneg(a): return [-v for v in a]
def psub(a, b): return padd(a, pneg(b))
def pmul(a, b):
    if not a or not b: return []
    r = [0]*(len(a)+len(b)-1)
    for i,u in enumerate(a):
        if u:
            for j,v in enumerate(b): r[i+j] += u*v
    while r and r[-1] == 0: r.pop()
    return r
def pshift(a, s): return [0]*s + a[:] if a else []
def pscale(a, c): return [c*v for v in a] if c else []
def qpow(m):  # (1-x)^m
    return [((-1)**k)*comb(m,k) for k in range(m+1)]

def bern_to_poly(delta, d):
    P = []
    for k, c in enumerate(delta):
        if c: P = padd(P, pscale(pshift(qpow(k), d-k), c))
    return P

def poly_to_bern(P, d):
    """cells delta_t = sum_j a_j binom(d-j, t) (P = sum a_j x^j, deg <= d)."""
    assert len(P) <= d+1, (len(P), d)
    return [sum(P[j]*comb(d-j, t) for j in range(len(P))) for t in range(d+1)]

def ballot(d):
    return [comb(d-1,k) - (comb(d-1,k-1) if k>=1 else 0) for k in range(d+1)]

def admissible(delta, d):
    if len(delta) != d+1: return False
    return all(abs(c) <= comb(d,k) and (c-comb(d,k)) % 2 == 0 for k,c in enumerate(delta))

def load_witnesses():
    W = {}
    d3 = json.load(open(os.path.join(HERE, "amm12592_floor_witnesses_R8_R16_R32.json")))
    for w in d3: W[w["R"]] = w
    for R, fn in ((64, "amm12592_floor_witness_R64_cellslim.json"),
                  (128, "amm12592_floor_witness_R128_direct.json")):
        W[R] = json.load(open(os.path.join(HERE, fn)))
    return W

def epoch_sum(R, blocks, prof):
    S = []
    for i in range(R):
        S = padd(S, pshift(bern_to_poly(blocks[i], prof[i]), i))
    return S

def residuals(R, blocks, prof):
    """sigma_{-1}=q^{R-1}; x*sigma_i = sigma_{i-1} - Delta_i. Returns [sigma_0..sigma_{R-1}]."""
    sig = qpow(R-1); out = []
    for i in range(R):
        t = psub(sig, bern_to_poly(blocks[i], prof[i]))
        assert (not t) or t[0] == 0, f"row {i}: residual not divisible by x"
        sig = t[1:] if t else []
        out.append(sig)
    return out

def closed_C(R):
    """C = x + q^{R-1} - sum_{j=2}^{R-1} x^j  (solves qC = p^R+q^R-p(p-q))."""
    C = padd([0,1], qpow(R-1))
    for j in range(2, R): C = psub(C, pshift([1], j))
    return C

if __name__ == "__main__":
    # sanity: closed_C satisfies q*C = p^R + q^R - p(p-q) for many R
    for R in [3,4,5,6,7,8,12,16,20,32,48,64,100,128,200,256]:
        lhs = pmul([1,-1], closed_C(R))
        rhs = psub(padd(pshift([1],R), qpow(R)), pmul([0,1], [-1,2]))  # p^R+q^R - p(p-q); p-q = 2x-1
        assert lhs == rhs, R
    print("closed_C identity  q*C == p^R + q^R - p(p-q)  verified for R in {3..256} sample: PASS")
