#!/usr/bin/env python3
"""
The refuel bill as a codelength (opus, 2026-10-05).

Setting (collatz_three_bit_sibling_flow_20261005.md, sections 5-6): bases B = odd b with v2(3b+1) in {1,2}; every odd
n = S^k(b) uniquely with S(u) = 4u+1 and k = floor((v2(3n+1)-1)/2); for a base b > 1 the actual successor is
U(b) = S^{k(b)}(G(b)) with G(b) in B (the base map).  The killed discounted flow K f <= rho f with f(S^j b) = r^j g(b)
is equivalent to the base inequalities
    g(b) <= d(b) g(G(b)),   d(b) = rho (1-r) r^{k(b)}           (C3),
and positive summable g exists iff every base reaches 1 (C4).  With R = -log_2 g the bill is
    R(b) - R(G(b)) >= log_2(1/(rho(1-r))) + k(b) log_2(1/r)   (= log_2(32/15) + 4k(b) for r = 1/16, rho = 1/2).

This script (all integer-exact unless stated):
A. verifies the decomposition and the base map, and the uniqueness of the child-by-depth map c -> b0(S^k c);
B. census of base edges by (a, k) with a = v2(3b+1) in {1,2}: exact fractions and the Haar prediction (2/3, 1/3) x (3/4)(1/4)^k;
C. the chain code: for every rooted base b, bill(b) = sum over its chain of [log_2(1/(rho(1-r))) + k_i log_2(1/r)] is the
   codelength of b in the base-tree code (child-by-depth sequence), so g*(b) = 2^-bill(b) is the pointwise-largest
   solution with g(1) = 1 and sum_B g* <= 1/(1-rho) (Kraft); measured: the partial Kraft sums, bill(b)/log_2 b
   (mean, quantiles, max), and the chain length T(b) = odd root time;
D. size payability: g = b^-s pays edge b iff s log_2(b/G(b)) >= bill; exact fractions of paid base edges for the
   codex parameters and for r = 1/4 with rho -> 1 (the critical bill 2 - log_2 3 + 2k), and the per-edge identity
   drop = a + 2k - log_2 3 + carry;
E. the codex kernel (88 bases): reproduce sum g < 1 and the bills of its bases.
Usage: python collatz_refuel_bill_code_20261005.py [N=1048576]
"""
import sys, math, json, os, time
from collections import Counter

def v2(n): return (n & -n).bit_length() - 1
def U(n):
    n = 3 * n + 1
    return n >> v2(n)
def S(u): return 4 * u + 1
def is_base(b): return b % 2 == 1 and v2(3 * b + 1) <= 2
def base_of(n):
    """(b, k) with n = S^k(b)"""
    k = 0
    while n % 8 == 5:
        n = (n - 1) // 4; k += 1
    return n, k
def b0(y):
    """unique base preimage of an odd target y with 3 not dividing y; None otherwise"""
    if y % 3 == 0: return None
    return (2 * y - 1) // 3 if y % 6 == 5 else (4 * y - 1) // 3

def main():
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 1 << 20
    out = []
    def P(*a):
        s = " ".join(str(x) for x in a); print(s, flush=True); out.append(s)
    t0 = time.time()
    L3 = math.log2(3)
    # ---- A
    bad = 0; cnt = 0; badk = 0
    for n in range(1, min(N, 1 << 18) + 1, 2):
        b, k = base_of(n)
        cnt += 1
        if not is_base(b) or S_pow(b, k) != n: bad += 1
        if k != (v2(3 * n + 1) - 1) // 2: badk += 1
    P(f"A. decomposition n = S^k(b), b a base (v2(3b+1) <= 2), k = floor((v2(3n+1)-1)/2): errors {bad}/{cnt}, k-formula errors {badk}")
    # child uniqueness: for bases c <= 2^14 and k <= 12: b0(S^k c) exists iff 3 does not divide S^k(c), and then U(b0(S^k c)) = S^k c with base map G = c, depth k
    bad = 0; cnt = 0; missing = 0
    for c in range(1, 1 << 14, 2):
        if not is_base(c): continue
        for k in range(0, 13):
            y = S_pow(c, k)
            b = b0(y)
            if b is None:
                missing += 1; continue
            cnt += 1
            if not is_base(b) or U(b) != y or base_of(U(b)) != (c, k): bad += 1
    P(f"   child-by-depth map c -> b0(S^k c): {cnt} children checked, errors {bad}; missing (3 | S^k c, exactly one residue of k mod 3) {missing}")
    # ---- B: base edges census
    edges = Counter(); nb = 0
    kmax = 0
    for b in range(3, N + 1, 2):
        if not is_base(b): continue
        a = v2(3 * b + 1); c, k = base_of(U(b))
        edges[(a, k)] += 1; nb += 1; kmax = max(kmax, k)
    P(f"B. base edges b <= {N}: {nb} bases; (a, k) fractions vs Haar (2/3, 1/3) x (3/4)(1/4)^k:")
    for a in (1, 2):
        row = " ".join(f"k={k}:{edges[(a,k)]/nb:.4f}({(2/3 if a == 1 else 1/3) * 0.75 * 0.25**k:.4f})" for k in range(0, 5))
        P(f"   a={a}: {row}  ...  max k {kmax}")
    # ---- C: chain code
    def chain(b):
        ks = []; aa = []
        while b != 1:
            a = v2(3 * b + 1); c, k = base_of(U(b))
            ks.append(k); aa.append(a); b = c
        return ks, aa
    def bill(ks, r, rho):
        c0 = math.log2(1 / (rho * (1 - r))); ck = math.log2(1 / r)
        return sum(c0 + k * ck for k in ks)
    params = {"codex (1/16, 1/2)": (1 / 16, 1 / 2), "critical (1/4, 1)": (1 / 4, 1.0 - 1e-12), "(1/4, 1/2)": (1 / 4, 1 / 2)}
    stats = {name: {"kraft": 0.0, "ratios": [], "paid": Counter(), "overdraft": []} for name in params}
    Tmax = (0, 0); billmax = {name: (0, 0) for name in params}
    chains_T = []
    sizes = {name: Counter() for name in params}
    for b in range(1, N + 1, 2):
        if not is_base(b): continue
        ks, aa = chain(b)
        T = len(ks); chains_T.append(T)
        if T > Tmax[0]: Tmax = (T, b)
        lb = math.log2(b) if b > 1 else 0.0
        # size drops along the chain
        drops = []
        x = b
        for i in range(T):
            c, k = base_of(U(x)); drops.append(math.log2(x) - math.log2(c)); x = c
        for name, (r, rho) in params.items():
            bl = bill(ks, r, rho)
            st = stats[name]
            st["kraft"] += 2.0 ** (-bl)
            if b > 1:
                st["ratios"].append(bl / lb)
                if bl > billmax[name][0]: billmax[name] = (bl, b)
            # per-edge size payability for s in 2,3,4,6,8 (first edge of b only, to count each base edge once)
            if b > 1:
                c0 = math.log2(1 / (rho * (1 - r))); ck = math.log2(1 / r)
                for s_ in (2, 3, 4, 6, 8):
                    if s_ * drops[0] >= c0 + ks[0] * ck: st["paid"][s_] += 1
                # chain overdraft with s = 6: min over prefixes of cumulative (s*drop - bill_i)
                bal = 0.0; mn = 0.0
                for i in range(T):
                    bal += 6 * drops[i] - (c0 + ks[i] * ck); mn = min(mn, bal)
                st["overdraft"].append(mn)
    nbases = sum(1 for b in range(1, N + 1, 2) if is_base(b))
    P(f"C. chain code over the {nbases} bases <= {N} (all rooted): chain length T = odd root time, mean {sum(chains_T)/len(chains_T):.2f}, max {Tmax}")
    for name, (r, rho) in params.items():
        st = stats[name]; rs = sorted(st["ratios"])
        q = lambda p: rs[min(len(rs) - 1, int(p * len(rs)))]
        P(f"   {name}: Kraft partial sum sum_(b<=N) 2^-bill(b) = {st['kraft']:.6f} (bound 1/(1-rho) = {1/(1-rho):.3f}); bill/log2(b): mean {sum(rs)/len(rs):.3f}, median {q(0.5):.3f}, 90% {q(0.9):.3f}, max {max(rs):.3f}; largest bill {billmax[name][0]:.1f} bits at b = {billmax[name][1]}")
    c0 = math.log2(32 / 15)
    P(f"   Haar prediction for the codex bill: per base change E[bill] = {c0:.3f} + 4 E[k] = {c0 + 4/3:.3f} bits, E[size drop] = E[a|base] + 2E[k] - log2 3 = {4/3 + 2/3 - L3:.3f} bits, so bill/log2 b -> {(c0 + 4/3)/(2 - L3):.2f}")
    # ---- D: size payability
    P("D. size weights g = b^-s: fraction of base edges (b <= N, first edge) paid, i.e. s*log2(b/G(b)) >= bill(edge):")
    for name in params:
        st = stats[name]
        P(f"   {name}: " + "  ".join(f"s={s_}: {st['paid'][s_]/(nbases-1):.4f}" for s_ in (2, 3, 4, 6, 8)) + f";  chains with no overdraft at any prefix (s = 6): {sum(1 for m in st['overdraft'] if m >= 0)/(nbases-1):.5f}, median worst overdraft {sorted(st['overdraft'])[len(st['overdraft'])//2]:.1f} bits")
    P(f"   exact per-edge drop: log2(b/G(b)) = a + 2k - log2 3 + log2((1 + 1/(3b)) / (1 + 2^(a+2k)(4^k-1)/(3*4^k... )))  [the carry; a in {{1,2}} on bases]; the critical bill (r=1/4, rho=1) is exactly 2 - log2 3 + 2k per edge, so its per-edge balance against s = 1 size is a - 2: zero on a = 2 edges, -1 bit on every a = 1 edge")
    # ---- E: kernel
    kp = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "05-knowledge", "results", "collatz_three_bit_sibling_flow_20261005.json")
    try:
        K = json.load(open(kp))
        bases = K["bases"]
        P(f"E. codex kernel: {len(bases)} entries; r = {K.get('r')}, rho = {K.get('rho')}; sample entry: {bases[0] if isinstance(bases, list) else list(bases.items())[0]}")
        # bills of the kernel bases under the codex parameters
        bl = []
        items = bases if isinstance(bases, list) else list(bases.values())
        for e in items:
            b = e["b"] if isinstance(e, dict) and "b" in e else (e["base"] if isinstance(e, dict) and "base" in e else (e[0] if isinstance(e, (list, tuple)) else None))
            if b is None: break
            ks, aa = chain(int(b)); bl.append((int(b), len(ks), bill(ks, 1/16, 1/2)))
        if bl:
            P(f"   kernel bases: chain lengths {min(x[1] for x in bl)}..{max(x[1] for x in bl)}; bills {min(x[2] for x in bl):.1f}..{max(x[2] for x in bl):.1f} bits; sum of 2^-bill over the kernel = {sum(2.0**-x[2] for x in bl):.6f} (the note reports sum_C g < 1 with its own path weights)")
            bl.sort(key=lambda x: -x[2]); P(f"   five largest kernel bills: {[(x[0], x[1], round(x[2],1)) for x in bl[:5]]}")
    except Exception as ex:
        P(f"E. kernel not parsed: {ex}")
    P(f"  [{time.time()-t0:.0f}s]")
    with open(__file__.replace(".py", ".out"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(out) + "\n")

def S_pow(b, k):
    for _ in range(k): b = 4 * b + 1
    return b

if __name__ == "__main__":
    main()
