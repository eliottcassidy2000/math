# deg_extended_check.py (lane energy, S94): exact brute-force checks for Proposition 7 (B2') of PROOFS.txt.
#  For w-bit windows: all odd ordered pairs (x, x'), X = x - x', Xi = rev x - rev x', height M = max(|X|,|Xi|)/gcd.
#  (i)  N(K) = #{odd pairs : M <= K}  against the bound  4 * 6^s1 * K (1 + ln K)(4 K 2^(-s1) + 1)  for K = 2^(i+1) >= 2^s2;
#  (ii) Deg_K(p_u), Deg'_K(p_v) against the extended divisor bounds, for every point and every dyadic K (w <= 10);
#  (iii) S_w and S_w / 2^(beta1 w), S_w / 6^(w/2).
import sys, math
from math import gcd, log

def rev(x, n):
    r = 0
    for i in range(n):
        if (x >> i) & 1:
            r |= 1 << (n - 1 - i)
    return r

def run(w, full_deg):
    s1 = (w + 1) // 2; s2 = w // 2
    R = [rev(x, w) for x in range(1 << w)]
    Ncount = [0] * (w + 2)     # Ncount[i] = # odd pairs with M in [2^i, 2^(i+1))
    S = 0.0
    degU = {}                  # (a, b) -> list over i of # (a', b') with M < 2^(i+1)   (only if full_deg)
    degV = {}
    for x in range(1 << w):
        rx = R[x]
        for xp in range(1 << w):
            X = x - xp
            if X % 2 == 0:
                continue
            Xi = rx - R[xp]
            if Xi % 2 == 0:
                continue
            g = gcd(abs(X), abs(Xi))
            M = max(abs(X), abs(Xi)) // g
            i = M.bit_length() - 1
            Ncount[i] += 1
            S += 1.0 / M
    ok = True
    cum = 0
    rows = []
    for i in range(w):
        cum += Ncount[i]
        K = 1 << (i + 1)
        if K >= (1 << s2):
            bound = 4 * 6 ** s1 * K * (1 + log(K)) * (4 * K / 2 ** s1 + 1)
            if cum > bound:
                ok = False
            rows.append(f"i={i}: N={cum} bound={bound:.3g} ratio={cum/bound:.4f}")
    beta1 = math.log2(6) / 2 + 0.2
    print(f"w={w} s1={s1} s2={s2}  S_w={S:.4f}  S_w/6^(w/2)={S/6**(w/2):.4f}  S_w/2^(beta1 w)={S/2**(beta1*w):.5f}  N(K) bound ok: {ok}")
    print("   " + " | ".join(rows))
    if full_deg:
        # exact Deg_K over all points p_u = (a, b) and all (a', b'), by direct enumeration of the definition
        okD = True
        worst = 0.0
        for a in range(-(1 << s1) + 1, 1 << s1):
            if a == 0:
                continue
            for b in range(-(1 << s1) + 1, 1 << s1, 2):
                cnt = [0] * (w + 2)
                for ap in range(-(1 << s2) + 1, 1 << s2):
                    for bp in range(-(1 << s2) + 1, 1 << s2, 2):
                        X = ap * (1 << s1) + b
                        Xi = a * (1 << s2) + bp
                        g = gcd(abs(X), abs(Xi))
                        M = max(abs(X), abs(Xi)) // g
                        cnt[M.bit_length() - 1] += 1
                cum = 0
                for i in range(w + 1):
                    cum += cnt[i]
                    K = 1 << (i + 1)
                    bound = (2 * K / 2 ** s1 + 1) * (2 ** (s2 + 1) * (1 + log(K)) + K)
                    worst = max(worst, cum / bound)
                    if cum > bound:
                        okD = False
        print(f"   extended Deg_K bound for all points and dyadic K: {okD} (max ratio {worst:.4f})")
    sys.stdout.flush()

if __name__ == "__main__":
    for w in range(4, int(sys.argv[1]) + 1):
        run(w, w <= int(sys.argv[2]))
