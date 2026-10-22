# Second-level congruence from the k-coloured Redei-Berge identity U_T(1^k) = sum_R (2k)^|R| k^(n-|V(R)|):
#   modulo 9, U_T(1,1,1) = 6 hc(T) for odd n (hc = directed Hamiltonian cycles) and 0 for even n, hence for the number
#   u_3(T) = c_3(T)/6 of UNORDERED covers by three directed paths:
#       u_3 = H + hc (mod 3)  (n odd),      u_3 = 0 (mod 3)  (n even).
# Exhaustive check over all tournaments n <= 8.  mac-mini-2026-10-22-sierp.
import subprocess
def tourns(n):
    for code in subprocess.run(['gentourng', str(n)], capture_output=True, text=True).stdout.split():
        if code.startswith('>'): continue
        out = [0] * n; p = 0
        for i in range(n):
            for j in range(i + 1, n):
                if code[p] == '1': out[i] |= 1 << j
                else: out[j] |= 1 << i
                p += 1
        yield out
for n in range(3, 9):
    bad = tot = 0
    for out in tourns(n):
        N = 1 << n
        dp = [[0] * n for _ in range(N)]
        for v in range(n): dp[1 << v][v] = 1
        for m in range(1, N):
            for v in range(n):
                c = dp[m][v]
                if not c: continue
                q = out[v] & ~m
                while q:
                    w = (q & -q).bit_length() - 1; q &= q - 1
                    dp[m | (1 << w)][w] += c
        Hs = [sum(dp[m]) for m in range(N)]; Hs[0] = 1
        full = N - 1
        c3 = 0
        for X in range(1, N):
            rest = full & ~X
            Y = rest
            while Y:
                Z = rest & ~Y
                if Z: c3 += Hs[X] * Hs[Y] * Hs[Z]
                Y = (Y - 1) & rest
        u3 = c3 // 6
        # Hamiltonian cycles: paths from vertex 0 covering all, closing back to 0
        hc = sum(dp[full][v] for v in range(n) if (out[v] >> 0) & 1 and v != 0) if False else 0
        # recompute paths that START at 0: run a DP restricted to start 0
        dp0 = [[0] * n for _ in range(N)]; dp0[1][0] = 1
        for m in range(1, N, 2):
            for v in range(n):
                c = dp0[m][v]
                if not c: continue
                q = out[v] & ~m
                while q:
                    w = (q & -q).bit_length() - 1; q &= q - 1
                    dp0[m | (1 << w)][w] += c
        hc = sum(dp0[full][v] for v in range(1, n) if (out[v] & 1))
        H = Hs[full]
        expect = (H + hc) % 3 if n % 2 else 0
        if u3 % 3 != expect: bad += 1
        tot += 1
    print(f"n={n}: {tot} tournaments, violations of u_3 = H + hc (odd n) / u_3 = 0 (even n) mod 3: {bad}", flush=True)
