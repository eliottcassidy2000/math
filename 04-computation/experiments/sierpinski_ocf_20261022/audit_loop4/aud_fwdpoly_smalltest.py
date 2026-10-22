# AUDIT: brute-force test of aud_fwdpoly.c on random tournaments n = 3..8 (all orderings), several x and primes.
import itertools, random, subprocess, sys
random.seed(20261022)
bad = 0; tests = 0
for n in range(3, 9):
    for rep in range(4):
        code = ''.join(random.choice('01') for _ in range(n * (n - 1) // 2))
        beats = [[0] * n for _ in range(n)]; k = 0
        for i in range(n):
            for j in range(i + 1, n):
                if code[k] == '1': beats[i][j] = 1
                else: beats[j][i] = 1
                k += 1
        cnt = [0] * n
        for s in itertools.permutations(range(n)):
            cnt[sum(beats[s[i]][s[i + 1]] for i in range(n - 1))] += 1
        for p in (1000000007, 2305843009213693951, 97):
            xs = [0, 1, 2, 5, 13]
            out = subprocess.run(['./aud_fwdpoly', str(p)] + [str(x) for x in xs], input=code + '\n', capture_output=True, text=True, check=True).stdout.split('\n')
            got = {int(a): int(b) for a, b in (l.split() for l in out if l.strip())}
            for x in xs:
                want = sum(c * x ** e for e, c in enumerate(cnt)) % p
                tests += 1
                if got[x % p] != want: bad += 1; print('MISMATCH', n, code, p, x, got[x % p], want)
print(f"aud_fwdpoly brute-force test: {tests} (tournament, p, x) cases, {bad} mismatches")
