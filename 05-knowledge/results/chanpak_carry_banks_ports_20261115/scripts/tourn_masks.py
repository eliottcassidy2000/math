
import sys
# gentourng default output (upper triangle row by row, '1' = arc i -> j for i < j) to lines of out-masks
n = int(sys.argv[1])
for line in sys.stdin:
    s = line.strip()
    if len(s) != n * (n - 1) // 2 or set(s) - set('01'): continue
    out = [0] * n; k = 0
    for i in range(n):
        for j in range(i + 1, n):
            if s[k] == '1': out[i] |= 1 << j
            else: out[j] |= 1 << i
            k += 1
    print(' '.join(map(str, out)))
