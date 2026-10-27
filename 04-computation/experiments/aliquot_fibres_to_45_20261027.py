import math
N = 2200
s = [0]*(N+1)
for d in range(1, N//2+1):
    for m in range(2*d, N+1, d): s[m] += d
# complete fibres s^{-1}(m) for m <= 45: sources n satisfy n <= (m-1)^2 <= 1936 < N
fib = {m: [n for n in range(2, N+1) if s[n] == m] for m in range(1, 46)}
print("m: sources (complete; n <= (m-1)^2)")
for m in range(1, 46): print(m, fib[m])
print("untouchable m <= 45:", [m for m in range(2, 46) if not fib[m]])
