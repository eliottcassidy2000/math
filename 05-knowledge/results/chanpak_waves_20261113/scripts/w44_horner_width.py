"""Width of the outer Horner digit.  Arrows N -> 3N + 2j from N = 1, j in a digit set, deleting every arrow INTO 7 or 21.
Which odd N <= LIMIT are unreachable?  (ARCHIVE.md of fresh_unit_prime_gluing_20261116 treats j in {0,1,2}.)"""
LIMIT = 200000
for digits in ([0, 1, 2], [0, 1, 2, 3], [0, 1, 2, 3, 4], [0, 1, 2, 3, 4, 5]):
    reach = {1}; stack = [1]
    while stack:
        n = stack.pop()
        for j in digits:
            m = 3 * n + 2 * j
            if m > LIMIT or m in (7, 21) or m in reach: continue
            reach.add(m); stack.append(m)
    miss = [n for n in range(1, LIMIT + 1, 2) if n not in reach]
    print(f"digits {digits}: unreachable odd N <= {LIMIT}: {len(miss)}; first {miss[:12]}")
# shield orbit of the hole 7 and cyclotomic values
print("7*3^k:", [7 * 3 ** k for k in range(5)], " Phi_3(q) = q^2+q+1 for q=2..9:", [q * q + q + 1 for q in range(2, 10)], " Phi_3(16) =", 16 * 16 + 16 + 1, "= 2730/10")
