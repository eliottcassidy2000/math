# 27, 63, 1473, 1729 in the signed descent, the Collatz merge at 182, and one-apex witnesses.
def D(h): return (h + 1) // 2 if h % 4 == 1 else ((3 * h + 3) // 4 if h % 8 == 3 else (3 * h - 1) // 4)
def orbit(h):
    o = [h]
    while o[-1] > 3: o.append(D(o[-1]))
    return o
for h in (27, 63, 47, 55, 1473, 1729, 2945):
    print("count %d, half-count %d, signed descent: %s" % (h, (h - 1) // 2, orbit(h)))
def T(n): return n // 2 if n % 2 == 0 else (3 * n + 1) // 2
def to(n, k):
    o = [n]
    for _ in range(k): o.append(T(o[-1]))
    return o
print("shortcut orbits:", to(27, 10), to(63, 8), to(62, 8))
assert all(to(27 + 1024 * j, 10)[-1] == 182 + 6561 * j == to(63 + 2304 * j, 8)[-1] for j in range(200))
print("27 + 1024 j and 63 + 2304 j both reach 182 + 6561 j: checked for j < 200")
def h_word(w):
    h, A, p = 1, 0, 0
    for c in w:
        if c == '1': A, p = 2 * A + 1, 1
        else: h, A, p = h + A + p, 2 * A, 0
    return h
for w in ("100010", "1001100", "11101001000"): print("one-apex word %s: H = %d" % (w, h_word(w)))
def J(a):
    b = [i for i in range(a.bit_length()) if a >> i & 1]
    return sum(1 if y - x == 1 else 2 ** (y - x - 2) for i, x in enumerate(b) for y in b[i + 1:])
print("one-apex half-counts a - J(a) equal to 736 (count 1473) for a < 2^13:", [a for a in range(1, 1 << 13) if a - J(a) == 736])
# the two rays from the triangle: 2^a + 1 = 3^b
print("2^a + 1 = 3^b for a, b <= 400:", [(a, b) for a in range(1, 401) for b in range(1, 401) if 2 ** a + 1 == 3 ** b])
