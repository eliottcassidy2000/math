#!/usr/bin/env python3
"""223/233 meeting = (v, v+1) coalescence at v = 334; lifted family; 105 vs 233."""
def T(x): return x // 2 if x % 2 == 0 else (3 * x + 1) // 2
def Tn(x, n):
    for _ in range(n): x = T(x)
    return x
def odd_count(x, n):
    c = 0
    for _ in range(n): c += x % 2; x = T(x)
    return c
print('T(223) =', T(223), ' T^9(233) =', Tn(233, 9), ' odd steps of 233 in 9:', odd_count(233, 9))
print('T^6(334) =', Tn(334, 6), ' T^6(335) =', Tn(335, 6), ' first equal time:',
      next(t for t in range(1, 50) if Tn(334, t) == Tn(335, t)))
# THM-4581 6(c): merge of (v, v+1) by Terras time 6 depends only on v mod 64
cls = [r for r in range(64) if any(Tn(r + 64 * 10**6, t) == Tn(r + 64 * 10**6 + 1, t) for t in range(1, 7))]
print('residues v mod 64 with (v, v+1) merged by time 6:', len(cls), '/64; 334 mod 64 =', 334 % 64, 'in list:', 334 % 64 in cls)
# minimal joint family n = 223 + 31104 s, m = 233 + 32768 s
ok = all(T(223 + 31104 * s) == Tn(233 + 32768 * s, 9) + 1 and Tn(223 + 31104 * s, 7) == Tn(233 + 32768 * s, 15)
         for s in range(0, 3000))
print('family (223 + 31104 s, 233 + 32768 s): T(n) = T^9(m) + 1 and T^7(n) = T^15(m) for s < 3000:', ok)
ok2 = all(Tn(223 + 62208 * t, 7) == Tn(233 + 65536 * t, 15) == 425 + 118098 * t for t in range(0, 3000))
print('prior family (223 + 62208 t, 233 + 65536 t) meets at 425 + 118098 t for t < 3000:', ok2)
# 105 vs 233: equal-time merge?
def terras_to_cycle(x):
    seq = [x]
    while x not in (1, 2):
        x = T(x); seq.append(x)
    return seq
a, b = terras_to_cycle(105), terras_to_cycle(233)
print('Terras times to {1,2}: 105 ->', len(a) - 1, 'ending at', a[-1], '; 233 ->', len(b) - 1, 'ending at', b[-1])
L = max(len(a), len(b)) + 4
eq = next((t for t in range(L) if Tn(105, t) == Tn(233, t)), None)
print('first t with T^t(105) = T^t(233):', eq)
print('T^7(105), T^7(233) =', Tn(105, 7), Tn(233, 7), ' diff =', Tn(233, 7) - Tn(105, 7))
