"""THM-4890 (D) arithmetic: m = floor(n/2) - 2, C = 2^m, P = 2^(n-m-1), B = P/C.
bracket(n) = 2 + 5B/C + 3(4 ln B + 10.58) + 3(B + 1 + ln B)/C  (claimed <= 76),
total(n)   = 1 + 2B + bracket(n)                                  (claimed <= 109; type 3 <= 2P^2 = 2B PC).
Also the strip-lemma per-pair bound for W = C: (C/beta + 1)(2 beta C/h' + 3) <= 2C^2/P + 3C/beta + 5 when h' >= P = BC."""
from fractions import Fraction
import math

worst = (0, None)
rows = []
for n in range(12, 81):
    m = n // 2 - 2
    C = 2 ** m; P = 2 ** (n - m - 1); B = P // C
    br = 2 + 5 * B / C + 3 * (4 * math.log(B) + 10.58) + 3 * (B + 1 + math.log(B)) / C
    tot = 1 + 2 * B + br
    rows.append((n, m, C, P, B, br, tot))
    if br > worst[0]:
        worst = (br, n)
for r in rows[:8]:
    print('n=%d m=%d C=%d P=%d B=%d bracket=%.4f total=%.4f' % r)
print('max bracket over n = 12..80: %.4f at n = %d (claimed <= 76); max total %.4f (claimed <= 109)'
      % (worst[0], worst[1], max(r[6] for r in rows)))
# per-pair strip bound algebra: (C/b + 1)(2 b C/h' + 3) = 2C^2/h' + 3C/b + 2bC/h' + 3 <= 2C^2/P + 3C/b + 2 + 3 for b <= B, h' >= P
ok = True
for B in (8, 16):
    for C in (16, 32, 64, 128):
        P = B * C
        for b in range(1, B + 1):
            for hp in (P, P + 1, 2 * P - 1):
                lhs = (Fraction(C, b) + 1) * (Fraction(2 * b * C, hp) + 3)
                rhs = Fraction(2 * C * C, P) + Fraction(3 * C, b) + 5
                if lhs > rhs:
                    ok = False
print('per-pair strip bound (C/beta+1)(2 beta C/h\'+3) <= 2C^2/P + 3C/beta + 5 for beta <= B, h\' >= P: %s' % ok)
