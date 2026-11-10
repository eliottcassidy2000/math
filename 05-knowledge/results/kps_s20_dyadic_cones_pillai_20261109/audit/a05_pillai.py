# a05: THM-5109 (E1).
#  (a) exhaustive exact search of 1 <= 2^x - 3^y <= 2^y for x <= 400 (all y with 3^y < 2^x);
#  (b) where does Ellison's inequality |2^x - 3^y| > 2^x e^(-x/10) (all y) fail, for 1 <= x <= 20000?
#      (quoted: for x >= 12 it fails only at x = 13, 14, 16, 19, 27);
#  (c) the constants of the inequality chain in PROOF.md (E1).
import math
sols = []
for x in range(1, 401):
    y = 0
    while 3 ** y < 2 ** x:
        d = 2 ** x - 3 ** y
        if 1 <= d <= 2 ** y: sols.append((x, y, d))
        y += 1
print('(a) solutions of 1 <= 2^x - 3^y <= 2^y with x <= 400:', sols)
strict = [(x, y, d) for (x, y, d) in sols if d < 2 ** y]
print('    of which strict (2^x - 3^y < 2^y):', strict)

fails = []; minmargin = None
for x in range(1, 20001):
    # closest powers of 3 on both sides
    y0 = int(x * math.log(2) / math.log(3))
    worst = None
    for y in range(max(0, y0 - 2), y0 + 3):
        d = abs(2 ** x - 3 ** y)
        if d == 0: continue
        # compare ln d - x ln 2  vs  -x/10   (math.log is accurate for big ints)
        margin = (math.log(d) - x * math.log(2)) + x / 10.0
        if worst is None or margin < worst: worst = margin
    if worst <= 0: fails.append(x)
    elif x >= 28 and (minmargin is None or worst < minmargin[1]): minmargin = (x, worst)
print('(b) x in [1, 20000] where Ellison\'s inequality fails for some y:', fails)
print('    restricted to x >= 12:', [x for x in fails if x >= 12])
print('    smallest log-margin for 28 <= x <= 20000: x = %d, ln(|2^x-3^y|/2^x) + x/10 = %.4f' % minmargin)
c1 = 1 - 1 / (10 * math.log(2)); c2 = math.log(2) / math.log(3)
print('(c) 1 - 1/(10 ln 2) = %.6f (PROOF uses > 0.855);  log 2 / log 3 = %.6f (PROOF uses < 0.631)' % (c1, c2))
