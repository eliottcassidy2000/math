#!/usr/bin/env python3
"""Referee: Antihydra odd-orbit counter forms (index range of the v_3 form) and A4 halting set on TM-verified rules."""
import sys, os
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from tm_rules_check import a4_rule, HALT

def v2(n):
    return (n & -n).bit_length() - 1
def v3(n):
    k = 0
    while n % 3 == 0:
        n //= 3; k += 1
    return k

# Antihydra (public rule): a=8, b=0; even: a=3a/2, b+=2; odd: halt if b==0 else a=(3a-1)/2, b-=1
a, b = 8, 0
odd_seen = []
steps = 0
while len(odd_seen) < 3000:
    if a % 2 == 0:
        a, b = 3 * a // 2, b + 2
    else:
        odd_seen.append((a, b))      # counter on arrival at the odd term
        if b == 0:
            print('halt'); break
        a, b = (3 * a - 1) // 2, b - 1
o = [x[0] for x in odd_seen]; c = [x[1] for x in odd_seen]
v = [v2(3 * x - 1) for x in o]
ok_twisted = all(o[i + 1] == 3 ** (v[i] - 1) * (3 * o[i] - 1) // 2 ** v[i] for i in range(len(o) - 1))
ok_counter = all(c[i] == 6 + sum(2 * v[l] - 3 for l in range(i)) for i in range(len(o)))
ok_v3_from1 = all(c[i] == 6 + sum(2 * v3(o[l]) - 1 for l in range(1, i + 1)) for i in range(len(o)))
ok_v3_from0 = all(c[i] == 6 + sum(2 * v3(o[l]) - 1 for l in range(0, i + 1)) for i in range(len(o)))
print('o_0 =', o[0], 'counter at o_0 =', c[0], '| twisted map ok:', ok_twisted, '| counter 6+sum_{l<i}(2v_l-3) ok:', ok_counter)
print('v3 form with l = 1..i:', ok_v3_from1, '| v3 form with l = 0..i (as written, l <= i):', ok_v3_from0,
      '| v3(o_0) =', v3(o[0]), ' c_1 =', c[1], ' 6+sum_{l=0}^{1} =', 6 + (2 * v3(o[0]) - 1) + (2 * v3(o[1]) - 1))

# A4 halting set on my rules: from C(a,1), follow rules to next visit or HALT
def follow(a):
    q = (a, 1)
    q = a4_rule(q)
    for _ in range(10 ** 6):
        if q == HALT:
            return 'HALT'
        if q[1] == 1:
            return q[0]
        q = a4_rule(q)
    return None
H = set()
m = 0
while (1 << (m + 2)) - m - 4 < 20000:
    H.add((1 << (m + 2)) - m - 4); m += 1
bad = [a for a in range(20000) if (follow(a) == 'HALT') != (a in H)]
print('A4: from C(a,1) the next event is HALT iff a = 2^(m+2)-m-4, a < 20000:', not bad, bad[:5])
# visits from the start and the Hydra form f(a) = floor(3a/2)+2 on C(a,3)
okf = all(a4_rule((x, 3)) == ((3 * x) // 2 + 2, 1 if x % 2 == 0 else 3) for x in range(1, 5000))
print('A4: C(a,3) -> C(floor(3a/2)+2, 1 if a even else 3) for 1 <= a < 5000:', okf)
