#!/usr/bin/env python3
"""forcing_table.py: merge strongB_3..9.out into the all-order table for counts h < 125 (one nontrivial strong component):
digit 0 / digit 1 at some order (2^k B = h + 1) / digit 2 / rooted doubling child at some order (charge 0 or 2^k B = h - 1,
k >= 1); the primes without a tower-forced 3p + 2; the counts without a rooted doubling child; wave 50's P111 in the
rooted-append sense (primes p = 1 mod 4 below 250: does h = (p+1)/2 have a rooted doubling child?)."""
import re
def is_prime(n): return n > 1 and all(n % d for d in range(2, int(n**0.5)+1))
rows = {}; Bsets = {}
for k in range(3, 10):
    for l in open(f'strongB_{k}.out'):
        if not l.startswith('  '): continue
        h, rest = l.split(':', 1); h = int(h); flags, bv = rest.split('|'); f = [int(x) for x in flags.split()]
        cur = rows.get(h, [0,0,0,0]); rows[h] = [max(a, b) for a, b in zip(cur, f)]
        Bsets.setdefault(h, set()).update(int(x) for x in bv.split() if x.isdigit())
# recompute the all-order flags from the B sets (digit1any: 2^k B = h+1; charge0: a_t = 0 known only from flags; 2^k B = h-1)
res = {}
for h in rows:
    Bs = Bsets[h]
    d1 = any((B << k) == h + 1 for B in Bs for k in range(0, 12))
    dbl = rows[h][3] or any((B << k) == h - 1 for B in Bs for k in range(1, 12))
    res[h] = (rows[h][0], int(d1), rows[h][2], int(dbl))
att = sorted(h for h in res if h < 125)
print("counts h < 125 with a strong realiser of order <= 9:", att)
print("missing below 125 (holes or only decomposable):", [h for h in range(1, 125, 2) if h not in res])
print("TABLE h: digit0 digit1any digit2 doublingany")
for h in att: print(f"  {h}: {res[h][0]} {res[h][1]} {res[h][2]} {res[h][3]}  B-values {sorted(Bsets[h])[:12]}{' ...' if len(Bsets[h]) > 12 else ''}")
print("primes p < 125 WITHOUT a tower-forced 3p+2 (no 2^k B_t = p+1):", [h for h in att if is_prime(h) and not res[h][1]])
print("primes p < 125 WITH a tower-forced 3p+2:", [h for h in att if is_prime(h) and res[h][1]])
print("counts h < 125 without digit 0 (no B_t = h):", [h for h in att if not res[h][0]])
print("counts h < 125 without digit 2:", [h for h in att if not res[h][2]])
print("counts h < 125 without a rooted doubling child at any order:", [h for h in att if not res[h][3]])
print("sterile forks h = 11 mod 18 below 125:", [h for h in range(11, 125, 18)], "; THM-5272 doubling failures: 11, 29, 35")
p1 = [p for p in range(5, 250) if is_prime(p) and p % 4 == 1]
print("P111 (rooted-append sense): primes p = 1 mod 4 < 250 whose (p+1)/2 has NO rooted doubling child:", [p for p in p1 if (p+1)//2 in res and not res[(p+1)//2][3]], "; (p+1)/2 not a count:", [p for p in p1 if (p+1)//2 not in res])
