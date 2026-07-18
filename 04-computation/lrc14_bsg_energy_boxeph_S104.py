#!/usr/bin/env python3
"""
BSG/PFR ATTACK on the inverse theorem -- the input is missing  (boxeph-2026-07-18-S104)

BSG/PFR/Freiman all need additive energy / small doubling as INPUT and output structure.
This script shows M<1/13 supplies NO such input:
 - additive energy DISCRIMINATES AP (E high, |C-C| minimal) from non-AP (E low, |C-C| large);
 - the full deep-well residue set has |R-R|=47 > 3k-4 (dimension 2 -- Freiman fails);
 - band-avoiding 12-sets (the LOCAL content of M<1/13) can have Sidon-like LOW energy.
So the crux factors: Half 1 (open) M<1/13 => high core energy [Diophantine->additive], and
Half 2 energy => AP [BSG/PFR]. BSG/PFR attack Half 2; Half 1 is the open crux (Tao n=12, S94).
"""
from collections import Counter


def energy(A):
    A = list(A)
    c = Counter(x + y for x in A for y in A)
    return sum(v * v for v in c.values())


def ddiff(A):
    A = list(A)
    return len(set(x - y for x in A for y in A))


n = 12
cube = n ** 3
print('additive energy / doubling DISCRIMINATES AP from non-AP cores (n=12, n^3=1728):')
cores = {
    'deep-well core {1..12} (AP)': list(range(1, 13)),
    'non-AP divisor-rich': [1, 2, 3, 4, 6, 8, 12, 24, 36, 48, 72, 144],
    'near-AP {1..11,13}': [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 13],
    'squares {1,4,...,144}': [i * i for i in range(1, 13)],
}
for name, C in cores.items():
    print(f'  {name:28}: |C-C|={ddiff(C):3}  E={energy(C):5}  E/n^3={energy(C)/cube:.3f}')

R = [14 * i for i in range(1, 13)] + [169]
print(f'\nfull deep-well residue set R: |R-R|={ddiff(R)} > 3k-4={3*13-4}  => DIMENSION 2 (Freiman fails on R)')

print('\nband-avoiding 12-sets inside [14,169] (LOCAL content of M<1/13) -- energy is NOT forced:')
for name, C in [('AP 14*{1..12}', [14 * i for i in range(1, 13)]),
                ('Sidon-like', [14, 15, 17, 21, 29, 45, 58, 73, 89, 120, 152, 169]),
                ('irregular', [14, 20, 28, 39, 55, 77, 90, 108, 130, 150, 160, 169])]:
    inband = all(14 <= x <= 169 for x in C)
    print(f'  {name:14} in[14,169]={inband}  |C-C|={ddiff(C):3}  E={energy(C):5}  E/n^3={energy(C)/cube:.3f}')
print('\n=> band-avoidance compatible with Sidon-like LOW energy => M<1/13 gives BSG/PFR no input.')
print('   The crux is Half 1 (M<1/13 => high energy), a Diophantine->additive bridge, OPEN (Tao n=12).')
