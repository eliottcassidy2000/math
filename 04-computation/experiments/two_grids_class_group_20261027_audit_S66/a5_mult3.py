"""Post-hoc check of HYP-9429's 'Pattern' sentence: how many misses are multiples of 3, and the naive tail probability."""
from fractions import Fraction
from math import comb
pop = [6] + [n for n in range(51, 201) if n % 8 != 5]
m3 = [n for n in pop if n % 3 == 0]
miss9 = [6, 54, 60, 75, 78, 102, 124, 140, 186]
miss7 = [6, 54, 60, 75, 78, 124, 140]
print('population (6 and non-obstructed 51..200):', len(pop), 'multiples of 3:', len(m3), f'({len(m3)/len(pop):.3f})')
for name, S in (('nine misses at B<=60', miss9), ('seven persistent misses', miss7)):
    k = sum(1 for n in S if n % 3 == 0)
    N, K, d = len(pop), len(m3), len(S)
    hyper = sum(Fraction(comb(K, j) * comb(N - K, d - j), comb(N, d)) for j in range(k, d + 1))
    p = Fraction(1, 3)
    binom = sum(comb(d, j) * p ** j * (1 - p) ** (d - j) for j in range(k, d + 1))
    print(f'{name}: {k} of {d} are multiples of 3 ({[n for n in S if n % 3 == 0]}); '
          f'hypergeometric tail P(>= {k}) = {float(hyper):.4f}; binomial(1/3) tail = {float(binom):.4f}')
