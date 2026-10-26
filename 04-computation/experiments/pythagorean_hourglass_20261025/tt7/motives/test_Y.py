import re, sys
exec(open('identify_new.py').read().split('results = []')[0].replace('print(f"{len(qs)} primes', 'pass  # print(f"{len(qs)} primes'))
Yv = {}
for l in open('triple_Y_7_299.out'):
    m = re.match(r'q=(\d+) Y=(-?\d+)', l); Yv[int(m.group(1))] = int(m.group(2))
qq = [q for q in qs if q in Yv]
idx = [qs.index(q) for q in qq]
B2 = [[c[i] for i in idx] for c in base]; N2 = [[c[i] for i in idx] for c in newcols]
Ycol = [Yv[q] for q in qq]
rb = rank(B2); rbn = rank(B2 + N2)
print(f"{len(qq)} primes; new dim {rbn - rb}; Y in base span: {rank(B2 + [Ycol]) == rb}; Y in base+new: {rank(B2 + N2 + [Ycol]) == rbn}; "
      f"new dim after adding Y: {rank(B2 + [Ycol] + N2) - rank(B2 + [Ycol])}")
H18 = [S[18][q] for q in qq]
print("S_H18 - Y in base span:", rank(B2 + [[a - b for a, b in zip(H18, Ycol)]]) == rb)
