import csv,math
from pathlib import Path
p=Path(__file__).parent
strata=[]
for n in range(1,6):
 h={int(r['H']):int(r['count']) for r in csv.DictReader(open(p/f'hist-base4-n{n}.csv'))}
 assert sum(h.values())==4**(n*(n-1)//2)
 assert sum(k*v for k,v in h.items())==math.factorial(n)*2**(n-1)*4**(n*(n-1)//2-n+1)
 for l in open(p/f'unrestricted-{n}.txt'):
  if l.startswith('stratum '):strata.append(tuple(map(int,l.split()[1].split(','))))
 print(f'PASS unrestricted n={n}: full graph and path-incidence totals; H7={h.get(7,0)}, H21={h.get(21,0)}')
ref=[tuple(map(int,r)) for r in list(csv.reader(open(p.parent/'census/general_census.csv')))[1:]]
assert sorted(strata)==sorted(ref)
print('PASS every scanner unrestricted target stratum (order,count,absent pairs,reciprocal pairs,strongness,frequency)')
