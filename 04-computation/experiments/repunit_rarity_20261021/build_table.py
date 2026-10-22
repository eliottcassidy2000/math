# Complete strong-class table sigma(j,d) = #strong classes of order j with H = d, all orders, odd d < 375.
# Orders 1..9 from the mac-mini census (tour_census output: H C strong aut anti c3); orders >= 10 from S50's
# strong_ears_375.out (complete for d < 375: f(10)=125, f(11)=225, f(12)=375).  Cross-check on orders 7..9.
import json, gzip, os
from collections import defaultdict
base='../tournament_pseudoprimes_20261021/'
sig=defaultdict(lambda: defaultdict(int))
sig[1][1]=1   # K1
for n in range(3,10):
    p=base+f'census_n{n}.txt'
    f=open(p) if os.path.exists(p) else gzip.open(p+'.gz','rt')
    for l in f:
        H,C,s,A,AA,c3=map(int,l.split())
        if s and H<375: sig[H][n]+=1
S50={int(k):{int(a):b for a,b in v.items()} for k,v in json.load(open('strong_polys_lt375.json')).items()}
mism=[]
for d,terms in S50.items():
    for j,c in terms.items():
        if j<=9:
            if sig[d][j]!=c: mism.append((d,j,sig[d][j],c))
        else: sig[d][j]=c
print('cross-check orders <= 9 vs S50 for d >= 125: mismatches', len(mism), mism[:5])
tab={d:dict(sig[d]) for d in range(1,375,2)}
json.dump({str(k):v for k,v in tab.items()},open('strong_table_lt375.json','w'))
print('P_d(1) for small d:',[(d,sum(tab[d].values())) for d in range(1,60,2)])
