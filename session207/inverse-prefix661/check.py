import itertools,json,pathlib,math
P=pathlib.Path(__file__).resolve().parent
rows=json.load(open(P/'census.json'))['orders']
m={};strong_min={}
for r in rows:
 for h in r['spectrum']:m.setdefault(h,r['n'])
 for h in r['strong_spectrum']:strong_min.setdefault(h,r['n'])
for line in (P/'fill.out').read_text().splitlines():
 h,code,base,pat=map(int,line.split());n=9;a=[[False]*n for _ in range(n)];k=0
 for i in range(n):
  for j in range(i+1,n):a[i][j]=bool(code>>k&1);a[j][i]=not a[i][j];k+=1
 brute=sum(all(a[p[i]][p[i+1]] for i in range(n-1)) for p in itertools.permutations(range(n)))
 assert brute==h and h not in m
 reach=[{i}|{j for j in range(n) if a[i][j]} for i in range(n)]
 for k in range(n):
  for i in range(n):
   if k in reach[i]:reach[i]|=reach[k]
 assert all(len(r)==n for r in reach)
 m[h]=9;strong_min[h]=9
 print('PASS strong n9 witness',h,code,'by all 9! permutations')
assert set(m)==set(range(1,662,2))-{7,21}
assert m==strong_min
result={'scope':'All integers h in [1,661]; even h and 7,21 are impossible. For each listed odd h, realization orders are exactly n>=m(h).','minimum_orders':{str(h):m[h] for h in sorted(m)}}
(P/'minimum_orders.json').write_text(json.dumps(result,indent=2)+'\n')
print('PASS complete exact inverse table for all 329 attainable values in [1,661].')
# Direct controls for the cyclic fan construction H=2^k+1.
for k in range(1,6):
 n=k+2;a=[[False]*n for _ in range(n)]
 for i in range(k):
  for j in range(i+1,k):a[i][j]=True
  a[i][k]=True;a[k+1][i]=True
 a[k][k+1]=True
 h=sum(all(a[p[i]][p[i+1]] for i in range(n-1)) for p in itertools.permutations(range(n)))
 assert h==2**k+1
print('PASS cyclic-fan identity controls k=1..5.')
