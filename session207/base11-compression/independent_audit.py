import json,pathlib,itertools
P=pathlib.Path(__file__).resolve().parent
for r in json.loads((P/'witnesses.json').read_text()):
 n=r['n'];c=r['code'];a=[[False]*n for _ in range(n)];k=0
 for i in range(n):
  for j in range(i+1,n):
   a[i][j]=bool(c>>k&1);a[j][i]=not a[i][j];k+=1
 h=sum(all(a[p[i]][p[i+1]] for i in range(n-1)) for p in itertools.permutations(range(n)))
 assert h==r['h']
print('PASS all13 physical witnesses by full vertex permutations')
for s,d in [(9,75),(10,125),(11,225),(5,11),(6,19),(7,29),(8,47),(5,15),(6,25),(6,27),(6,45),(8,289)]:
 assert d**5>=11**s,(s,d)
assert 5**5>=11**3
print('PASS large-order induction and all finite normalization thresholds')
costs={1:0,3:3,5:4,9:5,15:5,17:6,45:6,51:7,85:7,153:7,255:8,765:11}
assert set(costs)=={3**r*5**t*17**u for r in range(3) for t in range(2) for u in range(2)}
for d,s in costs.items():assert 3**5*11**s<=11**3*d**5,(d,s)
assert costs[765]==3+8 and 765==3*255
# Slope comparison without floating logs: 5/log11 < (5/2)/log3 iff 9<11.
assert 3**2<11
print('PASS all12 residual inequalities, physical join765, and exact slope comparison')
