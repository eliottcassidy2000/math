#!/usr/bin/env python3
"""Exact arithmetic and direct permutation checks; no floating-point inputs."""
import itertools,json,pathlib
P=pathlib.Path(__file__).resolve().parent
W=json.loads((P/'witnesses.json').read_text())
def brute(w):
 n=w['n'];a=[[False]*n for _ in range(n)];k=0
 for i in range(n):
  for j in range(i+1,n):
   a[i][j]=bool(w['code']>>k&1);a[j][i]=not a[i][j];k+=1
 return sum(all(a[x][y] for x,y in zip(p,p[1:])) for p in itertools.permutations(range(n)))
for w in W:assert brute(w)==w['h'],w
checks=[]
def check(label,lhs,rhs):
 assert lhs>=rhs,label
 checks.append(dict(label=label,lhs=lhs,rhs=rhs))
for n,d in [(9,75),(10,125),(11,225)]:check('Busch base '+str(n),d**5,11**n)
check('Busch step',5**5,11**3)
for n,d in [(5,11),(6,19),(7,29),(8,47)]:check('small good threshold '+str(n),d**5,11**n)
for n,d in [(5,15),(6,25),(6,27),(6,45),(8,289)]:check('replacement '+str(d),d**5,11**n)
cost={1:0,3:3,5:4,9:5,15:5,17:6,45:6,51:7,85:7,153:7,255:8,765:11}
assert set(cost)=={3**r*5**t*17**u for r in range(3) for t in range(2) for u in range(2)}
for d,n in cost.items():
 if d not in (1,765):assert any(w['h']==d and w['n']==n for w in W)
 check('residual '+str(d),11**3*d**5,3**5*11**n)
assert 765==3*255 and cost[765]==cost[3]+cost[255]
print(json.dumps({'status':'PASS','direct_permutation_witnesses':len(W),'integer_checks':checks,'residual_costs':cost},indent=2))
