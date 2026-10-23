#!/usr/bin/env python3
import json,itertools,hashlib,pathlib,math
from check import context_coeff,check,tour,profile

def build(B):
 r=len(B); q=[[0]*(2*r+1) for _ in range(2*r+1)]
 for a,b in itertools.combinations(range(1,2*r+1),2):
  if (a<=r)==(b<=r):q[a][b]=1
 for i in range(r):
  q[0][i+1]=1;q[r+1+i][0]=1
  for j in range(r):q[i+1][r+1+j]=B[i][j];q[r+1+j][i+1]=1-B[i][j]
 return q
bad=[[1,0,0],[1,0,0],[0,1,1]]
good=[[1,0,0],[0,1,0],[1,0,1]]
qs=[build(B) for B in [bad,good]]
cs=[context_coeff(q) for q in qs]
check([sum(row) for row in qs[0]]==[sum(row) for row in qs[1]],'labelled scores')
check([c[4] for c in cs]==[0,6],'top kernel coefficients')
triangle_checks=0
for n in range(3,6):
 for bits in range(1<<(n*(n-1)//2)):
  a=tour(n,bits);f=profile(a)
  c3=sum(all(sum(a[v][w] for w in tri if w!=v)==1 for v in tri) for tri in itertools.combinations(range(n),3))
  expected=math.comb(n,3)+2*c3+3*math.comb(n,4)
  check(f[n-2]==math.factorial(n-2)*expected,'penultimate profile triangle formula')
  lhs=12*f[n-2]-6*(n-2)*f[n-1]+3*math.comb(n-1,2)*f[n]
  check(lhs==24*math.factorial(n-2)*c3-math.factorial(n)*(n-2),'centered polynomial triangle lattice')
  triangle_checks+=1
# Exhaustive lower-order incidence-margin minimality, not generic quotient-score minimality.
margin_tables=0
for r in [1,2]:
 seen={}
 for bits in range(1<<(r*r)):
  B=[[bits>>(r*i+j)&1 for j in range(r)] for i in range(r)]
  key=(tuple(map(sum,B)),tuple(sum(B[i][j] for i in range(r)) for j in range(r)))
  per=sum(all(B[i][p[i]] for i in range(r)) for p in itertools.permutations(range(r)))
  check(key not in seen or seen[key]==per,'small shore margins determine permanent')
  seen[key]=per;margin_tables+=1
result={'triangle_checks':triangle_checks,'small_margin_tables':margin_tables,'status':'PASS','matrices':[bad,good],'labelled_scores':[sum(row) for row in qs[0]],'word_coefficient_vectors':cs,'top_a4':[c[4] for c in cs],'bad_Hall_set_rows':[0,1],'bad_Hall_neighborhood_columns':[0],'script_sha256':hashlib.sha256(pathlib.Path(__file__).read_bytes()).hexdigest()}
print(json.dumps(result,sort_keys=True,indent=2))
