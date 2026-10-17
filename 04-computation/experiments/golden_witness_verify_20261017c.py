#!/usr/bin/env python3
"""Independent positive-witness replay, without SAT or the producer imports."""
import hashlib,json,math
from pathlib import Path
R=Path('04-computation/experiments/golden_tail_20261017c')
checks=0
def check(p):
 global checks
 checks+=1
 if not p:raise RuntimeError('independent witness failure')
def add(a,b):return tuple((x+y)%19 for x,y in zip(a,b))
def mul(a,b):return ((a[0]*b[0]-a[1]*b[1])%19,(a[0]*b[1]+a[1]*b[0])%19)
def scale(a,b):return (a*b[0]%19,a*b[1]%19)
def norm(a):return (a[0]**2+a[1]**2)%19
def diff(a,b):return add(a,scale(-1,b))
d=json.loads((R/'boundary.json').read_bytes());s=json.loads((R/'circle.json').read_bytes())
check(s['source_sha256']==hashlib.sha256((R/'boundary.json').read_bytes()).hexdigest())
E=[(0,1),(1,2),(2,3),(3,4),(4,0),(0,5),(5,6),(6,7),(7,0)]
res=[(x,y) for x in range(19) for y in range(19)]
edges=[(a,b) for a in range(361) for b in range(a+1,361) if norm(diff(res[a],res[b]))==1]
check(len(edges)==3610)
for b in d['branches']:
 a=tuple(b['alpha']);a2=mul(a,a);a4=mul(a2,a2);phi=(1+2*a[0])%19;root5=(2*phi-1)%19
 check(mul(a4,a)==(1,0) and a!=(1,0))
 imag=(0,b['sqrt11']*root5%19)
 ports=[(0,0),(1,0),add((1,0),a),add(add((1,0),a),a2),scale(-1,a4),scale(12,add((3,0),scale(-1,imag))),(15,0),scale(12,add((3,0),imag))]
 check(ports==list(map(tuple,b['points'])))
seen=set();case_counts=[]
for obj in d['cyclotomic_templates']:
 pat=tuple(obj['pattern']);seen.add(pat);b=d['branches'][obj['branch']];aa,bb,t=obj['projection']
 c=[obj['permutation'][[4,3,2,4,3,2,1,0,4,1,0,4,3,2,1,0,2,1,0][(aa*x+bb*y+t)%19]] for x,y in res]
 check(all(c[a]!=c[b] for a,b in edges))
 for i,p in enumerate(b['points']):check(all(c[a]!=pat[i] for a,z in enumerate(res) if norm(diff(z,p))==1))
 check(all(pat[a]!=pat[b] for a,b in E))
for obj in s['answers']:
 pat=tuple(obj['pattern']);check(pat not in seen);seen.add(pat);w=obj['witness'];b=d['branches'][w['branch']];c=w['colours']
 check(len(c)==369 and c[361:]==list(pat))
 check(all(c[a]!=c[b] for a,b in edges))
 for i,p in enumerate(b['points']):check(all(c[a]!=pat[i] for a,z in enumerate(res) if norm(diff(z,p))==1))
 check(all(pat[a]!=pat[b] for a,b in E))
 case_counts.append([list(pat),obj['positive_templates']])
check(len(seen)==570)
for pat in seen:
 check(pat[0]==0 and all(pat[i]<=max(pat[:i])+1 for i in range(1,8)))
check(sum(math.prod(range(6-(max(p)+1)+1,7)) for p in seen)==(5**5-5)*(5**4+5)//6)
print(json.dumps(dict(status='all570 complete positive field-boundary witnesses independently replayed',checks=checks,case_counts=case_counts),sort_keys=True))
