#!/usr/bin/env python3
"""Independent audit: exhaustive binary words, signed edge reconstruction.
Does not import primary search code. Rebuilds exact realized signature tables,
checks all valuation identities on binary words, and verifies cycle outputs.
"""
import itertools,json,math,argparse
from fractions import Fraction

def primes(d):
 fs=[]
 for p in range(2,math.isqrt(d)+2):
  e=0
  while d%p==0:d//=p;e+=1
  if e:fs.append((p,e))
 if d>1:fs.append((d,1))
 return fs

def first_change(w,i,j):
 for t in range(len(w)):
  if w[(i+t)%len(w)]!=w[(j+t)%len(w)]:return t
 raise ValueError('identical rotations')

def word_tables(N):
 tables={};union=set();identity_checks=0
 for n in range(3,N+1):
  rows=set();primitive_count=0
  proper=[d for d in range(1,n) if n%d==0]
  for number in range(2**(n-1)):
   w=[0]+[(number>>j)&1 for j in range(n-1)]
   if any(all(w[i]==w[i%d] for i in range(n)) for d in proper):continue
   primitive_count+=1
   def q(i,j):return 0 if w[i%n]==w[j%n] else 1+first_change(w,(i+1)%n,(j+1)%n)
   rows.add(tuple(q(i,j) for i,j in ((0,1),(0,2),(1,2))))
   if n<=12:
    for k in range(1,n):
     assert sum(q(i,(i+k)%n) for i in range(n))==n
     identity_checks+=1
  union.update(rows)
  tables[n]={'words_up_to_complement':primitive_count,'triangle_signatures':len(rows)}
 return sorted(union),tables,identity_checks

def search(b,table,N):
 fs=primes(b)
 packets=[]
 for p,e in fs:
  scale=p**(e-(p==2))
  packets.append([tuple(scale**q for q in row) for row in table])
 found=set();admissible=0
 for choice in itertools.product(*packets):
  sizes=[math.prod(row[j] for row in choice) for j in range(3)]
  # Independent sign enumeration, direct recurrence equality, no unit equation.
  for signs in itertools.product((-1,1),repeat=3):
   edge=[sizes[j]*signs[j] for j in range(3)]
   a=[edge[0]+edge[1]-edge[2],edge[0]+edge[2]-edge[1],edge[1]+edge[2]-edge[0]]
   if len(set(a))<3:continue
   C=b*a[1]-a[0]**2
   if b*a[2]-a[1]**2!=C:continue
   if any(math.gcd(t,b)!=1 for t in a) or math.gcd(C,b)!=1:continue
   admissible+=1
   x=Fraction(a[0],b);c=Fraction(C,b*b);orbit=[]
   for _ in range(N+1):
    if x in orbit:
     if x==orbit[0] and len(orbit)>=3:
      nums=[int(b*z) for z in orbit]
      can=min(tuple(nums[i:]+nums[:i]) for i in range(len(nums)))
      found.add((C,can))
     break
    if len(orbit)==N or x.denominator!=b:break
    orbit.append(x);x=x*x+c
 return sorted(found),admissible

def parametric_controls(max_b):
 found={}
 for m in range(-max_b//2,max_b//2+1):
  if m==0:continue
  for n in range(1,max_b//2+1):
   if m+n==0 or math.gcd(m,n)!=1:continue
   den=2*m*n*(m+n);b=abs(den)
   if b>max_b:continue
   A=m**6+2*m**5*n+4*m**4*n*n+8*m**3*n**3+9*m*m*n**4+4*m*n**5+n**6
   numerator=m**3+2*m*m*n+m*n*n+n**3
   a0=numerator*(1 if den>0 else -1)
   C=-A;orbit=[a0]
   for _ in range(2):
    nxt=orbit[-1]**2+C
    assert nxt%b==0
    orbit.append(nxt//b)
   assert orbit[-1]**2+C==b*orbit[0]
   can=min(tuple(orbit[i:]+orbit[:i]) for i in range(3))
   found.setdefault(b,set()).add((C,can))
 return found

def main():
 p=argparse.ArgumentParser();p.add_argument('--primary',required=True);p.add_argument('--output',required=True);p.add_argument('--skip-four-prime',action='store_true');args=p.parse_args()
 primary=json.load(open(args.primary));maxN=max(r['period_bound'] for r in primary['results']);cache={};tables={};checks=0
 for N in sorted(set(r['period_bound'] for r in primary['results'])):
  cache[N],tab,num=word_tables(N);tables[N]=tab;checks+=num
  print('word bound',N,'actual union',len(cache[N]),flush=True)
 controls=parametric_controls(max(r['b'] for r in primary['results']));audits=[]
 for row in primary['results']:
  b=row['b'];N=row['period_bound']
  expected=sorted((c['C'],tuple(c['numerators'])) for c in row['cycles'])
  assert expected==sorted(controls.get(b,set())),('parametric mismatch',b)
  if args.skip_four_prime and row['omega']>=4:
   audits.append({'b':b,'parameterization_match':True,'independent_search':'NOT RUN'});continue
  result,count=search(b,cache[N],N)
  assert expected==result,('search mismatch',b)
  audits.append({'b':b,'parameterization_match':True,'independent_search':'PASS','admissible_triangles':count,'cycles':len(result)})
  print('independent PASS b',b,'triangles',count,'cycles',len(result),flush=True)
 out={'status':'PASS','word_tables':tables,'cyclic_valuation_identity_checks':checks,'audits':audits}
 open(args.output,'w').write(json.dumps(out,indent=2)+'\n')
 print('PASS',len(audits),'denominators',checks,'cyclic valuation checks')
if __name__=='__main__':main()
