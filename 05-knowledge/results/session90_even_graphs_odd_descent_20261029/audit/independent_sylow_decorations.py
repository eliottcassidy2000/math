#!/usr/bin/env python3
"""Independent n=5 Hall audit: literal S5 actions, not Burnside formulas."""
from independent_decorated_counts import setup, check
from collections import Counter
import json

def compose(p,q):return tuple(p[q[i]] for i in range(len(p)))
def inverse(p):return tuple(p.index(i) for i in range(len(p)))
def generated(g):
 e=tuple(range(len(g)));a=e;out=[]
 while a not in out:out.append(a);a=compose(g,a)
 return set(out)

def census(n,p):
 es,acts=setup(n);m=len(es);q=1<<m
 def pack(x,r,d):return x|(r<<m)|(d<<(2*m))
 gen=tuple(list(range(1,p))+[0]+list(range(p,n)));P=generated(gen)
 normal=[a for a in acts if {compose(compose(a[0],g),inverse(a[0])) for g in P}==P]
 gact=next(a for a in acts if a[0]==gen)
 tfix=[x for x in range(q) if gact[1][x]^gact[3]==x]
 xfix=[x for x in range(q) if gact[1][x]==x]
 dfix=[d for d in range(1<<n) if gact[2][d]==d]
 T=set();G=set();normalT=set();normalG=set();rawEven=0
 for x in tfix:
  for r in xfix:
   for d in dfix:
    T.add(min(pack(a[1][x]^a[3],a[1][r],a[2][d]) for a in acts))
    normalT.add(min(pack(a[1][x]^a[3],a[1][r],a[2][d]) for a in normal))
 for x in xfix:
  for r in xfix:
   for d in dfix:
    def even(group):
     return all((x&a[4]).bit_count()%2==0 for a in group if a[1][x]==x and a[1][r]==r and a[2][d]==d)
    e=even(acts);ne=even(normal)
    check(e==ne,'Frattini normalizer/global sign')
    if e:
     rawEven+=1
     G.add(min(pack(a[1][x],a[1][r],a[2][d]) for a in acts))
     normalG.add(min(pack(a[1][x],a[1][r],a[2][d]) for a in normal))
 check(len(T)==len(G)==len(normalT)==len(normalG),'Sylow class equality')
 return T,G,dict(prime=p,normalizer=len(normal),fixed_pairs=len(tfix)*len(xfix)*len(dfix),fixed_even_graphs=rawEven,classes=len(T))

def main():
 t3,g3,r3=census(5,3);t5,g5,r5=census(5,5)
 check(not(t3&t5),'odd source cannot contain both Sylows')
 check(bool(g3&g5),'target intersection nonempty')
 check(len(t3|t5)>len(g3|g5),'strict Hall obstruction')
 import independent_decorated_counts as lib
 print(json.dumps(dict(status='PASS',checks=lib.CHECKS,rows=[r3,r5],source_union=len(t3|t5),target_union=len(g3|g5),target_intersection=len(g3&g5)),indent=2))
if __name__=='__main__':main()
