from collections import Counter
from itertools import permutations
from pathlib import Path
import json

def need(x):
 if not x: raise AssertionError

def graph(word):
 g=[0]*len(word)
 for i,m in enumerate(word):
  need(m<1<<min(i,3))
  for j in range(i):
   u,v=(i,j) if i-j<=3 and (m>>(i-j-1))&1 else (j,i)
   g[u]|=1<<v
 return g

def hcount(g):
 n=len(g);f=(1<<n)-1;v=[[0]*n for _ in range(f+1)]
 for i in range(n):v[1<<i][i]=1
 for u in range(1,f+1):
  for i,z in enumerate(v[u]):
   if not z: continue
   bits=g[i]&~u
   while bits:
    b=bits&-bits;bits^=b
    v[u|b][b.bit_length()-1]+=z
 return sum(v[-1])

def strong(g):
 n=len(g)
 for root in range(n):
  seen=1<<root
  while True:
   nxt=seen
   for i in range(n):
    if seen>>i&1:nxt|=g[i]
   if nxt==seen:break
   seen=nxt
  if seen!=(1<<n)-1:return False
 return True

def main():
 allrows=json.loads(Path(__file__).with_name('STRONG_TRANSITIONS.json').read_text())
 found=set();D=set();levels=Counter()
 for r in allrows:
  w=r['word'];a,b=graph(w[:-1]),graph(w)
  h,k=hcount(a),hcount(b)
  need(len(w)<=7 and strong(a) and strong(b))
  need((h,k)==(r['H_old'],r['H_new']))
  pair=h%16,k%16
  need(pair==(r['old_mod16'],r['new_mod16']))
  found.add(pair);D.add((h%16,(k-h)//2%8));levels[len(w)]+=1
 need(len(found)==64 and len(D)==64)
 samples=((0,0),(0,0,2),(0,0,0,4),(0,0,2,1),
    (0,0,0,4,4),(0,0,3,6,3),(0,0,0,4,0,4),(0,0,0,4,6))
 need([hcount(graph(w))*pow(hcount(graph(w[:-1])),-1,16)%16 for w in samples]==list(range(1,16,2)))
 vals=set()
 for n in range(1,7):
  for p in permutations(range(n)):
   ix={v:i for i,v in enumerate(p)}
   if any(ix[i]>ix[j] for i in range(n) for j in range(i+4,n)): continue
   old=[sum(1<<v for v in p[ix[i]+1:]) for i in range(n)]
   for a in range(1<<min(n,3)):
    g=old[:]+[0]
    for i in range(n):
     if n-i<=3 and a>>(n-i-1)&1:g[n]|=1<<i
     else:g[i]|=1<<n
    vals.add(hcount(g))
 need(vals=={1,3,5,9,11,15,17,19,23,25,27,31,33,39})
 old=(0,1,3,3,2,5)
 need(strong(graph(old)) and hcount(graph(old))==25)
 vals=[hcount(graph(old+(a,))) for a in range(8)]
 need(vals==[25,61,43,53,103,115,95,81])
 need(all(strong(graph(old+(a,))) for a in range(1,8)))
 need(not strong(graph(old+(0,))))
 deltas=[(v-25)//2 for v in vals]
 need({v%8 for v in deltas}==set(range(8)))
 for m,V in enumerate(deltas):
  a,b,c=[m>>i&1 for i in range(3)]
  need(V==18*a+9*b+39*c-13*a*b-12*a*c-13*b*c)
 two=((2,4,73),(1,6,155),(1,1,109),(1,5,143),
      (2,1,81),(4,2,275),(1,4,133),(1,2,119))
 for r,(a,b,final) in enumerate(two):
  g1,g2=graph(old+(a,)),graph(old+(a,b))
  need(strong(g1) and strong(g2) and hcount(g2)==final)
  need((final-25)//2%8==r)
 fanrows=json.loads(Path(__file__).with_name('STRONG_FANS.json').read_text())
 need(sorted(row['old_mod16'] for row in fanrows)==list(range(1,16,2)))
 for row in fanrows:
  W=tuple(row['word']);g=graph(W);h=hcount(g)
  need(strong(g) and h==row['H_old'] and h%16==row['old_mod16'])
  gs=[graph(W+(a,)) for a in range(8)];vs=list(map(hcount,gs))
  need(vs==row['ends'] and {(v-h)//2%8 for v in vs}==set(range(8)))
  need(all(strong(g) for g in gs[1:]) and not strong(gs[0]))
 print('64 exact strong witnesses and eight full fans from strong prefixes, independently verified. Transitive prefix 7/8. Orders:',dict(sorted(levels.items())))

if __name__=='__main__':main()
