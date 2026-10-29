"""Literal cyclic words, golden CRT readers, and phase-compatible mirror charts."""
from pathlib import Path
from collections import Counter,defaultdict
from math import gcd,lcm,prod,isqrt
import json
HERE=Path(__file__).resolve().parent
CHECKS=0
def need(v,m):
 global CHECKS
 CHECKS+=1
 if not v:raise RuntimeError(m)
def fib(n):
 a,b=0,1
 for _ in range(n):a,b=b,a+b
 return a
def lucas(n):return fib(n-1)+fib(n+1)
def factors(n):
 out=[];p=2
 while p*p<=n:
  if n%p==0:
   e=0
   while n%p==0:n//=p;e+=1
   out.append((p,e))
  p+=1
 if n>1:out.append((n,1))
 return out
def words(n):
 def rec(i,w,previous,first):
  if i==n:
   if not(previous and first):yield w
   return
  yield from rec(i+1,w,0,first)
  if not previous:yield from rec(i+1,w|(1<<i),1,first if i else 1)
 yield from rec(0,0,0,0)
def rotate(w,n):return ((w<<1)&((1<<n)-1))|(w>>(n-1))
def reflect(w,n):return sum(((w>>i)&1)<<((-i)%n) for i in range(n))
def golden_mul(x,y,m):
 a,b=x;c,d=y
 return ((a*c+b*d)%m,(a*d+b*c+b*d)%m)
def goldpow(n,m):
 x=(1,0)
 for _ in range(n):x=golden_mul(x,(0,1),m)
 return x
def prime_order(x,p):
 a=1
 for k in range(1,p):
  a=a*x%p
  if a==1:return k
 raise RuntimeError('unit order')

def analyze(n):
 bank=list(words(n));need(len(bank)==lucas(n),'literal Lucas word count')
 if n%2:
  need(n%3==0,'declared odd CRT family')
  fs=factors(lucas(n)//4);need(all(e==1 and p%2 for p,e in fs),'squarefree odd norm factors')
  A=(fib(n-1)-1)//2;B=fib(n)//2
  roots=[]
  for p,e in fs:
   need(B%p!=0,'linear local ideal')
   r=-A*pow(B,-1,p)%p
   need((r*r-r-1)%p==0 and pow(r,n,p)==1,'local golden root')
   roots.append((p,r,prime_order(r,p)))
  periods=[3]+[d for p,r,d in roots]
  sizes=[4]+[p for p,r,d in roots]
  powers=[[(1,2,3)[i%3]]+[pow(r,i,p) for p,r,d in roots] for i in range(n)]
  def read(w):
   z=[0]*len(sizes)
   for i in range(n):
    if w>>i&1:
     z[0]^=powers[i][0]
     for j in range(1,len(z)):z[j]=(z[j]+powers[i][j])%sizes[j]
   return tuple(z)
  def turn(z):return ((0,2,3,1)[z[0]],)+tuple(z[j+1]*r%p for j,(p,r,d) in enumerate(roots))
  algebra={'type':'F4 times split prime fields','roots':roots,'local_orders':periods,'sizes':sizes}
 else:
  need(n==12,'declared even boundary')
  powers=[goldpow(i,8)+(pow(3,i,5),) for i in range(n)]
  def read(w):return tuple(sum(powers[i][j] for i in range(n) if w>>i&1)%m for j,m in enumerate((8,8,5)))
  def turn(z):return golden_mul(z[:2],(0,1),8)+(3*z[2]%5,)
  algebra={'type':'O/(8) times F5','ring_size':320,'phi_mod5':3}
 readers={w:read(w) for w in bank};fibres=defaultdict(list)
 for w,z in readers.items():
  fibres[z].append(w);need(readers[rotate(w,n)]==turn(z),'reader preserves rotation')
 need(len(fibres)==(lucas(n) if n%2 else 320),'reader image size')
 need(Counter(map(len,fibres.values()))==({1:lucas(n)} if n%2 else {1:319,3:1}),'reader fibre law')
 if n%2:
  need(prod(sizes)==lucas(n),'CRT cardinality')
  for z in readers.values():
   period=lcm(*(d for d,x in zip(periods,z) if x)) if any(z) else 1
   cur=z
   for _ in range(period):cur=turn(cur)
   need(cur==z,'compatible CRT clock returns')
 chart={};seeds=[]
 for w in sorted(bank):
  if w in chart:continue
  cyc=[];v=w
  while v not in chart:
   chart[v]=(w,len(cyc));cyc.append(v);v=rotate(v,n)
  need(v==w,'rotation orbit closes')
  seeds.append((w,len(cyc)))
 lengths=dict(seeds);mirror=[]
 for w,L in seeds:
  target,j=chart[reflect(w,n)]
  need(lengths[target]==L,'reflection preserves period')
  v=w
  for i in range(L):
   need(chart[reflect(v,n)]==(target,(j-i)%L),'chart mirror with retained phase')
   v=rotate(v,n)
  mirror.append([w,L,target,j])
 ambiguity=[]
 projections=([tuple(range(1,len(next(iter(fibres)))))]+[(j,) for j in range(len(next(iter(fibres))))]
              if n%2 else [(0,1),(2,),(0,),(1,)])
 for indices in projections:
  observations=defaultdict(set);witnesses=defaultdict(dict)
  for w,z in readers.items():
   key=tuple(z[j] for j in indices);out=tuple(readers[reflect(w,n)][j] for j in indices)
   observations[key].add(out);witnesses[key].setdefault(out,w)
  largest=max(observations,key=lambda key:len(observations[key]))
  ambiguity.append({'coordinates':indices,'output_count_histogram':sorted(Counter(map(len,observations.values())).items()),
                    'max_outputs':len(observations[largest]),'witness_input':largest,
                    'witnesses':[[w,out] for out,w in sorted(witnesses[largest].items())]})
 out={'period':n,'word_count':len(bank),'algebra':algebra,'reader_fibre_sizes':sorted(Counter(map(len,fibres.values())).items()),
      'chart_period_histogram':sorted(Counter(L for w,L in seeds).items()),'chart_count':len(seeds),
      'reflection_fixed_words':sum(reflect(w,n)==w for w in bank),'projection_ambiguities':ambiguity}
 if n%2:
  need(out['reflection_fixed_words']==fib((n+3)//2),'odd reflection count')
  units=sum(all(x for x in readers[w]) for w,L in seeds)
  expected=prod(m-1 for m in sizes)//lcm(*periods)
  need(units==expected,'unit chart count from compatible local phases')
  out.update(unit_charts=units,unit_relative_phase_classes=prod(periods)//lcm(*periods),
             coarse_unit_local_charts=prod((m-1)//d for m,d in zip(sizes,periods)))
 else:out['zero_words']=fibres[(0,0,0)]
 # Store one compact table at the new first odd composite, not every state.
 if n==15:out['mirror_chart_seeds']=mirror
 return out

def main(write=True):
 global CHECKS
 CHECKS=0
 rows=[analyze(n) for n in (9,12,15,21,27)]
 out={'status':'FINITE-EXACT literal words; stated CRT models only','rows':rows,'checks':CHECKS}
 if write:(HERE/'crt.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8')
 print(json.dumps({'checks':CHECKS,'summary':[{k:r[k] for k in ('period','word_count','chart_count','reflection_fixed_words')} for r in rows]}))
 return out
if __name__=='__main__':main()
