#!/usr/bin/env python3
"""Portable exact weighted equality-partition census, standard library only."""
import argparse, collections, itertools, json, math, pathlib, time
PAIRS=tuple(itertools.combinations(range(9),2)); ALL=(1<<36)-1
LINES=((0,1,2),(3,4,5),(6,7,8),(0,3,6),(1,4,7),(2,5,8),(0,4,8),(2,4,6))
def require(x,s):
 if not x: raise RuntimeError(s)
def mask(v): return sum(1<<k for k,(i,j) in enumerate(PAIRS) if v[i]==v[j])
def squares(n):return {x*x%n for x in range(n)}
def histogram(n):
 q=squares(n); h=collections.Counter()
 for c,a,b in itertools.product(q,repeat=3):
  v=tuple(x%n for x in (a,3*c-a-b,b,c-a+b,c,c+a-b,2*c-b,a+b-c,2*c-a))
  if all(x in q for x in v): h[mask(v)]+=1
 return h

def independent(n):
 q=squares(n); centers=collections.defaultdict(list); h=collections.Counter()
 for e in q: centers[3*e%n].append(e)
 for a,b,c in itertools.product(q,repeat=3):
  s=(a+b+c)%n
  for e in centers[s]:
   g=(s-c-e)%n;i=(s-a-e)%n;j=(s-b-e)%n
   d=(s-a-g)%n;f=(s-c-i)%n;v=(a,b,c,d,e,f,g,j,i)
   if all(x in q for x in v) and all(sum(v[x] for x in line)%n==s for line in LINES):h[mask(v)]+=1
 return h

def meet(h,k):
 out=collections.Counter()
 for a,x in h.items():
  for b,y in k.items():out[a&b]+=x*y
 return out

def factors(n):
 out=[];p=2
 while p*p<=n:
  if n%p==0:
   q=1
   while n%p==0:n//=p;q*=p
   out.append(q)
  p+=1
 if n>1:out.append(n)
 return out

def partitions():
 def rec(v):
  if len(v)==9:
   sizes=collections.Counter(v).values()
   yield (mask(v),(-1)**(9-len(set(v)))*math.prod(math.factorial(s-1) for s in sizes));return
  for t in range(max(v)+2):yield from rec(v+(t,))
 return list(rec((0,)))
def zeta(h,partitions):return {p:sum(c for m,c in h.items() if p&m==p) for p,_ in partitions}
def rgs(m):
 out=[]
 for i in range(9):
  prior=next((j for j in range(i) if m&(1<<PAIRS.index((j,i)))),None)
  out.append(out[prior] if prior is not None else max(out,default=-1)+1)
 return ''.join(map(str,out))
def serial(h):return [{'partition_rgs':rgs(m),'collision_mask':hex(m),'count':c} for m,c in sorted(h.items())]
def build(limit):
 parts=partitions();require(len(parts)==21147,'Bell number');cache={};rows=[]
 for n in range(1,limit+1):
  h=histogram(n);cache[n]=h
  require(h==independent(n),f'independent histogram n={n}')
  crt=collections.Counter({ALL:1})
  for q in factors(n):crt=meet(crt,cache.get(q) or histogram(q))
  require(crt==h,f'CRT histogram n={n}')
  require(h[0]%8==0,f'D4 n={n}')
  rows.append({'n':n,'total_arrays':sum(h.values()),'partition_support':len(h),'distinct_labeled':h[0],'distinct_D4_orbits':h[0]//8})
 examples={};checks=[]
 for a,b in [(5,7),(5,11),(3,67),(4,9),(8,25),(7,11),(16,27)]:
  ha=cache.get(a) or histogram(a);hb=cache.get(b) or histogram(b);h=histogram(a*b)
  require(h==independent(a*b),f'independent example {a*b}')
  require(meet(ha,hb)==h,f'example CRT {a*b}')
  za=zeta(ha,parts);zb=zeta(hb,parts);z=zeta(h,parts)
  require(all(z[p]==za[p]*zb[p] for p,_ in parts),f'zeta {a*b}')
  mobius=sum(mu*z[p] for p,mu in parts);require(mobius==h[0],f'Mobius {a*b}')
  examples[str(a*b)]={'factors':[a,b],'factor_histograms':{str(a):serial(ha),str(b):serial(hb)},'full_histogram':serial(h),'total_arrays':sum(h.values()),'distinct_labeled':h[0],'distinct_D4_orbits':h[0]//8,'mobius_distinct':mobius,'zeta_partitions_checked':len(parts)}
 for n in range(1,limit//3+1):
  if n%3:
   require(cache[3*n]==collections.Counter({m:2*c for m,c in cache[n].items()}),f'3-neutrality {n}')
   checks.append(n)
 return {'schema':'weighted-fixed-cell-CRT-v1','domain':'Z/nZ; zero included; value arrays not root arrays','pair_order':PAIRS,'partition_order':'finer <= coarser; discrete minimum','bell9':len(parts),'direct_independent_and_CRT_verified_range':[1,limit],'census':rows,'examples':examples,'three_neutrality_checked_n':checks}
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--limit',type=int,default=160);ap.add_argument('--output',default='certificate.json');ap.add_argument('--verify',action='store_true');a=ap.parse_args();t=time.perf_counter();d=build(a.limit);s=json.dumps(d,sort_keys=True,indent=2)+'\n';p=pathlib.Path(a.output)
 if a.verify:require(json.loads(p.read_text())==json.loads(s),'certificate differs')
 else:p.write_text(s)
 print(json.dumps({'status':'PASS','limit':a.limit,'elapsed_seconds':round(time.perf_counter()-t,3),'examples':{n:{k:v[k] for k in ['total_arrays','distinct_labeled','distinct_D4_orbits']} for n,v in d['examples'].items()}},indent=2))
if __name__=='__main__':main()
