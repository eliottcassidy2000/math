#!/usr/bin/env python3
"""Exact Z/nZ collision-partition census and independent full-modulus check."""
import argparse,itertools,json,math,pathlib
PAIRS=list(itertools.combinations(range(9),2))
LINES=[(0,1,2),(3,4,5),(6,7,8),(0,3,6),(1,4,7),(2,5,8),(0,4,8),(2,4,6)]
ALL=(1<<36)-1

def mask(m):return sum(1<<k for k,(i,j) in enumerate(PAIRS) if m[i]==m[j])
def roots(n):return {x*x%n:x for x in range(n)}
def magic(n,c,a,b):return tuple(x%n for x in [a,3*c-a-b,b,c-a+b,c,c+a-b,2*c-b,a+b-c,2*c-a])
def valid(n,m):
 return all(x in roots(n) for x in m) and len({sum(m[i] for i in L)%n for L in LINES})==1

def enumerate_masks(n):
 Q=roots(n);D={}
 for c,a,b in itertools.product(Q,repeat=3):
  m=magic(n,c,a,b)
  if all(x in Q for x in m):D.setdefault(mask(m),m)
 return D

def independent_masks(n):
 # Enumerate TOP ROW first, using only line constraints; no magic() calls.
 Q=roots(n);centers={}
 for e in Q:centers.setdefault(3*e%n,[]).append(e)
 D={}
 for a,b,c in itertools.product(Q,repeat=3):
  s=(a+b+c)%n
  for e in centers.get(s,[]):
   g=(s-c-e)%n;i=(s-a-e)%n;h=(s-b-e)%n
   d=(s-a-g)%n;f=(s-c-i)%n
   m=(a,b,c,d,e,f,g,h,i)
   if all(x in Q for x in m) and len({sum(m[j] for j in L)%n for L in LINES})==1:D.setdefault(mask(m),m)
 return D

def minimal(D):
 out={}
 for x in sorted(D,key=lambda x:(x.bit_count(),x)):
  if not any(y&x==y for y in out):out[x]=D[x]
 return out

def factors(n):
 out=[];p=2
 while p*p<=n:
  if n%p==0:
   q=1
   while n%p==0:q*=p;n//=p
   out.append(q)
  p+=1
 if n>1:out.append(n)
 return out

def crt(v,mods):
 N=math.prod(mods)
 return sum(x*(N//q)*pow(N//q,-1,q) for x,q in zip(v,mods))%N

def compose(states,local):
 D={}
 for x,w in states.items():
  for y,m in local.items():D.setdefault(x&y,w+[m])
 return minimal(D)

def build(limit):
 cache={};rows=[]
 for n in range(2,limit+1):
  fs=factors(n)
  for q in fs:
   if q not in cache:cache[q]=enumerate_masks(q)
  states={ALL:[]}
  for q in fs:states=compose(states,minimal(cache[q]))
  direct=enumerate_masks(n);ind=independent_masks(n)
  if set(direct)!=set(ind):raise RuntimeError(f"full-mask mismatch at {n}")
  if set(states)!=set(minimal(direct)):raise RuntimeError(f"CRT mismatch at {n}")
  witness=None
  if 0 in states:
   locals_=states[0];m=tuple(crt(v,fs) for v in zip(*locals_))
   if not valid(n,m) or len(set(m))!=9:raise RuntimeError(f"invalid CRT witness at {n}")
   witness={'entries':m,'square_roots':[roots(n)[x] for x in m],'local_entries':locals_,'local_masks':[hex(mask(m)) for m in locals_]}
  rows.append({'n':n,'prime_power_components':fs,'component_parker':[0 not in cache[q] for q in fs],'parker':0 not in states,'full_masks':len(direct),'minimal_masks':[hex(x) for x in states],'witness':witness})
  if n%25==0:print('checked',n,flush=True)
 components=[{'modulus':q,'full_mask_count':len(D),'minimal':[{'mask':hex(x),'entries':m} for x,m in minimal(D).items()]} for q,D in sorted(cache.items())]
 return {'schema':'fixed-cell-crt-collision-v1','domain':'Z/nZ (not extension fields)','pairs':PAIRS,'limit_inclusive':limit,'components':components,'rows':rows,'smallest_nonparker_with_only_parker_components':next((r['n'] for r in rows if not r['parker'] and all(r['component_parker'])),None)}

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--limit',type=int,default=160);ap.add_argument('--output',default='census.json');ap.add_argument('--verify',action='store_true');a=ap.parse_args()
 data=build(a.limit)
 if a.verify:
  if json.loads(pathlib.Path(a.output).read_text())!=json.loads(json.dumps(data)):raise RuntimeError('certificate mismatch')
 else:pathlib.Path(a.output).write_text(json.dumps(data,indent=2)+'\n')
 print('PASS: two independent full-modulus enumerations and CRT antichain composition agree through',a.limit)
if __name__=='__main__':main()
