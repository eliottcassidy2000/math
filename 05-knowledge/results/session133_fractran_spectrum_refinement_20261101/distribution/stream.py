"""Positive linear-map to literal ordered FRACTRAN; no compiler imports."""
from fractions import Fraction
from pathlib import Path
import json, math

ATOMS=[(3,3,1),(4,5,1),(5,9,2),(5,11,1),(5,13,1),(5,15,2),(6,15,2),(6,17,2),(6,19,2)]
H=list(range(1,22,2)); L=6

def primes(n):
 out=[];p=2
 while len(out)<n:
  if all(p%q for q in out if q*q<=p):out.append(p)
  p+=1
 return out

def matrix():
 d=len(H)*L; a=[[0]*d for _ in range(d)]
 ix=lambda h,k:H.index(h)*L+k
 for h in H:
  a[ix(h,0)][ix(h,0)]+=1
  for k,q,c in ATOMS:
   if h%q==0:a[ix(h,0)][ix(h//q,k-1)]+=c
  for k in range(1,L):a[ix(h,k)][ix(h,k-1)]=1
 return a

def compile_linear(a):
 d=len(a)
 if not d or any(len(r)!=d or any(type(v)!=int or v<0 for v in r) for r in a):raise ValueError('nonnegative square integer matrix required')
 ps=primes(6*d+2);x=ps[:d];y=ps[d:2*d];c=ps[2*d:6*d];Q,clock=ps[-2:]
 # Two drain phases per register, with an explicit return trampoline.
 rows=[]; factors=[]
 def add(num,den):
  nv=math.prod(p**e for p,e in num.items());dv=math.prod(p**e for p,e in den.items())
  f=Fraction(nv,dv); rows.append((f.numerator,f.denominator));factors.append((num,den))
 add({c[0]:1},{Q:1})
 for i in range(d):
  C,D=c[2*i:2*i+2];nxt=c[2*i+2]
  num={D:1};num.update({y[j]:a[j][i] for j in range(d) if a[j][i]})
  add(num,{C:1,x[i]:1});add({nxt:1},{C:1});add({C:1},{D:1})
 for i in range(d):
  C,D=c[2*d+2*i:2*d+2*i+2]
  add({D:1,x[i]:1},{C:1,y[i]:1})
  add(({c[2*d+2*i+2]:1} if i<d-1 else {Q:1,clock:1}),{C:1})
  add({C:1},{D:1})
 return {'dimension':d,'x':x,'y':y,'control':c,'checkpoint':Q,'clock':clock,'fractions':rows},factors

def run(a,v,rounds,literal=False):
 if len(v)!=len(a) or any(type(e)!=int or e<0 for e in v):raise ValueError('nonnegative vector of matrix dimension required')
 if type(rounds)!=int or rounds<0:raise ValueError('nonnegative integer rounds required')
 p,fs=compile_linear(a);ps=sorted(set(p['x']+p['y']+p['control']+[p['checkpoint'],p['clock']]))
 ex={q:0 for q in ps};ex[p['checkpoint']]=1
 for q,e in zip(p['x'],v):ex[q]=e
 N=math.prod(q**e for q,e in ex.items()) if literal else None
 result=[list(v)];steps=0;counts=[]
 # Index control only: every row has exactly one control in denominator.
 controls=set(p['control']+[p['checkpoint']]);index={q:[] for q in controls}
 for k,(_,den) in enumerate(fs):
  cs=set(den)&controls
  if len(cs)!=1:raise RuntimeError('bad control')
  index[next(iter(cs))].append(k)
 current=p['checkpoint'];prev=0
 while len(result)<=rounds:
  candidates=range(len(fs)) if literal else index[current]
  for k in candidates:
   num,den=fs[k]
   eligible=N%p['fractions'][k][1]==0 if literal else all(ex[q]>=e for q,e in den.items())
   if eligible:
    if literal:N=N//p['fractions'][k][1]*p['fractions'][k][0]
    for q,e in den.items():ex[q]-=e
    for q,e in num.items():ex[q]+=e
    current=next(q for q in num if q in controls)
    steps+=1;break
  else:raise RuntimeError('unexpected halt')
  if current==p['checkpoint']:
   out=[ex[q] for q in p['x']]
   if any(ex[q] for q in p['y']) or sum(ex[q] for q in controls)!=1:raise RuntimeError('bad checkpoint')
   if literal and N!=math.prod(q**e for q,e in ex.items()):raise RuntimeError('integer mismatch')
   result.append(out);counts.append(steps-prev);prev=steps
 return result,counts,p

def ordinary_table(nmax):
 out=[]
 for n in range(nmax+1):
  row={h:0 for h in H}
  if n==0:row[1]=1
  else:
   for h in H:
    row[h]=out[n-1][h]
    for k,q,c in ATOMS:
     if n>=k and h%q==0:row[h]+=c*out[n-k][h//q]
  out.append(row)
 return out

if __name__=='__main__':
 a=matrix();p,_=compile_linear(a);base=Path(__file__).parent
 (base/'bounded21.fractran').write_text('\n'.join(f'{n}/{d}' for n,d in p['fractions'])+'\n')
 (base/'program.json').write_text(json.dumps({k:v for k,v in p.items() if k!='fractions'},indent=2)+'\n')
 print(len(p['fractions']),p['dimension'],p['checkpoint'],p['clock'])

def tail_matrix():
 """All eleven N_n(h) as state, for n>=6; no lag registers."""
 d=len(H);a=[[int(i==j) for j in range(d)] for i in range(d)]
 for h,c in [(3,1),(5,1),(11,1),(13,1),(17,2),(19,2)]:a[H.index(h)][0]+=c
 a[H.index(9)][H.index(3)]+=1
 a[H.index(15)][H.index(5)]+=2
 return a

def write_tail():
 a=tail_matrix();p,_=compile_linear(a);v=[ordinary_table(6)[6][h] for h in H]
 seed=p['checkpoint']*p['clock']**6*math.prod(q**e for q,e in zip(p['x'],v));base=Path(__file__).parent
 (base/'tail21.fractran').write_text('\n'.join(f'{n}/{d}' for n,d in p['fractions'])+'\n')
 data={k:v for k,v in p.items() if k!='fractions'};data.update(initial_order=6,initial_vector=v,initial_integer=str(seed))
 (base/'tail_program.json').write_text(json.dumps(data,indent=2)+'\n')
 return data

if __name__=='__main__':write_tail()
