#!/usr/bin/env python3
import json, pathlib, argparse
P=pathlib.Path(__file__).resolve().parent
LINES=((0,1,2),(3,4,5),(6,7,8),(0,3,6),(1,4,7),(2,5,8),(0,4,8),(2,4,6))
def ck(x,s):
 if not x: raise RuntimeError(s)
def grid(c,a,b):return (a,3*c-a-b,b,c-a+b,c,c+a-b,2*c-b,a+b-c,2*c-a)
def rootsmod(x,m):return [r for r in range(m) if r*r%m==x%m]
def liftroot(r,x,p,k):
 m=p**k;choices=[r+t*m for t in range(p) if (r+t*m)**2%(p*m)==x%(p*m)]
 ck(len(choices)==1,'odd unit lift unique');return choices[0]
def census(m):
 Q={r*r%m for r in range(m)};count=0;first=None;patterns=set()
 for c in Q:
  for a in Q:
   for b in Q:
    v=tuple(z%m for z in grid(c,a,b))
    if not all(z in Q for z in v):continue
    if m in (3,4,8):patterns.add(tuple(z for z in v))
    if len(set(v))==9:
     count+=1
     if first is None:first=v
 return {'modulus':m,'square_residues':len(Q),'distinct_magic_arrays':count,'first':first,'all_magic_square_residue_arrays':sorted(patterns) if m in(3,4,8) else None}
def run():
 rows=[census(m) for m in (2,3,4,8,9,16,27,32,64,128)]
 ck(all(r['distinct_magic_arrays']==0 for r in rows if r['modulus'] not in(27,128)),'negative census')
 ck(all(r['distinct_magic_arrays']>0 for r in rows if r['modulus'] in(27,128)),'positive census')
 L=(8,1,6,3,5,7,4,9,2);families=[]
 for p,k0,scale,shift in [(3,3,3,-2),(2,7,8,-7)]:
  v=tuple(scale*t+shift for t in L);sums=tuple(sum(v[j] for j in l) for l in LINES)
  ck(len(set(sums))==1,'integer magic template')
  levels=[]
  for k in range(k0,k0+5):
   m=p**k; rr=[rootsmod(x,m) for x in v]
   ck(all(rr) and len(set(x%m for x in v))==9,'tower witness')
   levels.append({'k':k,'modulus':m,'least_roots':[r[0] for r in rr],'no_wrap':sums[0]<=m})
  families.append({'p':p,'entries':v,'sum':sums[0],'levels':levels})
 # Smooth p=59 first Hensel layer: p^3 parameter choices, unique lifts of specified unit roots.
 # Verify 101 deterministic parameter triples and several 3-adic layers; proof gives full count.
 p=59;base=(28,3,17,5,16,27,15,29,4);r0=[rootsmod(x,p)[0] for x in base];samples=0
 for t in range(101):
  c=16+p*(t%p);a=28+p*((3*t+1)%p);b=17+p*((t*t+2)%p)
  v=grid(c,a,b);rr=[liftroot(r,x,p,1) for r,x in zip(r0,v)]
  ck(all((r*r-x)%(p*p)==0 for r,x in zip(rr,v)),'lift correctness')
  ck(len({sum(rr[j]*rr[j] for j in l)%(p*p) for l in LINES})==1,'lift magic')
  samples+=1
 fixed=[]
 for x,r in zip(base,r0):
  chain=[r]
  for k in range(1,9):r=liftroot(r,x,p,k);chain.append(r)
  fixed.append(chain)
 return {'censuses':rows,'lo_shu_towers':families,'odd_smooth_samples':samples,'fixed_p59_root_towers':fixed}
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--write',action='store_true');a=ap.parse_args();out=json.dumps(run(),indent=2,sort_keys=True)+'\n'
 if a.write:(P/'certificate.json').write_text(out)
 else:ck((P/'certificate.json').read_text()==out,'certificate match')
 print('PASS: prime-power census, integer templates, compatible unit-root towers')
if __name__=='__main__':main()
