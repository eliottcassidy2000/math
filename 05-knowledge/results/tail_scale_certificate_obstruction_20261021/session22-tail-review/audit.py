#!/usr/bin/env python3
"""Independent supplied-subset audit. No producer import and no annulus sweep."""
import argparse, hashlib, itertools, json
from pathlib import Path
P=argparse.ArgumentParser();P.add_argument('--repo-root',type=Path,required=True);P.add_argument('--packet',type=Path,required=True);P.add_argument('--base-dir',type=Path,required=True);P.add_argument('--output',type=Path);P.add_argument('--verify',type=Path);args=P.parse_args()
checks=0
def ck(ok,label):
 global checks
 checks+=1
 if not ok:raise ValueError(label)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
pins={str(args.packet/'receipt.json'):'0d23c15d441b67986e7d03d64c679df9536e75c225c826bf03cd49ecf3e11cf2',str(args.packet/'check.py'):'213a6ca60c69125451f9fd483ce166ad62fca485a17eb7ca4a8b5205ac935544',str(args.packet/'README.md'):'fd5ea3c1002a6ff469c918bd330fbf91082a3c02d4179ecf1c6511c55175e22b',str(args.base_dir/'check.py'):'1d24188af100b3394b9e47b4b27e35e80b15863e9f033c05081868d7ff5a1632',str(args.base_dir/'receipt.json'):'e5ec575906451e4a6ee293a98ed85dc3d90a2a8b18b878a921df151166d9032a',str(args.repo_root/'05-knowledge/results/session18_native_annulus/frontier.json'):'391f6d3380bcae4afb9a1bbd9522b74b1e9006012fa10a92be5bf0bbf938754f'}
for p,h in pins.items():ck(sha(Path(p))==h,'source pin '+p)
r=json.loads((args.packet/'receipt.json').read_text());b=next(x for x in json.loads((args.base_dir/'receipt.json').read_text())['cases'] if x['N']==4096)
V=set(json.loads((args.repo_root/'05-knowledge/results/session18_native_annulus/frontier.json').read_text())['frontier'])
A=((1,0),(2,1),(1,1),(1,2),(0,1),(-1,1))
def q(v):return v[0]**2-v[0]*v[1]+v[1]**2
def beta(v,w):return q((v[0]+w[0],v[1]+w[1]))-q(v)-q(w)
def primitive(s):return s>0 and s%8 in (1,3,7)
def step(s):
 ck(primitive(s),'primitive source')
 u=3*s+1;t=0
 while u%2==0:u//=2;t+=1
 while u%8==5:u=(u-1)//4;t+=2
 ck(primitive(u),'primitive endpoint')
 return u,t
def charge(scale,s):
 x=y=0
 while s<4*scale:
  v=(((1,2),(-2,-1),(1,-1)) if s<2*scale else ((0,1),(-1,-1),(1,0)))[s%3]
  x+=v[0];y+=v[1];s=4*s+1
 return x,y
def energy(c):
 vals=[]
 for signs in itertools.product((-1,1),repeat=6):
  v=tuple(sum(signs[i]*c[i]*A[i][j] for i in range(6)) for j in (0,1));vals.append((q(v),v,signs))
 return max(vals)
ck(r['N']==262144 and r['horizon']==111,'literal scale cutoff');N=r['N'];h=r['horizon'];counts=[0]*12;axes=[0]*6;seen=set();digest=hashlib.sha256();clocks=[]
for row in r['pending_rows']:
 s,g15,g35,x,y,end,total=row
 ck(s not in seen,'distinct witness');seen.add(s)
 ck(128*N<=s<256*N and primitive(s),'upper source band')
 ck((g15,g35)==(s%15,s%35) and g15%5==g35%5,'joint guards')
 ck(charge(64*N,s)==(x,y) and q((x,y))==1,'independent charge')
 v=s;clock=0
 for j in range(h):
  ck(v not in V,'no early frontier membership');w,t=step(v)
  digest.update(f'{s},{j},{v},{v%15},{v%35},{w},{t}\n'.encode());v=w;clock+=t
 ck(v not in V and v==end and clock==total,'literal cutoff state and clock')
 clocks.append(clock)
 full=A+tuple((-a,-b) for a,b in A);i=full.index((x,y));counts[i]+=1;axes[i%6]+=1
ck(len(seen)==4096,'subset size');D,vertex,signs=energy(axes)
baseaxes=[b['pending_new_counts'][i]+b['pending_new_counts'][i+6] for i in range(6)]
ck(b['old_unarrived']==0 and b['horizon']==105,'inherited old-complete cut')
ck(energy(baseaxes)[0]==b['D_new']==91876,'independent base energy')
ck(counts==r['pending_direction_counts'] and axes==r['pending_subset_axis_totals'],'reported counts')
ck(D==r['pending_subset_D']==7499929 and D>64*b['D_new'],'strict D hostile')
ck(len(seen)>8*b['pending_nonzero'] and any(a>8*c for a,c in zip(axes,baseaxes)),'count hostile')
# Lift formula and exact periodic inverse classes, plus integer clock check.
for s in range(1,4097,2):
 if not primitive(s):continue
 u,t=step(s);a=2 if s%8==1 else 1;L=64*s+7*(2**a+3)
 ck(9*s==2**t*(3*u+1)-(2**a+3),'inverse edge identity')
 ck(step(L)==(u,t+6),'clock lift identity')
 ck(L%512=={1:113,3:227,7:483}[s%8],'lift residue')
for z in range(1,4097):
 if z%512 in (113,227,483):
  c=49 if z%512==113 else 35;s=(z-c)//64
  ck(z-c==64*s and primitive(s),'positive inverse residue image')
x,y=32881,247921
ck(step(x)==(385,8) and step(y)==(2905,8),'hostile actual edges')
ck(x not in V and y not in V and 385 in V and 2905 not in V,'hostile membership')
ck(charge(64*4096,x)==charge(64*4096,y)==(-1,-1),'hostile charges')
ck((x%15,x%35)==(y%15,y%35)==(1,16),'hostile guards')
for s,length,T,shift in r['translation_controls']:
 path=[s];ts=[]
 for j in range(length):
  u,t=step(path[-1]);path.append(u);ts.append(t)
 ck(T==sum(ts) and shift==840*2**T,'cylinder input')
 for k in (0,1,2,7):
  u=s+k*shift;cum=0
  for j,t in enumerate(ts):
   ck(u==path[j]+840*k*3**j*2**(T-cum) and u%105==path[j]%105,'cylinder state/guards')
   u,tt=step(u);ck(tt==t,'cylinder exact integer clock');cum+=t
  ck(u==path[-1]+840*k*3**length,'cylinder endpoint')
out={'status':'PASS: independent subset replay and scoped proof review','input_sha256':sorted(pins.values()),'subset_size':len(seen),'literal_edges':len(seen)*h,'prefix_rows_sha256':digest.hexdigest(),'base_axes':baseaxes,'subset_axes':axes,'base_D':b['D_new'],'subset_D':D,'maximizing_axis_vertex':vertex,'maximizing_axis_signs':signs,'candidate_bound':64*b['D_new'],'strict_gap':D-64*b['D_new'],'clock_range':[min(clocks),max(clocks)],'checks':checks,'scope':'Only supplied witnesses; no full enlarged annulus census, no claim on 153N failure or Collatz entry.'}
blob=(json.dumps(out,sort_keys=True,separators=(',',':'))+'\n').encode()
if args.output:args.output.write_bytes(blob)
if args.verify:ck(args.verify.read_bytes()==blob,'immutable replay')
print(blob.decode().strip());print('receipt_sha256='+hashlib.sha256(blob).hexdigest())
