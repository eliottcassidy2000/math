#!/usr/bin/env python3
"""600 independent code channels through the two protected Q49 interfaces."""
import hashlib,json
from math import isqrt
from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent;RESULTS=HERE.parent
sys.path.insert(0,str(RESULTS/'session124_mixed_radix_square_rays_20261031'))
from core import Constructor,prefix,threshold
OLD=RESULTS/'seeded_eventual_cycle_theorem_20261019'
NATIVE=RESULTS/'square_sum_spectrum_overnight_20261031a'
P=(1,8,17,19,6,10,15,21,4,12,24)
checks=0
def need(x,m):
 global checks
 checks+=1
 if not x:raise RuntimeError(m)
def edge(x,y):return(min(x,y),max(x,y))

def assignments(alphabet):
 # Deterministic unrelated choices; every one of 600 source addresses retained.
 return {(a,b):alphabet[int.from_bytes(hashlib.sha256(f'{a}:{b}'.encode()).digest()[:4],'big')%len(alphabet)]
         for a in range(1,25) for b in range(25)}

def decorated_A(source_pairs,codes,outer,family,n):
 m=len(next(iter(source_pairs.values()))[0]);d=n-49*m
 record=json.loads((OLD/'certified-pair-transitions'/f'{family}_p{m%2}_d{d}.json').read_text())
 result=[]
 for word in record['words']:
  out=[]
  for atom in word:
   kind,c=atom[:2]
   if kind=='single':out.append(c)
   elif kind=='head':out.append(49-c)
   else:
    _,c,i,o=atom
    code=codes.get((abs(c),outer),'')
    p=source_pairs[code][i]
    block=[49*x+(-1)**(j+1)*c for j,x in enumerate(p)]
    if kind=='tail':block=block[1:]
    out.extend(block if o==1 else block[::-1])
  result.append(out)
 colors={x:j%2 for j,x in enumerate(result[0])}
 for i,out in enumerate(result):
  need(len(out)==n+i and sorted(out)==list(range(1,n+i+1)),'decorated A ownership')
  need(out[:2]==[69,127],'fixed A prefix')
  need(abs(out.index(2)-out.index(7))==1 and out.index(2)%2==0,'fixed A cut ports')
  for x,y in zip(out,out[1:]):need(isqrt(x+y)**2==x+y,'decorated A edge')
  if i:need(all(colors[x]==j%2 for j,x in enumerate(out) if x<=n),'decorated A common parity')
 return result

def audit(N,alphabet):
 codes=assignments(alphabet);minimum=max(threshold(w) for w in codes.values())
 need(N>=2401*minimum+1200,'joint coded threshold')
 n=(N-24)//49;d=N-49*n;m=(n-24)//49
 ctor=Constructor();pairs={w:ctor.pair(m,w) for w in sorted(set(codes.values())|{''})}
 for w,pair in pairs.items():need(all(p[:len(prefix(w))]==prefix(w) for p in pair),'input code prefix')
 family='Aplus' if (n+int(d>=57))%2==0 else 'Aminus'
 banks=[decorated_A(pairs,codes,b,family,n) for b in range(25)]
 orders=[];cuts=[]
 for pair in banks:
  orders.append(tuple(1 if p.index(2)<p.index(7) else -1 for p in pair))
  cuts.append(tuple(min(p.index(2),p.index(7))+1 for p in pair))
 need(len(set(orders))==1 and len(set(cuts))==1,'all final cut positions and orientations retained')
 table=json.loads((NATIVE/'n24_forest_table.json').read_text())
 record=next(r for r in table if (r['parity'],r['residue'],r['order0'],r['order1'])==(n%2,d,*orders[0]))
 required={edge(a,b) for a,b in zip(P,P[1:])}|{edge(482,3362),edge(3362,6242)}
 used=set();ray_count=0;length_hist={}
 for u in range(-24,25):
  if u in (0,2):continue
  for v in range(-24,25):
   w=codes[(abs(u),abs(v))];t=49*u+v
   ray=[2401*x+(-1)**j*t for j,x in enumerate(prefix(w))]
   need(not used.intersection(ray),'independent rays are disjoint')
   used.update(ray);required.update(edge(a,b) for a,b in zip(ray,ray[1:]))
   ray_count+=1;length_hist[len(ray)]=length_hist.get(len(ray),0)+1
 need(ray_count==2303 and not used.intersection(P),'ray/native separation')
 missing=set(required);seen=bytearray(N+1);first=last=None;count=0;digest=hashlib.sha256()
 def vertices():
  for atom in record['word']:
   if atom[0]=='single':yield atom[1]
   else:
    _,c,i,part,o=atom;p=banks[abs(c)][i];cut=cuts[abs(c)][i]
    lo,hi=(0,cut) if part==0 else(cut,len(p))
    for k in(range(lo,hi) if o==1 else range(hi-1,lo-1,-1)):
     yield 49*p[k]+(-1)**(k+1)*c
 for x in vertices():
  need(1<=x<=N and not seen[x],'final interval ownership');seen[x]=1;count+=1
  digest.update((str(x)+'\n').encode())
  if first is None:first=x
  if last is not None:
   need(isqrt(last+x)**2==last+x,'final square edge');missing.discard(edge(last,x))
  last=x
 need(count==N and sum(seen)==N,'complete interval')
 need(isqrt(first+last)**2==first+last,'closing edge');missing.discard(edge(first,last))
 need(not missing,'every prescribed edge retained')
 return {'N':N,'alphabet':alphabet,'source_index':m,'independent_channels':600,'rays':ray_count,
         'length_histogram':length_hist,'prescribed_edges':len(required),'threshold':2401*minimum+1200,
         'assignments_sha256':hashlib.sha256(json.dumps(sorted((a,b,w) for (a,b),w in codes.items())).encode()).hexdigest(),
         'cycle_decimal_lines_sha256':digest.hexdigest(),'constructor_checks':ctor.checks,'base_checks':ctor.bases.checks}

def hostile():
 ctor=Constructor();a=ctor.pair(1062,'A')[0];b=ctor.pair(1062,'B')[0]
 pa={x:i%2 for i,x in enumerate(a)};pb={x:i%2 for i,x in enumerate(b)}
 bad=sorted(x for x in pa if pa[x]!=pb[x]);need(bool(bad),'different phase hostile exists')
 x=bad[0]
 # Zero-based signs: label =49x-(-1)^position*c.
 left=49*x-(-1)**pa[x];right=49*x+(-1)**pb[x]
 need(left==right,'unpaired offset collision')
 return {'source_index':1062,'codes':['A','B'],'parity_disagreements':len(bad),
         'first_source_label':x,'zero_based_parities':[pa[x],pb[x]],'offsets':[1,-1],
         'duplicated_lifted_label':left,'missing_band_partner':98*x-left}

def main():
 h=hostile()
 controls=[audit(2551062,'AB'),audit(2551063,'AB'),audit(5000082,'ABC')]
 print(json.dumps({'status':'passed','checks':checks,'hostile':h,'controls':controls,
                   'scope':'Uniform 600-channel theorem uses opposite-offset parity agreement and the inherited complete interfaces; these controls expand three complete cycles.'},indent=2,sort_keys=True))

if __name__=='__main__':main()
