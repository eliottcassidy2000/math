#!/usr/bin/env python3
"""Universal integer replay of whole-block or arbitrary single-prefix-split transitions."""
import json,sys
from math import isqrt
from pathlib import Path
checks=0
def check(ok,msg):
 global checks
 checks+=1
 if not ok:raise RuntimeError(msg)
p=Path(sys.argv[1]);a=json.loads(p.read_text());check(a['status']=='witness','witness status')
par=a['parity'];d0=a['residue'];family=a['family']
check(par in (0,1) and 24<=d0<=72 and family in ('B','Aplus','Aminus'),'state')
split=a.get('split_prefix_blocks',[])
check(isinstance(split,list) and len(split)<=1 and all(isinstance(c,int) and -24<=c<=24 for c in split) and len(a['words'])==2,'split declaration and word count')
check(family=='B' or a.get('cut27') is True,'A transition must protect2–7')
profiles=[]
for d,word in zip((d0,d0+1),a['words']):
 seen_c=set();seen_s=set();heads=0;profile={};position=1;last=None
 for atom in word:
  kind=atom[0]
  if kind=='single':
   _,v=atom;check(1<=v<=24 and v not in seen_s,'singleton coverage');seen_s.add(v)
   profile[('s',v)]=position;left=right=v;length=1
  elif kind=='head':
   check(len(atom)==2 and atom[1] in split,'split head');heads+=1;profile[('h',atom[1])]=position
   left=right=49-atom[1];length=1
  else:
   check(kind in ('block','tail'),'atom type')
   _,c,i,o=atom;check(-24<=c<=24 and c not in seen_c,'block coverage');seen_c.add(c)
   check((kind=='tail')==(c in split),'only declared block is split')
   check(i==int(49+(-1)**(par+1)*c<=d),'block selection');check(o in (-1,1),'orientation')
   ell=(par+i)%2;left=147+c if c in split else 49-c;right=49*(8 if ell==0 else 10)+(-1)**ell*c
   original_start_parity=0 if c in split else 1
   if o==1:shift=(position-original_start_parity)%2
   else:left,right=right,left;shift=(position-ell)%2
   profile[('b',c)]=shift;length=(ell-(c in split))%2
  if last is not None:check(isqrt(last+left)**2==last+left,'square connector')
  last=right;position=(position+length)%2
 check(seen_c==set(range(-24,25)) and seen_s==set(range(1,25)) and heads==len(split),'exact partition')
 if split:check(profile[('h',split[0])]==(1+profile[('b',split[0])])%2,'common head-tail shift')
 check(position==(1+par+d)%2,'expanded length parity')
 check(word[:2]==[['single',1],['single',3]] if family=='B' else word[0]==['block',-20,int(49+(-1)**(par+1)*(-20)<=d),1],'required prefix')
 terminal=(10 if (par+d)%2==0 else 8) if family=='Aplus' else (8 if (par+d)%2==0 else 10)
 check(word[-1]==['single',terminal],'required endpoint')
 if a.get('cut27'):
  check(any(word[k:k+2] in ([['single',2],['single',7]],[['single',7],['single',2]]) for k in range(len(word)-1)),'protected edge')
  check(profile[('s',2)]==1 and profile[('s',7)]==0,'protected parity')
 profiles.append(profile)
check(profiles[0]==profiles[1],'common parity profiles')
base=json.loads((Path(__file__).parent/'low_prefix_pair_test.json').read_text())
if base['n']%2==par:
 expanded=[]
 for d,word in zip((d0,d0+1),a['words']):
  out=[]
  for atom in word:
   if atom[0]=='single':out.append(atom[1])
   elif atom[0]=='head':out.append(49-atom[1])
   else:
    kind,c,i,o=atom;block=[49*v+(-1)**(j+1)*c for j,v in enumerate(base['paths'][i])]
    if kind=='tail':block=block[1:]
    out.extend(block if o==1 else block[::-1])
  check(sorted(out)==list(range(1,49*base['n']+d+1)),'expanded coverage')
  check(all(isqrt(x+y)**2==x+y for x,y in zip(out,out[1:])),'expanded square edges')
  check(out[:2]==([1,3] if family=='B' else [69,127]),'expanded prefix')
  terminal=(10 if len(out)%2==0 else 8) if family=='Aplus' else (8 if len(out)%2==0 else 10)
  check(out[-1]==terminal,'expanded endpoint')
  expanded.append(out)
 profile={v:(j+1)%2 for j,v in enumerate(expanded[0])}
 check(all(profile[v]==(j+1)%2 for j,v in enumerate(expanded[1]) if v in profile),'expanded nice-pair compatibility')
print(json.dumps({'status':'passed','checks':checks,'family':family,'parity':par,'residue':d0,'split_prefix_blocks':split,'witness':p.name},indent=2,sort_keys=True))
