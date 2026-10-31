#!/usr/bin/env python3
"""Universal terminal-cycle replay. Input compatibility retained; no output-pair claim."""
import json,sys
from pathlib import Path
from math import isqrt

def verify(a):
 checks=0
 def require(x,msg):
  nonlocal checks
  checks+=1
  if not x:raise RuntimeError(msg)
 require(a['status']=='witness' and a['family']=='B-final-cycle-prefix4','certificate type')
 p=a['parity'];d=a['residue'];cuts=a['split_prefix_blocks'];word=a['word']
 require(type(p)is int and p in(0,1) and type(d)is int and 24<=d<=72,'state')
 require(isinstance(cuts,list) and all(type(c)is int and -24<=c<=24 for c in cuts),'cuts range')
 require(len(cuts)==len(set(cuts)) and 4 not in cuts,'distinct unprotected cuts')
 singles=set();blocks=set();heads=set();length_a=length_b=0;first=last=None
 for atom in word:
  kind=atom[0]
  if kind=='single':
   require(len(atom)==2 and type(atom[1])is int,'singleton format');v=atom[1]
   require(1<=v<=24 and v not in singles,'singleton ownership');singles.add(v);left=right=v;length_b+=1
  elif kind=='head':
   require(len(atom)==2 and type(atom[1])is int,'head format');c=atom[1]
   require(c in cuts and c not in heads,'head ownership');heads.add(c);left=right=49-c;length_b+=1
  else:
   require(kind in('block','tail') and len(atom)==4,'block format');_,c,i,o=atom
   require(all(type(x)is int for x in(c,i,o)) and -24<=c<=24 and c not in blocks,'block ownership');blocks.add(c)
   require((kind=='tail')==(c in cuts),'declared split semantics')
   require(i==int(49+(-1)**(p+1)*c<=d) and o in(-1,1),'forced member/orientation')
   ell=(p+i)%2;left=147+c if kind=='tail' else 49-c;right=49*(8 if ell==0 else 10)+(-1)**ell*c
   if o==-1:left,right=right,left
   length_a+=1;length_b+=i-(kind=='tail')
  if first is None:first=left
  if last is not None:require(isqrt(last+left)**2==last+left,'square connector')
  last=right
 require(singles==set(range(1,25)) and blocks==set(range(-24,25)) and heads==set(cuts),'exact atom partition')
 require((length_a,length_b)==(49,d),'symbolic length49n+d')
 require(word[:4]==[['single',v] for v in[1,3,6,19]],'fixed low prefix')
 require(word[4]==['block',4,int(49+(-1)**(p+1)*4<=d),1],'intact recursive block')
 terminal=8 if(p+d)%2==0 else 15
 require(word[-1]==['single',terminal],'closing port')
 require(first==1 and isqrt(last+first)**2==last+first,'closing square')
 return {'status':'passed','checks':checks,'parity':p,'residue':d,'prefix_splits':len(cuts),'symbolic_length':[length_a,length_b],'terminal':terminal}

if __name__=='__main__':print(json.dumps(verify(json.loads(Path(sys.argv[1]).read_text())),indent=2))
