#!/usr/bin/env python3
import itertools,json,sys
from pathlib import Path
p=Path(sys.argv[1]) if len(sys.argv)>1 else Path('../path-techniques/WITNESSES.json')
out=[]
for w in json.loads(p.read_text()):
 a=w['rows'];n=len(a);mark=w['mark'];h=e=0
 for word in itertools.permutations(range(n)):
  if all(a[u]>>v&1 for u,v in zip(word,word[1:])):h+=1;e+=word[-1]==mark
 if not (a[mark]==1<<w['sole_successor'] and h==w['H'] and e==w['e_a'] and h-2*e==3):raise RuntimeError(w)
 out.append({'n':n,'H':h,'endpoint':e,'defect':3})
print(json.dumps({'status':'PASS','method':'literal permutations','cases':out},indent=2,sort_keys=True))
