"""Rebuild the finite universal-pair receipt from complete filtered root-cut files."""
import json
from pathlib import Path
seen={};counts={};matches=[]
for n in range(1,10):
 path=Path('CANDIDATES10.txt') if n==9 else Path('replay'+str(n))/'CANDIDATES10.txt'
 count=0
 for line in path.read_text().splitlines():
  v=list(map(int,line.split()));h,s,e,L,R,d,C,b,T=v[2:];r=s%3;k=(h,e,L,R,d,C,b,T,s+2*r)
  item={'n':n+1,'ix':v[0],'mask':v[1],'r':r}
  if k in seen and seen[k]['r']!=r:matches.append([seen[k],item])
  seen[k]=item;count+=1
 counts[n+1]=count
print(json.dumps({'status':'COMPLETE_CONDITIONAL_ON_RECIPROCITY_THEOREM','max_order':10,'candidate_counts_by_order':counts,'keys':len(seen),'universal_pair_matches':matches},indent=2))
