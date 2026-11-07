from probe import hp
import json
checks=[]
for p,s,c,A,C in [((1,1,2,2),(1,),0,60,-1),((1,1,0,2,2),(1,),0,120,1),((2,2,0),(3,2,2),1,72,1)]:
 for m in range(5):
  actual=hp(p+(c,)*m+s)
  if actual!=A*3**m+C:raise RuntimeError((p,s,c,m,actual))
  checks.append({'prefix':p,'suffix':s,'letter':c,'m':m,'H':actual})
print(json.dumps({'independent_subset_checks':len(checks),'checks':checks},indent=2))
