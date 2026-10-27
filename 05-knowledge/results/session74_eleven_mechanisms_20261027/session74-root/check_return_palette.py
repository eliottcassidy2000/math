import json

def need(b,m):
    if not b:raise RuntimeError(m)
p=11;H={1,3,4,5,9};N=set(range(1,p))-H
need({2*x%p for x in H}==N,'2 switches')
need({3*x%p for x in H}==H,'3 preserves')
need(pow(3,5,p)==1 and all(pow(3,k,p)!=1 for k in range(1,5)),'order3')
need(pow(4,4,p)==3,'phase4')
checks=0
for x in range(p):
 for y in range(p):
  if x==y:continue
  arc=(y-x)%p in H
  need(((2*(y-x))%p in H)!=arc,'reverse')
  need(((3*(y-x))%p in H)==arc,'preserve')
  checks+=2
mu={x:x+1 for x in range(p)}
for r in (2,3,11):
 pushed={}
 for x,w in mu.items():pushed[r*x%p]=pushed.get(r*x%p,0)+w
 need(sum(pushed.values())==66,'mass')
 need(len(pushed)==(1 if r==11 else 11),'image')
 if r==11:need(pushed=={0:66},'collapse')
print(json.dumps({'status':'PASS','arc_checks':checks,'repeat_multipliers':[2,3,11],'presentation_clocks':[10,15,55],'mass':66},sort_keys=True,indent=2))
