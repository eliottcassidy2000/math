"""Matched C0,C1 workloads; descriptive single-run timings, not complexity proof."""
from backbone import *
rows=[]
for n in [32,48,64,80]:
 t=time.perf_counter();u,c,stats=envelope(n,(3,5),True);tdp=time.perf_counter()-t
 t=time.perf_counter();ur=F(0);cr=F(0);count=0
 for p in parts(n):
  w,qs=sector(p,(3,5));ur+=w;cr+=w*sum(qs,F(0));count+=1
 tref=time.perf_counter()-t
 assert (u,c)==(ur,cr)
 rows.append({'n':n,'backbone':[3,5],'marked_dp_seconds':tdp,'partition_seconds':tref,'partition_count':count,**stats})
print(json.dumps({'status':'PASS','rows':rows},indent=2))
