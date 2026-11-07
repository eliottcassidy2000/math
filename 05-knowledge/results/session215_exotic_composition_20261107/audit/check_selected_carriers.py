from audit_support import *
results=[]
for w in json.loads((CARRIERS/'WITNESSES.json').read_text()):
 rows=w['rows'];start,h=subset_h(rows);n=len(rows);full=(1<<n)-1
 rev=[sum((not(rows[i]>>j&1))<<j for j in range(n) if i!=j) for i in range(n)]
 end,_=subset_h(rev)
 def C(mask):
  ans=0;s=mask
  while True:
   ans+=h(s)*h(mask^s)
   if not s:break
   s=(s-1)&mask
  return ans
 L=sum(end(s,0)*h(full^s) for s in range(1,full+1,2))
 R=sum(start(s,0)*h(full^s) for s in range(1,full+1,2))
 T=sum(start(s,0)*C(full^s) for s in range(1,full+1,2))
 p=[h(full),start(full,0),end(full,0),L,R,h(full^1),C(full),C(full^1),T]
 require(p==w['profile'],'nine profile')
 _,ch=subset_h(w['compiled_rows']);actual=ch((1<<len(w['compiled_rows']))-1)
 require(actual==w['compiled_count']==3*w['X']+2*w['r'],'compiled graph')
 results.append({'seed':w['seed'],'side':w['side'],'word':w['word'],'profile':p,'X':w['X'],'compiled_count':actual})
print(json.dumps({'status':'PASS','selected_independent_complete_profiles':results},indent=2,sort_keys=True))
