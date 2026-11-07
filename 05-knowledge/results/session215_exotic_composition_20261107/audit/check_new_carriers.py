from audit_support import *
results=[]
spec=[(4225,6,'JO',(2,20,41,67,77,27,167,63)),(4225,4,'JO',(6,12,24,49,195,23,47,111)),(4357,3,'ZJO',(10,20,65,22,37,15,59)),(4357,2,'ZJO',(2,48,11,67,45,77,23))]
records=[]
for x,root,word,seed in spec:
 rows=car.remark(seed,root);p=car.profile(rows)
 for c in word:rows=car.graph_move(rows,c);p=car.move(p,c)
 n=len(rows);core=[2,8,3,5];out=list(rows)+[0]*3
 for a in range(4):
  av=0 if a==0 else n+a-1
  for b in range(4):
   if core[a]>>b&1:out[av]|=1<<(0 if b==0 else n+b-1)
  if a:out[av]|=((1<<n)-1)^1
 r=p[1]%3;j=(p[1]-r)//3
 require(r in [1,2] and j<=8 and car.base(p,j)==x,'stagger labels')
 records.append(dict(seed=seed,r=r,side=r,word=word,rows=rows,profile=list(p),compiled_rows=out,compiled_count=3*x+2*r,X=x))
for w in records:
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
