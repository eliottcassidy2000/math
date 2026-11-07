from inherited import *
from verify_witnesses import hp,compiler
import json
spec=[(4225,1531,6,'JO',(2,20,41,67,77,27,167,63)),(4225,2803,4,'JO',(6,12,24,49,195,23,47,111)),(4357,460,3,'ZJO',(10,20,65,22,37,15,59)),(4357,258,2,'ZJO',(2,48,11,67,45,77,23))]
results=[]
for x,ix,root,w,a in spec:
 rows=remark(a,root);p=profile(rows)
 for c in w:rows=graph_move(rows,c);p=move(p,c)
 need(p==profile(rows),'nine-state literal');r=p[1]%3;j=(p[1]-r)//3
 need(1<=r<=2 and 0<=j<=8 and base(p,j)==x,'labels/base')
 out=compiler(rows);h=hp(out);need(h==3*x+2*r,'compiler')
 results.append(dict(X=x,bank_index=ix,root=root,seed_rows=a,word=w,r=r,j=j,profile=p,rows=rows,compiled_rows=out,compiled_count=h))
for i in (0,2):
 a,b=results[i:i+2];p,q=tuple(a['profile']),tuple(b['profile']);js=a['j'],b['j']
 need((p[2],p[4]+p[7],p[5])==(q[2],q[4]+q[7],q[5]),'synchronization')
 V=p[4]+p[7];d=p[5];offset=a['X']-V+2*d-p[2];coef=V-d
 for k in range(1,31):
  p,q=move(p,'Z'),move(q,'Z')
  need(base(p,js[0])==base(q,js[1])==offset+coef*(1<<k),'ray')
 a['pure_Z_ray']={'coefficient':coef,'offset':offset,'k_min':1}
 # All literal common words to length2, including terminal sinks.
 def check(p,q,aa,bb,depth):
  need(profile(aa)==p and profile(bb)==q,'literal common word')
  need(base(p,js[0])==base(q,js[1]),'equal bases')
  if depth==2:return
  for c in 'OZ':check(move(p,c),move(q,c),graph_move(aa,c),graph_move(bb,c),depth+1)
 check(tuple(a['profile']),tuple(b['profile']),tuple(a['rows']),tuple(b['rows']),0)
# All old O/Z words capable of realizing X<=4357 have <=6 letters after
# deleting terminal O, by inherited order cap floor(log2 M)+2 =14.
old=set()
for pair,js,root,aa,bb in PAIRS:
 p=profile(remark(aa,root))
 def oldwalk(p,depth):
  old.add(base(p,js[0]))
  if depth==6:return
  for c in 'OZ':oldwalk(move(p,c),depth+1)
 oldwalk(p,0)
need(4225 not in old and 4357 not in old,'old orbit exclusion')
print(json.dumps({'status':'PASS','carriers':results,'old_OZ_bases_excluded':[4225,4357],'old_word_bound':6},indent=2))
