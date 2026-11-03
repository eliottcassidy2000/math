import bank65_check as s
import frontier as f
from collections import Counter
c=s.c
D=c.W['L2']+('b',)+c.W['L2']
def count(rows):
 n=len(rows);width=max([v-u for v in range(n) for u in range(v) if rows[v]>>u&1]+[1]);states=Counter({():1})
 for v,row in enumerate(rows):
  mask=sum(1<<(v-u-1) for u in range(v) if row>>u&1);nxt=Counter()
  for st,num in states.items():
   for ns,k in f.step(st,mask,width).items():nxt[ns]+=num*k
  states=nxt
 return sum(num for st,num in states.items() if len(st)==1),width
for word,expected in [(('q1',),3),(('q5','a'),13),(('q5','a')+D,12784595),(('q1','b')+D,38353833)]:
 rows=c.graph(word);h,width=count(rows);ha,_=count(c.add_a(rows));haa,_=count(c.add_a(c.add_a(rows)));state=c.vm((1,0,0),c.mat(word));c.check((h,ha,haa)==(state[0],state[0]+2*state[1],state[0]+2*(state[1]+state[2])),'full-state literal frontier');c.check(h==expected,'frontier H');c.check(all(not(rows[i]>>i&1) and all(((rows[i]>>j)&1)+((rows[j]>>i)&1)==1 for j in range(len(rows)) if j!=i) for i in range(len(rows))),'literal tournament');print('word',word,'vertices',len(rows),'port',len(rows)-1,'H',h,'width',width,'rows',rows)
w2=c.observed(('b','a')); v=c.vm((13,5,9),c.D)
c.check(w2==(5,2,3),'second actual predecessor')
c.check(tuple(w2[i]-(3,3,3)[i] for i in range(3))==c.t,'predecessor state difference')
lhs=c.vm(c.vm(w2,c.M['b']),c.D)
c.check(lhs==tuple(3*v[i]+28*c.t[i] for i in range(3)),'second mixed physical return')
for word,expected in [(('b','a'),5),(('b','a','b')+D,38353841)]:
 rows=c.graph(word);h,width=count(rows);ha,_=count(c.add_a(rows));haa,_=count(c.add_a(c.add_a(rows)));state=c.vm((1,0,0),c.mat(word));c.check((h,ha,haa)==(state[0],state[0]+2*state[1],state[0]+2*(state[1]+state[2])),'full-state literal frontier 2');c.check(h==expected,'frontier H 2');c.check(all(not(rows[i]>>i&1) and all(((rows[i]>>j)&1)+((rows[j]>>i)&1)==1 for j in range(len(rows)) if j!=i) for i in range(len(rows))),'literal tournament 2');print('word',word,'vertices',len(rows),'port',len(rows)-1,'H',h,'width',width,'rows',rows)
