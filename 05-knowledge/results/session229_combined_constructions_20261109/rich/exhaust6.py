from search import *
best=0;hits=[]
for code in range(1<<15):
 a=[0]*6;k=0
 for i in range(6):
  for j in range(i+1,6):
   if code>>k&1:a[i]|=1<<j
   else:a[j]|=1<<i
   k+=1
 b=inspect(a)
 if b:
  score=b['distinct']/b['ratio']
  if score>best:
   best=score; print(json.dumps(dict(code=code,**b)),flush=True)
  if b['distinct']>=4:hits.append(dict(code=code,**b))
Path(__file__).with_name('exhaust6_hits.json').write_text(json.dumps(hits))
