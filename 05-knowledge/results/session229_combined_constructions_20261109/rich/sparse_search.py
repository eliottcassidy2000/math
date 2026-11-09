from search import *
rng=random.Random(2291210);hits=[];best=0
for n in range(7,13):
 for t in range(600):
  a=[sum(1<<j for j in range(i+1,n)) for i in range(n)]
  for i,j in rng.sample([(i,j) for i in range(n) for j in range(i+1,n)],1+t%6):a[i]^=1<<j;a[j]^=1<<i
  b=inspect(a)
  if b:
   b.update(n=n,trial=t);hits.append(b);score=b['distinct']/b['ratio']
   if score>best:
    best=score;print(json.dumps(b),flush=True);Path(__file__).with_name('sparse_best.json').write_text(json.dumps(b,indent=2))
Path(__file__).with_name('sparse_hits.json').write_text(json.dumps(hits,indent=2))
