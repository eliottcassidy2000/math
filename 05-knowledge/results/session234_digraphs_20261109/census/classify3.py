import itertools,json,collections
out={}
for n in range(1,5):
 cs={}
 for code in range(3**(n*(n-1)//2)):
  a=[[0]*n for _ in range(n)];z=code
  for i in range(n):
   for j in range(i+1,n):
    s=z%3;z//=3;a[i][j]=int(s!=1);a[j][i]=int(s!=0)
  h=sum(all(a[p[i]][p[i+1]] for i in range(n-1)) for p in itertools.permutations(range(n)))
  if h!=3:continue
  r=[x[:] for x in a]
  for k in range(n):
   for i in range(n):
    for j in range(n):r[i][j]|=r[i][k]&r[k][j]
  if not all(r[i][j] for i in range(n) for j in range(n)):continue
  canon=min(''.join(str(a[p[i]][p[j]]) for i in range(n) for j in range(n)) for p in itertools.permutations(range(n)))
  if canon not in cs:cs[canon]={'code':code,'matrix':a,'labelled_count':0}
  cs[canon]['labelled_count']+=1
 out[n]=cs;print(n,cs)
json.dump(out,open('strong_h3.json','w'),indent=2)
