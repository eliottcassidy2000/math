import json,itertools
from pathlib import Path
base=Path(__file__).resolve().parent.parent
bank=json.loads((base/'session30-almost-regular/census8_bank.json').read_text())['almostregular8'];out=[]
for row in bank:
 n=8;s=[[0]*n for _ in range(n)]
 for b,(i,j) in zip(row['bits'],itertools.combinations(range(n),2)):s[i][j]=1 if b=='1' else -1;s[j][i]=-s[i][j]
 d=list(map(sum,s));lo=[i for i in range(n) if d[i]<0];hi=[i for i in range(n) if d[i]>0];sig=None
 for p in itertools.permutations(hi):
  a=list(range(n))
  for i,j in zip(lo,p):a[i]=j;a[j]=i
  if all(s[a[i]][a[j]]==-s[i][j] for i in range(n) for j in range(i+1,n)):sig=a;break
 if sig is None:continue
 c=0
 for pl in itertools.permutations(lo):
  for ph in itertools.permutations(hi):
   a=list(range(n))
   for old,new in zip(lo+hi,pl+ph):a[old]=new
   if all(s[a[i]][a[j]]==s[i][j] for i in range(n) for j in range(i+1,n)) and all(a[sig[i]]==sig[a[i]] for i in range(n)):c+=1
 pairs=[(i,sig[i]) for i in range(n) if i<sig[i]];even=odd=flux=0;ends=[0]*n
 for ch in itertools.product(*pairs):
  for p in itertools.permutations(ch):
   if all(s[p[i]][p[i+1]]==1 for i in range(3)):
    v=p[-1];ends[v]+=1;even+=s[v][sig[v]]>0;odd+=d[v]<0;flux+=s[v][sig[v]]+d[v]
 if even%c or odd%c or flux!=2*(even-odd):raise RuntimeError('integrality/flux')
 out.append({'bits':row['bits'],'sigma':sig,'centralizer':c,'endpoint_halfpaths':ends,'f_even':even,'f_completion':odd,'blue_even':even//c,'blue_completion':odd//c,'flux':(even-odd)//c})
if len(out)!=15:raise RuntimeError('self-converse census')
print(json.dumps({'status':'PASS','universe':'all15 self-converse almost-regular order8 classes','changed_classes':sum(r['flux']!=0 for r in out),'flux_distribution':{str(v):sum(r['flux']==v for r in out) for v in sorted({r['flux'] for r in out})},'records':out},indent=2))
