import itertools,json
BITS=['000110001000000','000110011000000']
def statistics(s):
 n=len(s);aut=[];anti=[];hp=blue=0
 for p in itertools.permutations(range(n)):
  if all(s[p[i]][p[j]]==s[i][j] for i in range(n) for j in range(i+1,n)):aut.append(p)
  if all(s[p[i]][p[j]]==-s[i][j] for i in range(n) for j in range(i+1,n)):anti.append(p)
  if all(s[p[i]][p[i+1]]==1 for i in range(n-1)):
   hp+=1;blue+=all(s[p[n-1-i]][p[n-1-j]]==-s[p[i]][p[j]] for i in range(n) for j in range(i+1,n))
 return {'H':hp,'Aut':len(aut),'b':blue//len(aut)},aut,anti
out=[]
for bits in BITS:
 n=6;s=[[0]*n for _ in range(n)]
 for b,(i,j) in zip(bits,itertools.combinations(range(n),2)):s[i][j]=1 if b=='1' else -1;s[j][i]=-s[i][j]
 d=list(map(sum,s));r=[row+[-d[i]] for i,row in enumerate(s)]+[d+[0]]
 st,aut,anti=statistics(s);sr,ar,br=statistics(r)
 sig=next(p for p in anti if all(p[p[i]]==i for i in range(n)))
 c=sum(all(g[sig[i]]==sig[g[i]] for i in range(n)) for g in aut)
 pairs=[(i,sig[i]) for i in range(n) if i<sig[i]]
 f_even=f_odd=flux=0;endpoint=[0]*n
 for choice in itertools.product(*pairs):
  for p in itertools.permutations(choice):
   if all(s[p[i]][p[i+1]]==1 for i in range(len(p)-1)):
    v=p[-1];endpoint[v]+=1;f_even+=s[v][sig[v]]==1;f_odd+=d[v]==-1;flux+=s[v][sig[v]]+d[v]
 if f_even!=c*st['b'] or f_odd!=c*sr['b'] or flux!=2*c*(st['b']-sr['b']):raise RuntimeError('gate identity')
 out.append({'bits':bits,'T':st,'completion':sr,'anti_involution':sig,'centralizer':c,'halfpath_endpoint_counts':endpoint,'even_gate':f_even,'completion_gate':f_odd,'gate_flux':flux})
print(json.dumps({'status':'PASS','records':out},indent=2))
