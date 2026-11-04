import json,itertools,subprocess,hashlib
from pathlib import Path
P=Path(__file__).resolve().parent;G=P.parent/'closure_generation'
def need(x):
 if not x:raise RuntimeError('generation audit check failed')
def run(rows):
 n=len(rows);out=subprocess.run([str(P/'joint')],input=' '.join(map(str,[n,n,*rows,*range(n)])),capture_output=True,text=True,check=True).stdout
 h=None;fields={};j={}
 for line in out.splitlines():
  a=line.split()
  if a[0]=='H':h=int(a[1])
  elif a[0]=='PORT':fields[int(a[1])]=list(map(int,a[2:]))
  else:j[int(a[1]),int(a[2])]=int(a[3])
 return h,fields,j
def terminal(rows):
 n=len(rows);all=(1<<n)-1;ans=all
 for v in range(n):
  q=1<<v
  while True:
   new=q
   for u,r in enumerate(rows):
    if q>>u&1:new|=r
   if new==q:break
   q=new
  ans&=q
 return ans
out=[]
for word in ['22','24','42']:
 c=json.loads((G/f'candidate{word}.json').read_text());s=json.loads((P.parent/'closure'/f'candidate{word[0]}.json').read_text());ns=c['names'];rows=c['rows'];n=len(ns)
 need(ns==s['names']);need(n==18);need(c['prefix']==s['prefix']);need(c['target_H']==3*s['target_H']+int(word[1]))
 changed=[(ns[i],ns[j]) for i in range(n) for j in range(i+1,n) if (rows[i]^s['rows'][i])>>j&1]
 need({frozenset(z) for z in changed}=={frozenset(z) for z in c['all_changed_pairs']})
 for i in range(n):
  need(not(rows[i]>>i&1));need(rows[i]>>n==0)
  for j in range(i):need(((rows[i]>>j)&1)+((rows[j]>>i)&1)==1)
 need(terminal(rows)==(1<<n)-1)
 root=ns.index(39);ports=[ns[i]for i in range(n)if rows[i]>>root&1];need(sorted(ports)==sorted(c['in_ports']))
 h,f,j=run(rows);need(h==c['target_H']);need(f[root][3]==c['root_E2']);need(f[root][3]-h==c['root_defect'])
 sn=c['suffix_names'];sr=[sum(((rows[ns.index(u)]>>ns.index(v))&1)<<k for k,v in enumerate(sn)) for u in sn]
 hs,fs,js=run(sr);W=[sn.index(v)for v in ports];F=sum(fs[w][3]for w in W);J=sum(js.get((u,v),0)for u in W for v in W)
 need((hs,F,J)==(c['h'],c['F'],c['J']));need(h==hs-J+8*F);need(terminal(sr)&sum(1<<w for w in W))
 out.append({'word':word,'H':h,'all_source_names_preserved':True,'prefix':c['prefix'],'strong':True,'h':hs,'F':F,'J':J,'c':hs-J,'root_E2':f[root][3],'in_ports':ports,'changed_pairs':changed,'all_named_fields':{str(ns[i]):f[i]for i in range(n)},'joint_count':len(j),'joints':[{'end':ns[a],'start':ns[b],'value':v} for (a,b),v in j.items()],'suffix_port_joints':[{'end':sn[a],'start':sn[b],'value':js.get((a,b),0)} for a in W for b in W],'candidate_sha256':hashlib.sha256((G/f'candidate{word}.json').read_bytes()).hexdigest()})
print(json.dumps(out,indent=2));(P/'receipt.json').write_text(json.dumps(out,indent=2)+'\n')
