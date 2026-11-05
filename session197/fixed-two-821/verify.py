from pathlib import Path
from itertools import permutations
import subprocess,tempfile,json
from twoport_core import shapes,transitions,brute,physical,step,H
P=Path(__file__).resolve().parent

def need(ok,message):
 if not ok:raise RuntimeError(message)
def cycle_readout(v,c):
 return sum(v[i]for i in ((),(1,2),(3,5),(2,5))[c])
def cpp_table():
 seed=brute([2,0],[0,1]);s='const int hc[]={'+','.join(str(i)for i,a in enumerate(shapes(2))if len(a)==1)+'};\n'
 s+='const State seed={{'+','.join(map(str,seed))+'}};\nconst vector<vector<pair<int,int>>> tr[12]={\n'
 for c in range(4):
  for d in range(3):s+='{'+','.join('{'+','.join('{'+str(j)+','+str(m)+'}'for j,m in row)+'}'for row in transitions(2,c,d))+'},\n'
 return s+'};\n'
def python_closure(bound,limit):
 seed=brute([2,0],[0,1]);Q=[seed];seen={seed};values=set();edges=reject=0
 for v in Q:
  for c in range(4):
   C=cycle_readout(v,c)
   if C>limit:reject+=1;continue
   w=[0]*27
   for i,x in enumerate(v):
    for j,m in transitions(2,c,2)[i]:w[j]+=m*x
   w=tuple(min(bound+1,x)for x in w)
   h=H(w,2)
   if h>bound:continue
   edges+=1
   if C:values.add(h)
   if w not in seen:seen.add(w);Q.append(w)
 return dict(states=len(Q),admissible_edges=edges,cycle_rejected_edges=reject,support=sorted(values))
def graph(W):
 T=[2,0];ports=[0,1]
 for c in W:T,ports=physical(T,ports,c,2)
 need(ports==[0,1],'physical port identity');return T
def good(T,p):return all(T[a]>>b&1 for a,b in zip(p,p[1:]))
def literal_profile(T,r):
 n=len(T);h=e=s=R=0
 for p in permutations(range(n)):
  if good(T,p):h+=1;e+=p[-1]==r;s+=p[0]==r
  for cut in range(n):
   if p[cut]==r and good(T,p[:cut])and good(T,p[cut:]):R+=1
 V=[v for v in range(n)if v!=r];d=sum(good(T,p)for p in permutations(V))
 return h,e,R,s,d

def counts(T):
 n=len(T);N=1<<n;D=[[0]*n for _ in range(N)]
 for v in range(n):D[1<<v][v]=1
 for M in range(1,N):
  for v in range(n):
   if D[M][v]:
    for w in range(n):
     if not M>>w&1 and T[v]>>w&1:D[M|1<<w][w]+=D[M][v]
 h=sum(D[-1]);D=[[0]*n for _ in range(N)];D[1][0]=1
 for M in range(1,N):
  for v in range(n):
   if D[M][v]:
    for w in range(n):
     if not M>>w&1 and T[v]>>w&1:D[M|1<<w][w]+=D[M][v]
 return h,sum(D[-1][v]for v in range(n)if T[v]&1)

def main():
 need(cpp_table()==(P/'automaton.inc').read_text(),'serialized physical transfer differs')
 saved=json.loads((P/'SUPPORT.json').read_text());small=[]
 with tempfile.TemporaryDirectory()as tmp:
  for flag in ('-O2','-O0'):
   exe=Path(tmp)/('support'+flag)
   subprocess.run(['g++','-std=c++17',flag,str(P/'support.cpp'),'-o',str(exe)],check=True)
   out=json.loads(subprocess.check_output([str(exe),'821','5000000','3'],text=True))
   need(out==saved and out['complete']and out['processed']==out['states'],'complete source closure mismatch')
   if flag=='-O2':
    for limit in (1,2,3):
     got=json.loads(subprocess.check_output([str(exe),'21','5000000',str(limit)],text=True));py=python_closure(21,limit)
     need(all(got[k]==py[k]for k in ('states','admissible_edges','cycle_rejected_edges')),'separate Python closure')
     need(sorted(map(int,got['witnesses']))==py['support'],'separate Python support');small.append(dict(cycle_limit=limit,**py))
 need(saved['missing_odds']==[7,21,821]and saved['states']==374692,'frozen target exclusion')
 witness_checks=0;max_order=0
 for h,w in saved['witnesses'].items():
  v=brute([2,0],[0,1]);W=[];C=0
  for c,d in w['word']:
   need(d==2,'fixed physical ports');C=cycle_readout(v,c);need(C<=3,'forbidden witness prefix');v=step(v,2,c,d);W.append(c)
  T=graph(W);max_order=max(max_order,len(T))
  need(H(v,2)==int(h)and C==w['C']and 1<=C<=3,'uncapped physical witness replay')
  need(all(T[a]>>b&1 for a in range(2,len(T))for b in range(a+1,len(T))),'feedback deletion')
  witness_checks+=1
 seed=graph([0,0,3,1]);need(seed==[14,44,56,48,35,1],'literal seed rows')
 need(literal_profile(seed,0)==(37,5,49,4,5),'literal seed path/cut profile')
 need(counts(seed)==(37,4),'seed graph counts')
 family=[]
 for k in range(0,31):
  W=[0,0,3,1]+[0]*k+[1];v=brute([2,0],[0,1]);C=0
  for c in W:C=cycle_readout(v,c);v=step(v,2,c,2)
  want=81 if k==0 else 37+49*2**k
  need(H(v,2)==want and C==4,'uncapped infinite-family controls')
  if k<=6:
   T=graph(W);need(counts(T)==(want,4),'direct subset family controls');family.append(dict(k=k,H=want,C=4,rows=T))
 target=graph([0,0,3,1,0,0,0,0,1]);need(counts(target)==(821,4),'literal H821 graph')
 out=dict(status='PASS',closure={k:v for k,v in saved.items()if k!='witnesses'},python_small_closures=small,uncapped_witness_replays=witness_checks,maximum_stored_witness_order=max_order,literal_seed_profile=[37,5,49,4,5],family_uncapped_controls=31,family_direct_subset_controls=family,target_rows=target,target_H=821,target_C=4)
 s=json.dumps(out,sort_keys=True,indent=2)+'\n';(P/'VALIDATION.json').write_text(s);print(s,end='')
if __name__=='__main__':main()
