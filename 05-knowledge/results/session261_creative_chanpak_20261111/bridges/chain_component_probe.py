"""Focused controls of chain-run lifting, factorial components, and a stronger scout."""
from path_cover_probe import H,sub
import random,itertools,json,math
from pathlib import Path

def path_counts(r,A,B):
 n=len(r); dp=[[0]*n for _ in range(1<<n)]
 for a in A:dp[1<<a][a]=1
 total=0
 for s in range(1,1<<n):
  if s.bit_count()%2==0:total+=sum(dp[s][b] for b in B)
  for v in range(n):
   if dp[s][v]:
    q=r[v]&~s
    while q:
     bit=q&-q;q-=bit;w=bit.bit_length()-1;dp[s|bit][w]+=dp[s][v]
 return total

def build(words,rng):
 n=sum(map(len,words));color=[];comp=[]
 for j,w in enumerate(words):color+=list(w);comp +=[j]*len(w)
 r=[0]*(n+1)
 for u,v in itertools.combinations(range(n),2):
  if color[u]==color[v]:win=u if rng.randrange(2) else v
  elif comp[u]==comp[v]:win=u
  else:win=u if color[u]=='B' else v
  loser=v if win==u else u;r[win]|=1<<loser
 for v,c in enumerate(color):
  if c=='A':r[n]|=1<<v
  else:r[v]|=1<<n
 return tuple(r),color

def E(w):
 return sum(s.bit_count()%2==0 and w[(s&-s).bit_length()-1]=='A' and w[s.bit_length()-1]=='B' for s in range(1,1<<len(w)))

def main():
 rng=random.Random(2611707); cases=[]
 shapes=[['ABAB'],['AABABB'],['BAABABB'],['ABAB','AB'],['ABAB','ABAB'],['AB','AB','AB'],['BAAB','AABB']]
 for words in shapes:
  es=list(map(E,words));bound=sum(2*math.factorial(s.bit_count())*math.prod(es[i] for i in range(len(es)) if s>>i&1) for s in range(1,1<<len(es)))
  vals=[]
  for _ in range(12):
   r,c=build(words,rng);u=sub(r,list(range(len(c)))); delta=H(r)-H(u)
   if delta<bound:raise RuntimeError('chain block bound')
   for start,w in [(sum(map(len,words[:i])),w) for i,w in enumerate(words)]:
    rr=sub(r,list(range(start,start+len(w))));pc=path_counts(rr,[i for i,x in enumerate(w) if x=='A'],[i for i,x in enumerate(w) if x=='B'])
    if pc<E(w):raise RuntimeError('chain run lifting')
   vals.append(delta)
  cases.append({'words':words,'E':es,'factorial_bound':bound,'deltas':vals})
 return {'seed':2611707,'trials_per_word_list':12,'cases':cases}
if __name__=='__main__':
 out=main();Path(__file__).with_name('chain_component_receipt.json').write_text(json.dumps(out,indent=2)+'\n');print('Chain-component controls passed: 84 rooted constructions and their component path counts.')

def single_cycle_scout():
 from path_cover_probe import rows,graph
 rng=random.Random(2611708);ans=[]
 for n in [4,5,6,7,8]:
  slacks=[]
  for _ in range(100):
   r=rows(n,rng.getrandbits(n*(n-1)//2));z=n-1;A,B,es,d=graph(r,z);u=sub(r,list(range(n-1)));p=path_counts(u,A,B);slacks.append(1+2*p-d)
   if slacks[-1]<0:return {'status':'COUNTEREXAMPLE','rows':r,'root':z,'d':d,'P_even':p}
  ans.append({'rooted_order':n,'trials':100,'minimum_slack':min(slacks)})
 return {'status':'FINITE-EXACT evidence only; conjectural outside proved classes','seed':2611708,'cases':ans}
if __name__=='__main__':
 Path(__file__).with_name('single_cycle_scout.json').write_text(json.dumps(single_cycle_scout(),indent=2)+'\n')
