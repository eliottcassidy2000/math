from core import *
import json

def main():
 n=len(SEED);F=(1<<n)-1
 for i in range(n):
  require(not SEED[i]>>i&1,'loop')
  require(not SEED[i]&~F,'outside bit')
  for j in range(i):require(bool(SEED[i]>>j&1)!=bool(SEED[j]>>i&1),'tournament')
 require(SEED[ROOT]==sum(1<<v for v in PORTS),'physical ports')
 a=run_dp(SEED,ROOT);literal=literal_runs(SEED,ROOT)
 require(a[1:]==literal[1:],'independent full run enumeration')
 require(a==[59,179,180,60,0,0,0,0,0],'finite run certificate')
 O=PORTS;I=[2,3,5,6]
 cross=[(u,v)for u in O for v in I if SEED[u]>>v&1]
 require(cross==[(1,2),(1,3),(4,5),(4,6)],'triangle crossing certificate')
 require(strong(SEED),'strong seed')
 controls=[]
 for m in range(1,7):
  T=clone(SEED,ROOT,m);U,ps=birth(T,PORTS,'R',7,3)
  require(U==clone(SEED,ROOT,m+1),'literal returning graph')
  require(ps==PORTS,'same retained ports')
  require(strong(T),'strong clone')
  H=hp(T);require(H==60*3**m-1,'direct count')
  controls.append({'m':m,'order':len(T),'H':H,'ports':ps})
 # Independently test that swapping two internal fragments preserves each
 # literal three-run skeleton, and has no fixed points.
 count=0
 for w in permutations(range(7)):
  for gaps in combinations(range(8),3):
   word=[]
   for i in range(8):
    if i in gaps:word.append(ROOT)
    if i<7:word.append(w[i])
   if not all(SEED[u]>>v&1 for u,v in zip(word,word[1:])):continue
   p=[i for i,v in enumerate(word)if v==ROOT];u,v,z=p
   other=word[:u+1]+word[v+1:z]+[ROOT]+word[u+1:v]+word[z:]
   require(other!=word,'free internal-fragment involution')
   require(sorted(other)==sorted(word),'same labels')
   require(all(SEED[x]>>y&1 for x,y in zip(other,other[1:])),'preserved directed word')
   count+=1
 require(count==60,'all involution controls')
 return {'status':'PASS','rows':SEED,'root':ROOT,'ports':PORTS,'run_coefficients':a,
         'cross_edges':cross,'literal_skeleton_candidates':5040*255,
         'involution_controls':count,'physical_controls':controls}
if __name__=='__main__':print(json.dumps(main(),sort_keys=True,indent=2))
