from core import require,hp,strong
import json,itertools
Q=[2,4,41,19,39,11]
def graph(p,q,m):
 sizes=[1,1,p,1,q,m];groups=[];s=0
 for k in sizes:groups.append(list(range(s,s+k)));s+=k
 rows=[0]*s
 for i,G in enumerate(groups):
  for u in G:
   for v in G:
    if u<v:rows[u]|=1<<v
   for j,H in enumerate(groups):
    if Q[i]>>j&1:
     for v in H:rows[u]|=1<<v
 return rows,groups

def formula(p,q,m):
 A=4*(3**p-2**p)*(2**q-1)
 B=2*3**p+3*2**p-6-(2*3**p-3*2**p)*2**q
 C=(2**p-3)*(2**q-3)-2
 return A*3**m+B*2**m+C

def main():
 controls=[]
 # Exactly the tensor interpolation grid in the proof.
 for p,q,m in itertools.product(range(1,4),range(1,3),range(1,4)):
  rows,G=graph(p,q,m);H=hp(rows)
  require(H==formula(p,q,m),('grid',p,q,m,H))
  require(strong(rows),'strong family')
  controls.append([p,q,m,H])
 extra=[]
 for p,q,m in [(1,4,2),(2,4,3),(3,3,3),(4,1,4),(4,3,2),(5,1,2),(2,2,6)]:
  rows,G=graph(p,q,m);H=hp(rows)
  require(H==formula(p,q,m),('outside grid',p,q,m,H));extra.append([p,q,m,H])
 rows,G=graph(2,2,1)
 # Seed uses each of the two length-two chains in the reverse label order.
 original=[2,12,145,149,99,143,175,19];mapping=[0,1,3,2,4,6,5,7]
 relabeled=[sum(((original[u]>>v)&1)<<mapping[v] for v in range(8))for u in range(8)]
 target=[0]*8
 for u in range(8):target[mapping[u]]=relabeled[u]
 require(rows==target,'literal audited clock seed quotient')
 # Uniform matching upper bounds: exhibit covers of the crossing edges.
 # For U (module2): right side singleton b (module1), left side singleton c(module3).
 # For V (module4): its only incoming external module is singleton c(module3).
 # For R (module5): every crossing edge has left endpoint b or c.
 certificates=[]
 for p,q,m in [(1,1,1),(2,2,3),(3,4,2)]:
  rows,G=graph(p,q,m)
  for module,bound,leftmods,rightmods in [(2,2,[3],[1]),(4,1,[],[3]),(5,2,[1,3],[])]:
   v=G[module][0];outside=[u for i,H in enumerate(G)if i!=module for u in H]
   O=[u for u in outside if rows[v]>>u&1];I=[u for u in outside if rows[u]>>v&1]
   left={u for i in leftmods for u in G[i]};right={u for i in rightmods for u in G[i]}
   require(len(left)+len(right)==bound,'cover size')
   require(all(u in left or w in right for u in O for w in I if rows[u]>>w&1),'uniform matching cover')
   certificates.append([p,q,m,module,bound])
 for p in range(1,25):
  for q in range(1,25):
   B=2*3**p+3*2**p-6-(2*3**p-3*2**p)*2**q
   require((B==0)==((p,q)==(2,2)),'cancellation controls')
 return {'status':'PASS','quotient':Q,'interpolation_grid':controls,'outside_grid':extra,
         'crossing_cover_controls':certificates,'unique_cancellation':[2,2]}
if __name__=='__main__':print(json.dumps(main(),sort_keys=True,indent=2))
