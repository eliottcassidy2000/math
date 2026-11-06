"""D-interface controls independent of the main readout formulas."""
import hashlib,json
from pathlib import Path
from verify import bank,state,F,graph,counts,require,step

def main():
 cap=5000;limit=3333;B=bank(cap)
 plus={D for D,E,U,d in B if 3<=D<=limit and U-D-E==1}
 minus={D for D,E,U,d in B if 3<=D<=limit and U-D-E==-1}
 expected_plus=set();expected_minus=set()
 for k in range(15):
  for base in [6,10]+[6*((1<<j)+1)for j in range(1,15)]:
   n=base*(1<<k)-1
   if n<=limit:expected_plus.add(n)
  if k>=1:
   for n in[(1<<k)+1,3*((1<<k)+1)]:
    if n<=limit:expected_minus.add(n)
 require(plus==expected_plus,'complete D-interface plus parents')
 require(minus==expected_minus,'complete D-interface minus parents')
 require(plus&minus=={5,9,17},'D-interface both children')
 for D,E,U,d in B:
  require(2*(D-E+U)>=3*D+1,'D-changing growth')
 for n in range(3,100000,2):require(9*n+5>4*(2*n+1),'two D-changing steps exceed both children')
 literals=[]
 for n,wm,wp in[(5,'Z','O'),(9,'ZZ','ZO'),(17,'ZZZ','ZJO')]:
  for sign,w in[(-1,wm),(1,wp)]:
   parent=w+'Z';child=parent+'J'
   require(state(parent)[0]==n,'D parent')
   require(state(child)[0]==2*n+sign,'D child')
   parentHC=counts(graph(parent));childHC=counts(graph(child))
   require(parentHC[0]==n and childHC==(2*n+sign,1),'literal D interface graphs')
   literals.append(dict(parent=n,sign=sign,parent_word=parent,parent_counts=parentHC,child_word=child,child_counts=childHC))
 print(json.dumps(dict(status='PASS',physical_cover_cap=cap,complete_parent_count_cap=limit,
  profiles=len(B),plus_count=len(plus),minus_count=len(minus),both=sorted(plus&minus),literal_graphs=literals,
  source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()),indent=2,sort_keys=True))
if __name__=='__main__':main()
