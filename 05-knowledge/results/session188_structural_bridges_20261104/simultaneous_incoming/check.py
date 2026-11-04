"""Literal integration checks for incoming disjoint isolation and mixed terms."""
import sys,itertools,json
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent.parent/'response'))
from verify import hp,extend,need

def rev(a,ms):
 b=list(a)
 for m in ms:
  for u,v in itertools.combinations(range(len(a)),2):
   if (m>>u&1) and (m>>v&1):b[u]^=1<<v;b[v]^=1<<u
 return b

def main():
 T=[14,12,56,48,35,3];mods=[3,12,48];states=[rev(T,[m for i,m in enumerate(mods) if bits>>i&1]) for bits in range(8)]
 decoders=[]
 for bits,a in enumerate(states):
  for i,(x,c) in enumerate([(48,49),(3,7),(12,28)]):
   F=hp(extend(a,x,c));need(F==(483 if bits>>i&1 else 393),'incoming fixed decoder');decoders.append([bits,i,x,c,F])
 # Four-state quartic interaction, using only actual finite-difference queries.
 mon=5|(18<<6);s=mon;total=0;terms=[]
 while True:
  x=s&63;c=s>>6;vals=[hp(extend(states[k],x,c)) for k in range(4)];interaction=vals[0]-vals[1]-vals[2]+vals[3];sign=(-1)**(4-s.bit_count());total+=sign*interaction;terms.append({'a':x,'c':c,'four_counts':vals,'sign':sign})
  if s==0:break
  s=(s-1)&mon
 need(total==-2,'incoming mixed quartic')
 # Smaller actual-query nonadditivity controls supplied by the parent.
 A=[0,1,3,7];small=[]
 for m,n,x,expected in [(3,12,5,[19,11,15,9]),(3,7,9,[15,9,5,9])]:
  vals=[hp(extend(rev(A,ops),x,0)) for ops in ([],[m],[n],[m,n])];need(vals==expected,'small nonadditive query');mixed=vals[0]-vals[1]-vals[2]+vals[3];need(mixed%4==2,'small mixed residue');small.append({'T':A,'modules':[m,n],'query':[x,0],'four_counts':vals,'mixed_difference':mixed})
 return {'status':'PASS','incoming_commit':'b6135263ba67ee2908df7625c67a88a3df8fbe03','decoder_queries':decoders,'mixed_quartic':{'a_mask':5,'c_mask':18,'coefficient':total,'sixteen_terms':terms},'small_nonadditivity':small,'scope':'Exact integration controls; disjoint isolation is credited to incoming S1-S3. No general additivity or covered-laminar theorem is asserted.'}
if __name__=='__main__':print(json.dumps(main(),sort_keys=True,indent=2))
