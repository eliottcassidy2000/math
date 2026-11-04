"""Exact strictly positive formal-profile coarse blindness; no tournament claim."""
from pathlib import Path
from fractions import Fraction
from collections import defaultdict
import json,hashlib,importlib.util
ROOT=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('rank',ROOT.parent/'coarse_response_rank'/'verify.py');v=importlib.util.module_from_spec(spec);spec.loader.exec_module(v)
def need(ok,msg):
 if not ok:raise RuntimeError(msg)
def words(profiles,a,c,reverse=False):
 Q=[4,1,2] if reverse else [2,4,1];R=Q+[1<<4,0]
 for u in range(3):
  if a>>u&1:R[u]|=1<<3
  else:R[3]|=1<<u
  if c>>u&1:R[u]|=1<<4
  else:R[4]|=1<<u
 cap=[len(q)-1 for q in profiles];states={};answer=0
 for u in range(5):
  k=[0]*3;b=0
  if u<3:k[u]=1
  else:b=1<<(u-3)
  states[(*k,b,u)]=1
 for _length in range(1,sum(cap)+3):
  new=defaultdict(int)
  for (i,j,k,b,last),count in states.items():
   ks=[i,j,k]
   if b==3 and min(ks):answer+=count*profiles[0][i]*profiles[1][j]*profiles[2][k]
   nxt=R[last]
   while nxt:
    bit=nxt&-nxt;nxt-=bit;u=bit.bit_length()-1;kk=ks[:];bb=b
    if u<3:
     if kk[u]==cap[u]:continue
     kk[u]+=1
    else:
     bit=1<<(u-3)
     if bb&bit:continue
     bb|=bit
    new[(*kk,bb,u)]+=count
  states=new
 return answer
def main():
 raw=json.loads((ROOT/'profile_data.json').read_text());A=raw['A'];B=raw['B'];C=[Fraction(x) for x in raw['C']];D=raw['C_scaling_denominator'];Ci=raw['C_integer'];qs=[A,B,Ci]
 need(all(q[0]==0 and all(x>0 for x in q[1:]) for q in [A,B,C]),'strict positive support')
 need(A[1]==B[1]==C[1]==1,'normalized first entries')
 need(Ci==[int(x*D) for x in C] and all((x*D).denominator==1 for x in C),'denominator clearing')
 need(B[2]%2==1 and A[3]%6!=0,'explicit failure of tournament profile divisibility')
 Gt=list(map(list,zip(*v.G)));records=[]
 for a in range(8):
  for c in range(8):
   f=v.functional(v.numerator(v.G,a,c),qs);g=v.functional(v.numerator(Gt,a,c),qs)
   need(f==g and f>0,'exact all-query formal equality')
   ff=words(qs,a,c);gg=words(qs,a,c,True);need((f,g)==(ff,gg),'independent run-word model agreement')
   records.append({'block_query':[a,c],'integer_profile_counts':[f,g],'normalized_counts':[str(Fraction(f,D)),str(Fraction(g,D))]})
 return {'status':'PASS','scope':'Strictly positive finite rational/integer formal profiles with all64 coarse contrasts zero; these are not tournament profiles and this is not a tournament or full named-response collision.','orders':[len(A)-1,len(B)-1,len(C)-1],'normalized_first_entries':[1,1,1],'all64_equal':True,'independent_run_word_counts':128,'records':records,'source_sha256':{str(p.relative_to(ROOT.parent)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [ROOT/'PROOF.md',ROOT/'verify.py',ROOT/'profile_data.json',ROOT.parent/'coarse_response_rank'/'verify.py']}}
if __name__=='__main__':print(json.dumps(main(),sort_keys=True,indent=2))
