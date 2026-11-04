"""Exact finite profile certificate for prime-five quotients with triple children."""
from collections import defaultdict
from itertools import permutations
from pathlib import Path
import hashlib,json,subprocess,sys,tempfile
sys.path.insert(0,str(Path(__file__).resolve().parent.parent/'marked_pair_quotient'))
from verify import rows,prime,strong,induced,blow,opposite,hp,q2,need


def response(Q,a,c):
 n=5;R=list(Q)+[1<<6,0]
 for u in range(5):
  if a>>u&1:R[u]|=1<<5
  else:R[5]|=1<<u
  if c>>u&1:R[u]|=1<<6
  else:R[6]|=1<<u
 D={(1<<(2*u),u):1 for u in range(5)};D[1<<10,5]=1;D[1<<11,6]=1
 out=[0]*32
 for length in range(1,18):
  E=defaultdict(int)
  for (S,u),z in D.items():
   if S>>10==3:
    ks=[S>>(2*i)&3 for i in range(5)]
    if min(ks):
     mask=sum(1<<i for i,k in enumerate(ks) if k==1)
     out[mask]+=z*6**(5-mask.bit_count())
   vs=R[u]
   while vs:
    bit=vs&-vs;vs-=bit;v=bit.bit_length()-1
    if v<5:
     if S>>(2*v)&3==3:continue
     U=S+(1<<(2*v))
    else:
     b=1<<(10+v-5)
     if S&b:continue
     U=S|b
    E[U,v]+=z
  D=E
 return out

def evalp(p,mask):return sum(z*3**((s&mask).bit_count()) for s,z in enumerate(p))

def main():
 root=Path(__file__).resolve().parent
 triple_profiles={}
 for code in range(8):
  A=rows(3,code);h=hp(A);need(h in (1,3) and q2(A)==6,'complete triple profile classification')
  triple_profiles[code]=[0,h,6,6]
 representatives={};counts={};labelled=0
 for code in range(1024):
  Q=rows(5,code)
  if not prime(Q) or not strong(Q):continue
  key=min(tuple(induced(Q,p)) for p in permutations(range(5)))
  representatives.setdefault(key,Q);counts[key]=counts.get(key,0)+1;labelled+=1
 need(labelled==264 and len(representatives)==3,'complete prime-five quotient census')
 records=[]
 for key,Q in representatives.items():
  U=opposite(Q);a,c=1,2;P=response(Q,a,c);R=response(U,a,c)
  values=[]
  for mask in range(32):
   f,g=evalp(P,mask),evalp(R,mask);need(f!=g,'fixed query fails for a triple profile pattern')
   values.append({'cyclic_child_mask':mask,'T':f,'U':g,'difference':f-g})
  records.append({'quotient':Q,'labelled_quotient_count':counts[key],'a_blocks':a,'c_blocks':c,'T_profile_polynomial':P,'U_profile_polynomial':R,'profile_pattern_queries':values,'minimum_absolute_gap':min(abs(v['difference']) for v in values)})
 literal=[]
 with tempfile.TemporaryDirectory(prefix='prime5_triples_') as work:
  exe=Path(work)/'literal'
  subprocess.run(['g++','-std=c++17','-O3','-Wall','-Wextra','-Werror',str(root/'literal.cpp'),'-o',str(exe)],check=True)
  for record in records:
   Q=record['quotient']
   for mask in (0,31,21):
    blocks=[[2,4,1] if mask>>i&1 else [6,4,0] for i in range(5)]
    T,_=blow(Q,blocks);U,_=blow(opposite(Q),blocks);a,c=7,56;values=[]
    for A in (T,U):
     data='15 '+str(a)+' '+str(c)+' '+' '.join(map(str,A))+'\n'
     values.append(int(subprocess.check_output([str(exe)],input=data.encode())))
    expected=record['profile_pattern_queries'][mask]
    need(values==[expected['T'],expected['U']],'independent physical subset DP')
    literal.append({'quotient':Q,'cyclic_child_mask':mask,'T_rows':T,'U_rows':U,'a_mask':a,'c_mask':c,'T_count':values[0],'U_count':values[1]})
 print(json.dumps({'status':'PASS','scope':'Finite exact classification: every strong prime five-vertex quotient with every child of order three; no other child orders or quotient orders inferred.','labelled_prime_quotients':labelled,'triple_profile_classification':triple_profiles,'representative_certificates':records,'independent_literal_controls':literal,'source_sha256':{str(p.relative_to(root.parent)):hashlib.sha256(p.read_bytes()).hexdigest() for p in (root/'PROOF.md',root/'verify.py',root/'literal.cpp',root.parent/'marked_pair_quotient'/'verify.py')}},indent=2,sort_keys=True))
if __name__=='__main__':main()
