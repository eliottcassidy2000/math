#!/usr/bin/env python3
"""Physical stationary one-core inverse. No exploratory imports; -O-safe."""
from verify_retained_root_sewing import *
from pathlib import Path
from array import array
from functools import lru_cache

@lru_cache(None)
def F(q):
 require(q>=0,'nonnegative binary index')
 if not q:return 1
 x=q//2
 return F(x)+(0 if q&1 else x+(x&1))

def f_table(limit):
 out=[1]*(limit+1)
 for q in range(1,limit+1):
  x=q//2;out[q]=out[x]+(0 if q&1 else x+(x&1))
 return out

def step(st,letter):
 h,e,R,d=st
 return (h,0,2*R,d) if letter=='O' else (h+R-e,d,2*R-d,d)

def birth(A,letter):
 n=len(A);old=tuple(z|(1<<n)for z in A)
 if letter=='O':return old+(0,)
 return (old[0]^(1<<n),)+old[1:]+(1,)

def prefix_word(q):return '' if not q else ''.join('O' if b=='1' else 'Z' for b in bin(q)[2:])
def suffix_word(t):return '' if not t else ''.join('Z' if b=='1' else 'O' for b in bin(t)[2:][::-1])
def word_graph(A,word):
 for c in word:A=birth(A,c)
 validate(A);return A

def old_state(A):
 z,_=state(A,0);return (z[0],z[1],z[3],z[5])

def core_graph(z,q,t):
 A=word_graph((0,),prefix_word(q));A=sew(A,0,z['rows'],z['root'])
 return word_graph(A,suffix_word(t))

def scalar(sig,q,t,f=F):
 H,E,S,Q,L,C=sig
 return H*(f(q)+f(t)-1)+(L+C*t)*(q+1)-S*((q+1)&1)-t*(Q+H)+(t&1)*(H-E)

def suffix_scalar(st,t):
 h,e,R,d=st
 return h+t*(R-d)+d*(F(t)-1)+(t&1)*(d-e)

def weighted(q):
 b=[(q>>j)&1 for j in range(q.bit_length())];z=0
 for i in range(len(b)):
  for j in range(i):z+=b[i]*(1-b[j])*(1 if i-j==1 else 1<<(i-j-2))
 return 1+2*z

def digest(x):return hashlib.sha256(json.dumps(x,separators=(',',':')).encode()).hexdigest()

def h_flat(A):
 # Flat unsigned table avoids per-cell Python objects for the 19-vertex witness.
 n=len(A);N=1<<n;dp=array('Q',[0])*(N*n)
 for v in range(n):dp[(1<<v)*n+v]=1
 for mask in range(1,N):
  base=mask*n;V=mask
  while V:
   bit=V&-V;V-=bit;v=bit.bit_length()-1;c=dp[base+v]
   if not c:continue
   U=A[v]&~mask
   while U:
    b=U&-U;U-=b;w=b.bit_length()-1;dp[(mask|b)*n+w]+=c
 return sum(dp[(N-1)*n:])

def census(bank,M):
 ft=f_table(2*M);hits=bytearray(M//2+1);chosen={};targets={67,163,173,179,181,247,511,539,541,821,4279,4999,6377}
 for q in range(0,2*M,2):
  n=ft[q]
  if n<=M:
   hits[n//2]=1
   if n in targets:chosen.setdefault(n,{'kind':'plain','q':q})
 plain_count=sum(hits);tests=0;per=[]
 for z in bank:
  sig=z['signature'];H,E,S,Q,L,C=sig
  if H==1:continue
  count=tested=0
  for q in range(M//L):
   h=H*ft[q]-S*((q+1)&1)+L*(q+1);a=C*(q+1)-Q-H
   require(a>0,'strict inverse denominator')
   if h>M:continue
   for t in range((M-h)//a+1):
    tested+=1
    n=h+t*a+H*(ft[t]-1)+(t&1)*(H-E)
    if n<=M:
     require(n>0 and n&1,'physical odd output')
     hits[n//2]=1;count+=1
     if n in targets:chosen.setdefault(n,{'kind':'core','index':z['index'],'q':q,'t':t})
  tests+=tested;per.append({'index':z['index'],'tested':tested,'accepted_pairs':count})
 support=[n for n in range(1,M+1,2)if hits[n//2]]
 return {'bound':M,'nontransitive_signatures':len(per),'plain_values':plain_count,'core_pairs_tested':tests,'odd_outputs':len(support),'missing_odds':[n for n in range(1,M+1,2)if not hits[n//2]],'support_sha256':digest(support),'per_core':per},chosen

def run(M):
 deps={p:hashlib.sha256(Path(__file__).with_name(p).read_bytes()).hexdigest()for p in ['verify.py','verify_retained_root_sewing.py','RETAINED_ROOT_SEWING.md','RETAINED_SEWING_SOURCE_SCOPE.md']}
 letters,graphs,roots=make_bank(5)
 bank=[{'index':k,'signature':sig,'rows':A,'root':a}for k,(sig,(A,a))in enumerate(letters)]
 require((graphs,roots,len(bank),sum(z['signature'][0]>1 for z in bank))==(1099,5405,65,50),'complete literal core inventory')
 pref=0
 for q in range(128):
  A=word_graph((0,),prefix_word(q));expected=(F(q),(q+1)&1,q+1,1)
  require(old_state(A)==expected,'literal prefix full-state');pref+=1
  require(old_state(word_graph((0,),'ZZ'+prefix_word(q)))==expected,'leading zero full-state padding')
 for q in range(4096):require(F(q)==weighted(q),'literal weighted inversion sum')
 reciprocity=0
 for width in range(11):
  for q in range(1<<width):
   rev=int(format(q,'0%db'%width)[::-1],2)if width else 0
   require(F(q)==F((1<<width)-1-rev),'fixed-width reverse complement');reciprocity+=1
 small=[tournament(n,c)for n in range(1,4)for c in range(1<<(n*(n-1)//2))]
 suffixes=0;words=0
 for A in small:
  st=old_state(A)
  for t in range(32):
   C=word_graph(A,suffix_word(t));require(h_subsets(C)[-1]==suffix_scalar(st,t),'literal arbitrary-input suffix');suffixes+=1
  for width in range(7):
   for bits in range(1<<width):
    word=''.join('Z'if bits>>i&1 else'O'for i in range(width));s=st
    for c in word:s=step(s,c)
    require(s[0]==suffix_scalar(st,bits),'full word suffix normalization including trailing O');words+=1
 physical=0
 for z in bank:
  if z['signature'][0]==1:continue
  H,E,S,Q,L,C=z['signature']
  for q in range(4):
   mid=sew(word_graph((0,),prefix_word(q)),0,z['rows'],z['root'])
   expected=(H*F(q)-S*((q+1)&1)+L*(q+1),E,C*(q+1)-Q,H)
   require(old_state(mid)==expected,'actual core middle full-state')
   for t in range(4):
    A=word_graph(mid,suffix_word(t));require(h_subsets(A)[-1]==scalar(z['signature'],q,t),'actual core binary graph');physical+=1
 receipt,chosen=census(bank,M)
 if M==1000000:
  require((receipt['plain_values'],receipt['core_pairs_tested'],receipt['odd_outputs'],receipt['missing_odds'])==(249077,14197377,499997,[7,21,462661]),'complete million-bound inverse')
 elif M==200000:require((receipt['plain_values'],receipt['core_pairs_tested'],receipt['odd_outputs'],receipt['missing_odds'])==(50006,2539905,99998,[7,21]),'complete 200000-bound inverse')
 witnesses=[]
 for target,z in sorted(chosen.items()):
  A=word_graph((0,),prefix_word(z['q']))if z['kind']=='plain'else core_graph(bank[z['index']],z['q'],z['t'])
  require(h_subsets(A)[-1]==target,'selected literal target witness')
  witnesses.append(dict(z,target=target,rows=A,order=len(A)))
 # A literal second nontransitive core repairs the one-core gap.
 z1,z2=bank[17],bank[31];q,t=2,112
 A=word_graph((0,),prefix_word(q));A=sew(A,0,z1['rows'],z1['root']);A=sew(A,0,z2['rows'],z2['root'])
 middle=old_state(A);require(middle==(415,0,4128,15),'two-core physical middle state')
 A=word_graph(A,suffix_word(t));require(len(A)==19 and suffix_scalar(middle,t)==462661 and h_flat(A)==462661,'actual second-core escape')
 escape={'target':462661,'core_indices':[17,31],'q':q,'t':t,'middle_state':middle,'rows':A,'order':len(A)}
 print(json.dumps({'status':'PASS','dependencies':deps,'core_graphs':graphs,'rooted_core_candidates':roots,'inventory':bank,'literal_prefixes':pref,'weighted_indices':4096,'reverse_complements':reciprocity,'actual_suffixes':suffixes,'suffix_words':words,'actual_one_core_graphs':physical,'census':receipt,'witnesses':witnesses,'two_core_escape':escape},sort_keys=True,indent=2))

if __name__=='__main__':
 parser=argparse.ArgumentParser();parser.add_argument('--bound',type=int,default=1000000);args=parser.parse_args();require(args.bound>=1,'positive target bound');run(args.bound)
