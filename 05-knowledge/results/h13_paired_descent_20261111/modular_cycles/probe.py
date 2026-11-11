"""Native paired-extension probes; inherited compiler supplies literal graphs."""
from functools import lru_cache
from collections import defaultdict
from pathlib import Path
import importlib.util,json,hashlib,itertools
ROOT=Path(__file__).resolve().parents[4]
s=importlib.util.spec_from_file_location('c',ROOT/'05-knowledge/results/session251_creative_completion_20261111/ternary/compiler.py')
c=importlib.util.module_from_spec(s);s.loader.exec_module(c)

def state(rows,r):
 n=len(rows);full=(1<<n)-1;d=[[0]*n for _ in range(full+1)];z=[[0]*n for _ in range(full+1)]
 for v in range(n):d[1<<v][v]=1
 z[1<<r][r]=1
 for mask in range(1,full+1):
  for v in range(n):
   if not d[mask][v] and not z[mask][v]:continue
   av=rows[v]&~mask&full
   while av:
    b=av&-av;av-=b;u=b.bit_length()-1;d[mask|b][u]+=d[mask][v];z[mask|b][u]+=z[mask][v]
 h=[sum(a) for a in d];h[0]=1
 return h[-1],d[-1][r],sum(sum(z[m])*h[full^m] for m in range(1,full+1)),h[full^(1<<r)]

def add(rows,o):
 n=len(rows);r=list(rows)+[o]
 for v in range(n):
  if not o>>v&1:r[v]|=1<<n
 return r

def paired_births(twins,root,delta=(2,0,0,0)):
 n=len(twins[0]);arr=[[state(add(t,o),root) for o in range(1<<n)] for t in twins]
 lookup=defaultdict(list)
 for o,a in enumerate(arr[1]):lookup[tuple(x-y for x,y in zip(a,delta))].append(o)
 return [(o,o2,a) for o,a in enumerate(arr[0]) for o2 in lookup[a]]

def word_state(st,B):
 for j in range(B.bit_length()):
  h,E,U,d=st
  st=(h+U-E,d,2*U-d,d) if B>>j&1 else (h,0,2*U,d)
 return st

def letter_state(st,letter):
 h,E,U,d=st
 if letter=='O':return h,0,2*U,d
 if letter=='Z':return h+U-E,d,2*U-d,d
 if letter=='J':return h+U-E,h-E,h+U,h
 raise ValueError(letter)

def physical_word(rows,root,word):
 for letter in word:
  rows=add(rows,0 if letter=='O' else 1<<root)
  if letter=='J':root=len(rows)-1
 return rows,root

def deletion_table(core,cut):
 result=[]
 for v in range(len(core)):
  vs=[u for u in range(len(core)) if u!=v]
  rows=[sum(((core[u]>>w)&1)<<j for j,w in enumerate(vs)) for u in vs]
  newcut=sum(((cut>>w)&1)<<j for j,w in enumerate(vs))
  d,kp,km,ip,im,op,om,jp,jm=c.signature(rows,newcut)
  result.append((d,op,-om,jp+jm))
 return result

if __name__=='__main__':
 Q=[18,8,3,5,14]; C0=add(Q,14);C1=add(Q,8)
 report={'status':'FINITE-EXACT new-operation controls; all-height claims proved in PROOF.md; Chan-Pak OPEN',
         'base':[C0,C1], 'signatures':[c.signature(C0,38),c.signature(C1,35)],
         'deletion_tables':[deletion_table(C0,38),deletion_table(C1,35)],'controls':[]}
 for B in range(16):
  ts=[c.construct(C,cut,0,B) for C,cut in ((C0,38),(C1,35))];r=len(ts[0])-1
  states=[state(t,r) for t in ts]
  expected=25*c.g(B)+122*B+16*(B&1)+48
  if states[0][0]!=expected or states[1][0]!=expected+2*(B&1):raise RuntimeError(('H',B,states))
  if B:
   b=B.bit_length();R=int(f'{B:0{b}b}'[::-1],2)
   side=(25,147*(1<<b)-25*R,25)
   if states[0][1:]!=side or states[1][1:]!=side:raise RuntimeError(('sidecar',B,states))
   for t in ts:
    hs=c.hp_subsets(t);full=(1<<len(t))-1
    if min(hs[full^(1<<v)] for v in range(len(t)))!=25:raise RuntimeError(('deletion floor',B))
  report['controls'].append({'B':B,'states':states})
 ts=[c.construct(C,cut,0,1) for C,cut in ((C0,38),(C1,35))];r=len(ts[0])-1
 found=paired_births(ts,r)
 report['next_pair_births']=[{'mask_lower':a,'mask_upper':b,'lower_state':st,'delta':st[2]-st[0]-st[3]+st[1]} for a,b,st in found]
 old=[c.construct(Q,cut,0,4) for cut in (6,3)]
 report['same_count_source_reset']={'old_rows':old,'new_rows':ts,
  'old_states':[state(t,len(t)-1) for t in old], 'new_states':[state(t,len(t)-1) for t in ts]}
 report['original_pair_births']={}
 report['hidden_endpoint_pair_births']={}
 for B in range(4):
  ori=[c.construct(Q,cut,0,B) for cut in (6,3)]
  report['original_pair_births'][B]=paired_births(ori,len(ori[0])-1)
  report['hidden_endpoint_pair_births'][B]=paired_births(ori,len(ori[0])-1,(0,-2,0,0))
 equal=[c.construct(C,cut,0,0) for C,cut in ((C0,38),(C1,35))]
 report['equal73_pair_births']=paired_births(equal,len(equal[0])-1)
 for B in range(1024):
  a=word_state((211,13,347,13),B);b=word_state((211,25,269,25),B)
  gap=90*B-12*(c.g(B)-1)
  if a[0]-b[0]!=gap or gap<78*B+12*(B&1):raise RuntimeError(('reset descent',B,a,b))
 report['reset_word_controls']=1024
 pair=[letter_state(st,'J') for st in ((31,7,45,13),(33,7,45,13))]
 checks=0
 for length in range(8):
  for word in itertools.product('OZJ',repeat=length):
   a,b=pair;unit=(1,1,1,1)
   for letter in word:
    a=letter_state(a,letter);b=letter_state(b,letter);unit=letter_state(unit,letter)
   if tuple(y-x for x,y in zip(a,b))!=tuple(2*x for x in unit):raise RuntimeError(('difference compiler',word))
   checks+=1
 report['difference_compiler_word_checks']=checks
 literal=[]
 for word in ('J','JOZ','JOJ','JZJ'):
  sts=[]
  for cut in(6,3):
   t=c.construct(Q,cut,0,0);rr=len(t)-1;t,rr=physical_word(t,rr,word)
   sts.append(state(t,rr))
  literal.append({'word':word,'states':sts,'gap':sts[1][0]-sts[0][0]})
 report['literal_J_word_controls']=literal
 report['universe']={'literal_family':'B0..15, both25-core charts; all vertex deletions forB>=1',
 'source_search':'all64/128/256/256 masks on each originalB0/1/2/3 graph, separately, retaining root; full matching fibre product',
 'next_pair_search':'all256 masks on each211/213 graph, retainingroot',
 'equal73_search':'all128 masks on each73graph, retainingroot',
 'word_descent':'all0<=B<1024 from bothactual211states'}
 report['universe']['hidden_endpoint_search']='fullpairedcutsetsat originalB0..3, targetdifference(0,-2,0,0)'
 report['universe']['difference_compiler']='allO/Z/Jwordsoflength0..7; fourliteralJ-prefixprograms from31/33'
 report['compiler_sha256']=hashlib.sha256((ROOT/'05-knowledge/results/session251_creative_completion_20261111/ternary/compiler.py').read_bytes()).hexdigest()
 print(json.dumps({'family_controls':len(report['controls']),'next_pair_births':len(found),'equal73_pair_births':report['equal73_pair_births'],'source_reset':report['same_count_source_reset']['new_states']},indent=2),flush=True)
 Path(__file__).with_name('receipt.json').write_text(json.dumps(report,indent=2)+'\n')
