"""Independent capped grammar fixed point. Small states use literal covers;
large C3 values use separately generated quotient words, no producer imports.
"""
from kernel import *
from itertools import combinations_with_replacement, product
from functools import lru_cache
import json,hashlib,pathlib
B=127;N=7;cap=lambda x:min(x,B+1)
CYCLE=(2,4,1)

@lru_cache(None)
def words(counts):
 # Generate the only3 possible directed C3 words of prescribed total length.
 total=sum(counts);result=0
 for start in range(3):
  c=[0,0,0]
  for k in range(total):c[(start+k)%3]+=1
  result+=tuple(c)==counts
 return result

def cyclic_h(profiles):
 value=0
 for counts in product(*[range(1,len(p))for p in profiles]):
  w=words(counts)
  if w:
   value+=w*prod(factorial(c)*p[c] for p,c in zip(profiles,counts))
 return value

def prod(xs):
 v=1
 for x in xs:v*=x
 return v

def state(a):
 if len(a)<N:return ('s',len(a),tuple(cap(v) for v in path_cover_literal(a)))
 return ('L',cap(h_recursive(a)))

def val(s):return s[2][1] if s[0]=='s' else s[1]

def closure(seed_order):
 witnesses={}
 for n in range(1,seed_order+1):
  for mask in range(1<<(n*(n-1)//2)):
   a=tournament(n,mask);witnesses.setdefault(state(a),a)
 rounds=0
 while True:
  old=list(witnesses);new={};rounds+=1
  for t in combinations_with_replacement(old,2):
   if all(s[0]=='s' for s in t) and sum(s[1]for s in t)<N:
    a=ordinal(*(witnesses[s]for s in t));new.setdefault(state(a),a)
   else:new.setdefault(('L',cap(prod(map(val,t)))),None)
  for t in combinations_with_replacement(old,3):
   if any(s[0]=='L' for s in t):new.setdefault(('L',B+1),None)
   elif sum(s[1]for s in t)<N:
    a=substitute(CYCLE,[witnesses[s] for s in t]);new.setdefault(state(a),a)
   else:new.setdefault(('L',cap(cyclic_h(tuple(s[2]for s in t)))),None)
  new={s:a for s,a in new.items() if s not in witnesses}
  if not new:break
  witnesses.update(new)
 return {'seed_order':seed_order,'states':len(witnesses),'rounds_including_terminal_check':rounds,'attained':sorted({val(s)for s in witnesses if val(s)<=B}),'small_by_order':{n:sum(s[0]=='s'and s[1]==n for s in witnesses)for n in range(1,N)}}

ans=[closure(1),closure(5)]
require(ans[0]['states']==58 and ans[1]['states']==65,'closure state counts')
require(13 not in ans[0]['attained'] and 13 in ans[1]['attained'] and 19 not in ans[1]['attained'],'grammar holes')
require(sorted(set(ans[1]['attained'])-set(ans[0]['attained']))==[13,39,117],'seed expansion exact delta')
a=(18,8,3,5,14);validate(a);require(h_literal(a)==13,'outside grammar actualH13')
print(json.dumps({'status':'PASS','B':B,'N':N,'closures':ans,'script_sha256':hashlib.sha256(pathlib.Path(__file__).read_bytes()).hexdigest()},sort_keys=True,indent=2))
