"""Independent physical-state and backward-query replay of reversal cubes."""
from pathlib import Path
from functools import lru_cache
from itertools import combinations
import hashlib,json,subprocess,sys
HERE=Path(__file__).resolve().parent;SOURCE=HERE.parent/'unequal_reversal_cube'
def need(x,msg):
 if not x:raise RuntimeError(msg)
def opposite(R):return tuple(((1<<len(R))-1)^r^(1<<i) for i,r in enumerate(R))
def build(blocks,theta,bits):
 active=[i for i,B in enumerate(blocks) if len(B)>1];offsets=[0]
 for B in blocks:offsets.append(offsets[-1]+len(B))
 R=[0]*offsets[-1]
 for i,B in enumerate(blocks):
  if i in active and bits>>active.index(i)&1:B=opposite(B)
  j=(i+(-1 if theta else 1))%3
  for u,row in enumerate(B):R[offsets[i]+u]=(row<<offsets[i])|(((1<<len(blocks[j]))-1)<<offsets[j])
 return tuple(R)
def attach(R,a,c):
 n=len(R);U=list(R)+[1<<(n+1),0]
 for u in range(n):
  for mask,port in ((a,n),(c,n+1)):
   if mask>>u&1:U[u]|=1<<port
   else:U[port]|=1<<u
 return tuple(U)
@lru_cache(None)
def hp(R):
 N=(1<<len(R))-1
 @lru_cache(None)
 def walk(left,u):
  if not left:return 1
  nxt=left&R[u];out=0
  while nxt:
   bit=nxt&-nxt;nxt-=bit;out+=walk(left^bit,bit.bit_length()-1)
  return out
 return sum(walk(N^(1<<u),u) for u in range(len(R)))
def main():
 pins={n:hashlib.sha256((SOURCE/n).read_bytes()).hexdigest() for n in ('PROOF.md','verify.py')}
 saved=json.loads((SOURCE/'normal.json').read_text());states=queries=decoder_checks=0
 for case in saved['cases']:
  blocks=case['block_rows'];r=sum(len(B)>1 for B in blocks);need(case['states']==2**(1+r),'nontrivial state count');need(case['adaptive_query_bound']==3+r and case['fixed_query_count']==3+2*r,'literal query bounds')
  signatures=set();root={0:set(),1:set()};records={}
  for rec in case['records']:
   theta,bits=rec['theta'],rec['bits'];R=build(blocks,theta,bits);need(R==tuple(rec['rows']),'literal fixed physical names and orientations')
   vals=tuple(hp(attach(R,*q)) for q in case['queries']);need(vals==tuple(rec['fixed_values']),'backward-DP complete fixed signature');need(vals not in signatures,'fixed signature collision');signatures.add(vals)
   root[theta].add(vals[:3]);records[theta,bits]=vals;states+=1;queries+=len(vals)
  need(len(records)==2**(1+r) and all(len(root[t])==1 for t in range(2)),'exhaustive cube and root independence');need(sum(next(iter(root[0])))<sum(next(iter(root[1]))),'root sum orientation')
  for d in case['child_decoders']:
   need((d['values'][0]-d['values'][1])%4==2,'decoder residue two');idx=case['queries'].index(d['query'])
   for bits in range(1<<r):need(records[d['theta'],bits][idx]==d['values'][(bits>>d['bit_index'])&1],'decoder independence');decoder_checks+=1
 normal=subprocess.check_output([sys.executable,str(SOURCE/'verify.py')]);optimized=subprocess.check_output([sys.executable,'-O',str(SOURCE/'verify.py')]);need(normal==optimized and json.loads(normal)==saved,'producer replays and saved data agree')
 for name,data in [('producer.normal.stdout',normal),('producer.optimized.stdout',optimized)]: (HERE/name).write_bytes(data)
 need(pins=={n:hashlib.sha256((SOURCE/n).read_bytes()).hexdigest() for n in pins},'source freeze')
 result={'status':'PASS','independent_states':states,'independent_fixed_query_values':queries,'independent_decoder_bit_checks':decoder_checks,'producer_normal_optimized_equal':True,'source_sha256':pins,'audit_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
 (HERE/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,sort_keys=True))
if __name__=='__main__':main()
