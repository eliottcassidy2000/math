#!/usr/bin/env python3
import json,itertools,hashlib
from pathlib import Path
HERE=Path(__file__).resolve().parent

def need(x,msg):
 if not x:raise RuntimeError(msg)
def tt(n,p=None):
 p=list(range(n))if p is None else p
 return [sum(1<<v for v in p[p.index(u)+1:])for u in range(n)]
def op(r):
 full=(1<<len(r))-1
 return [full^s^(1<<i)for i,s in enumerate(r)]
def sub(Q,D):
 blocks=[];pos=0
 for d in D:blocks.append(list(range(pos,pos+len(d))));pos+=len(d)
 out=[0]*pos
 for i,d in enumerate(D):
  for u,row in enumerate(d):
   for v in range(len(d)):
    if row>>v&1:out[blocks[i][u]]|=1<<blocks[i][v]
  for j,B in enumerate(blocks):
   if Q[i]>>j&1:
    mask=sum(1<<v for v in B)
    for u in blocks[i]:out[u]|=mask
 return out,blocks

def dp(r):
 n=len(r);N=1<<n;end=[[0]*n for _ in range(N)];h=[0]*N;h[0]=1
 for i in range(n):end[1<<i][i]=1
 for s in range(1,N):
  for i,z in enumerate(end[s]):
   if z:
    bits=r[i]&~s
    while bits:
     b=bits&-bits;bits-=b;end[s|b][b.bit_length()-1]+=z
  h[s]=sum(end[s])
 return h,end

def actual(r,a,c):
 n=len(r);z=list(r)+[1<<(n+1),0]
 for i in range(n):
  if a>>i&1:z[i]|=1<<n
  else:z[n]|=1<<i
  if c>>i&1:z[i]|=1<<(n+1)
  else:z[n+1]|=1<<i
 return dp(z)[0][-1]

def prepare():
 cycle=[2,4,1];prime5=[6,12,24,17,3]
 definitions=[('singleton',[ [0],cycle,tt(2) ],[[0],cycle,tt(2)]),('cyclic_triple',[cycle,tt(3),tt(2)],[op(cycle),tt(3,[1,2,0]),tt(2)]),('transitive_five',[tt(5),tt(2),tt(2)],[op(tt(5)),tt(2),tt(2)]),('prime_five',[prime5,tt(2),tt(2)],[op(prime5),tt(2),tt(2)])]
 controls=[];graphs=[]
 for name,D,E in definitions:
  T,blocks=sub(cycle,D);U,_=sub(op(cycle),E);V=op(T);core_n=len(T)
  src,_=sub(cycle,[T,[0],[0]]);tgt,_=sub(cycle,[U,[0],[0]]);baseline,_=sub(cycle,[V,[0],[0]])
  # The exterior name in block2 points toward core block0.
  outerbits=1<<(core_n+1);A=blocks[0];v=A[0];compressed=[[u]for u in A]+blocks[1:];others=[B for B in compressed if B!=[v]]
  d=sum(bool(T[v]>>B[0]&1)for B in others);battery=[]
  for B in others:
   a=sum(1<<i for i in B)
   if T[v]>>B[0]&1:a|=1<<v
   battery.append([a|outerbits,(1<<v)|outerbits])
  if d%2:battery.append([(1<<v)|outerbits,(1<<v)|outerbits])
  if (len(compressed)-1+d)%2:battery.append([outerbits,(1<<v)|outerbits])
  battery.append([(1<<v)|outerbits,outerbits])
  need(len(set(tuple(z)for z in battery))==len(battery),'duplicate battery')
  for dd,ee in zip(D,E):need(dp(dd)[0]==dp(ee)[0],'hereditary profile equality not witnessed')
  h,ends=dp(T);_,starts=dp(op(T));bv=(ends[-1][v]+starts[-1][v])%2
  records=[]
  for a,c in battery:
   vals=[actual(r,a,c)for r in (src,tgt,baseline)]
   need(vals[1]==vals[2],'actual target/restored converse mismatch')
   records.append({'a_mask':a,'c_mask':c,'source':vals[0],'target':vals[1],'whole_core_converse':vals[2],'difference':vals[0]-vals[1],'residue_mod4':(vals[0]-vals[1])%4})
  need(any(z['residue_mod4']==2 for z in records),'anchor not separating')
  if bv:need(records[-1]['residue_mod4']==2,'endpoint parity branch')
  else:need(sum(z['difference']for z in records[:-1])%4==2,'compressed battery parity sum')
  ids=[]
  for kind,r in [('source',src),('target',tgt),('whole_core_converse',baseline)]:
   ids.append(len(graphs));graphs.append({'id':len(graphs),'case':name,'kind':kind,'rows':r})
  controls.append({'name':name,'core_rows':T,'core_target_rows':U,'core_converse_rows':V,'child_blocks':blocks,'anchor':v,'compressed_blocks':compressed,'ambient_quotient_rows':cycle,'outside_one_mask':outerbits,'b_anchor':bv,'queries':records,'graph_ids':ids})
 (HERE/'controls.json').write_text(json.dumps({'cases':controls,'graphs':graphs},indent=2)+'\n')
 with (HERE/'graphs.txt').open('w')as f:
  f.write(str(len(graphs))+'\n')
  for z in graphs:f.write(str(len(z['rows']))+' '+' '.join(map(str,z['rows']))+'\n')
 return controls,graphs

def check_bank():
 z=json.loads((HERE/'controls.json').read_text());bank=json.loads((HERE/'bank.json').read_text());need(len(bank)==len(z['graphs']),'bank size')
 for B,G in zip(bank,z['graphs']):
  need(B['rows']==G['rows'],'physical bank rows');m=B['masks'];L=len(m);co=B['coefficients'];need(len(co)==L*L,'incomplete named polynomial')
  for case in z['cases']:
   if G['id']in case['graph_ids']:
    for q in case['queries']:
     val=sum(co[i*L+j]for i,a in enumerate(m)if a&q['a_mask']==a for j,c in enumerate(m)if c&q['c_mask']==c)
     need(val==q[G['kind']],'full named polynomial fails actual anchor query')
 return {'cases':len(z['cases']),'graphs':len(bank),'physical_query_triples':sum(len(c['queries'])for c in z['cases']),'complete_named_coordinates':sum(len(b['coefficients'])for b in bank),'old_orders':sorted({len(b['rows'])for b in bank}),'all_cubic_quartic_zeros_retained':True,'bank_sha256':hashlib.sha256((HERE/'bank.json').read_bytes()).hexdigest()}
if __name__=='__main__':
 import sys
 if len(sys.argv)>1 and sys.argv[1]=='prepare':
  a,b=prepare();print(json.dumps({'prepared_cases':len(a),'graphs':len(b)}))
 else:print(json.dumps(check_bank(),sort_keys=True,indent=2))
