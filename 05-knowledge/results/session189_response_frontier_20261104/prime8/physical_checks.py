from pathlib import Path
from itertools import product
import json,subprocess,math
ROOT=Path(__file__).resolve().parent

def need(c,m):
 if not c:raise RuntimeError(m)
def op(T):return [((1<<len(T))-1)^(1<<i)^r for i,r in enumerate(T)]
def extend(T,a,c):
 n=len(T);return [r|(((a>>i)&1)<<n)|(((c>>i)&1)<<(n+1)) for i,r in enumerate(T)]+[((1<<n)-1)^a|(1<<(n+1)),((1<<n)-1)^c]
def tt(n):return [sum(1<<j for j in range(i+1,n)) for i in range(n)]
def blow(Q,blocks):
 T=[];offsets=[]
 for B in blocks:offsets.append(len(T));T.extend([0]*len(B))
 for i,B in enumerate(blocks):
  for u,row in enumerate(B):T[offsets[i]+u]=(row<<offsets[i])|sum(((1<<len(C))-1)<<offsets[j] for j,C in enumerate(blocks) if Q[i]>>j&1)
 return T,offsets

def physical(Q,blocks,a,c,ambient=False):
 T,ofs=blow(Q,blocks);aa=sum(((1<<len(B))-1)<<ofs[i] for i,B in enumerate(blocks) if a>>i&1);cc=sum(((1<<len(B))-1)<<ofs[i] for i,B in enumerate(blocks) if c>>i&1)
 if ambient:
  n=len(T);T=[row|(1<<n) for row in T]+[1<<(n+1),(1<<n)-1];aa|=1<<(n+1);cc|=1<<(n+1)
 return extend(T,aa,cc),aa,cc

def subs(a):
 s=a
 while True:
  yield s
  if not s:break
  s=(s-1)&a

def main():
 A=json.loads((ROOT/'anchor_exploration.json').read_text());need(len(A)==2,'two obstruction representatives');half=[m for m in range(256) if m.bit_count()<=2];qs=list(product(half,repeat=2));cases=[];where=[]
 def add(T,meta):need(len(T)<=20,'20! integer bound');where.append(meta);cases.append(T)
 for ri,entry in enumerate(A):
  for side in range(2):
   Q=entry['rows'] if side==0 else op(entry['rows'])
   for a,c in qs:add(extend(Q,a,c),dict(kind='named_basis',ri=ri,side=side,query=[a,c]))
 probes=[]
 for ri,entry in enumerate(A):
  for pattern in entry['patterns']:
   V=pattern['V'];first=next((i for i,v in enumerate(pattern['first']) if v==2),None)
   if first is not None:q=(0,1<<first);mod=4;gap=2
   else:
    second=next((i for i,v in enumerate(pattern['second']) if v==4),None);need(second is not None,'all512 patterns covered');q=(1<<second,0);mod=8;gap=4
   blocks=[tt(2 if V>>i&1 else 1) for i in range(8)];record=dict(ri=ri,V=V,query=q,modulus=mod,predicted=gap,physical_queries=[],augmented_rows=[],child_sets=[])
   offset=0
   for block in blocks:record['child_sets'].append(list(range(offset,offset+len(block))));offset+=len(block)
   for side in range(2):
    Q=entry['rows'] if side==0 else op(entry['rows']);T,a,c=physical(Q,blocks,*q);record['physical_queries'].append([a,c]);record['augmented_rows'].append(T);record['port_names']=[len(T)-2,len(T)-1];add(T,dict(kind='anchor',record=len(probes),side=side))
   probes.append(record)
 high=[(255,255),(255,0),(0,255),(85,170),(170,85),(199,43),(37,212),(172,97)]
 for ri,entry in enumerate(A):
  for side in range(2):
   Q=entry['rows'] if side==0 else op(entry['rows'])
   for q in high:add(extend(Q,*q),dict(kind='named_high',ri=ri,side=side,query=q))
 ambient=[]
 for ri,entry in enumerate(A):
  for V in [0,37,126,255]:
   record=next(p for p in probes if p['ri']==ri and p['V']==V);blocks=[tt(2 if V>>i&1 else 1) for i in range(8)];meta=dict(record);meta['physical_queries']=[];meta['augmented_rows']=[]
   for side in range(2):
    Q=entry['rows'] if side==0 else op(entry['rows']);T,a,c=physical(Q,blocks,*record['query'],ambient=True);meta['physical_queries'].append([a,c]);meta['augmented_rows'].append(T);meta['port_names']=[len(T)-2,len(T)-1];meta['exterior_names']=[len(T)-4,len(T)-3];add(T,dict(kind='ambient',record=len(ambient),side=side))
   ambient.append(meta)
 inp=str(len(cases))+'\n'+''.join(str(len(T))+'\n'+' '.join(map(str,T))+'\n' for T in cases)
 subprocess.run(['g++','-std=c++17','-O3','-Wall','-Wextra','-Werror',str(ROOT/'physical_h.cpp'),'-o',str(ROOT/'physical_h')],check=True)
 vals=list(map(int,subprocess.check_output([str(ROOT/'physical_h')],input=inp.encode()).split()));need(len(vals)==len(cases),'all literal counts')
 named={(ri,side):{} for ri in range(2) for side in range(2)};highvalues=[]
 for meta,value in zip(where,vals):
  kind=meta['kind'];need(value%2==1,'physical Hamilton parity')
  if kind=='named_basis':named[meta['ri'],meta['side']][tuple(meta['query'])]=value
  elif kind=='anchor':probes[meta['record']].setdefault('values',[0,0])[meta['side']]=value
  elif kind=='ambient':ambient[meta['record']].setdefault('values',[0,0])[meta['side']]=value
  else:highvalues.append((meta,value))
 records=[];full={}
 for ri in range(2):
  sides=[]
  for side in range(2):
   v=named[ri,side];coef={(a,c):sum((-1)**((a^x).bit_count()+(c^y).bit_count())*v[x,y] for x in subs(a) for y in subs(c)) for a,c in qs};f=[0]*65536
   for (a,c),z in coef.items():f[a+(c<<8)]=z
   for bit in range(16):
    b=1<<bit
    for mask in range(65536):
     if mask&b:f[mask]+=f[mask^b]
   need(all(0<x<=math.factorial(10) and x%2 for x in f),'full reconstructed cube positivity and bound')
   need(all(f[a+(c<<8)]==z for (a,c),z in v.items()),'full cube basis reconstruction')
   full[ri,side]=f;sides.append(dict(values=[v[q] for q in qs],coefficients=[coef[q] for q in qs],complete_query_values_a_low_c_high=f))
  gap=[x-y for x,y in zip(full[ri,0],full[ri,1])];first=next(k for k,z in enumerate(gap) if z);records.append(dict(rows=A[ri]['rows'],sides=sides,first_physical_separator=[first&255,first>>8,gap[first]]))
 for meta,value in highvalues:
  a,c=meta['query'];need(full[meta['ri'],meta['side']][a+(c<<8)]==value,'nonbasis physical complete bank check')
 for p in probes+ambient:need((p['values'][0]-p['values'][1])%p['modulus']==p['predicted'],'actual anchor or ambient residue')
 result=dict(status='PASS',basis_queries=qs,records=records,anchor_probes=probes,ambient_probes=ambient,nonbasis_literal_values=[dict(**m,value=v) for m,v in highvalues],physical_cases=len(cases),physical_order_max=max(map(len,cases)))
 (ROOT/'physical_records.json').write_text(json.dumps(result,separators=(',',':'))+'\n');print(json.dumps(dict(status='PASS',physical_cases=len(cases),physical_order_max=max(map(len,cases)),full_named_values=262144,anchor_patterns=len(probes),ambient_probes=len(ambient))))
if __name__=='__main__':main()
