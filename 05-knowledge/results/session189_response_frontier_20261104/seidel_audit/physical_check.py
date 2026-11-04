from pathlib import Path
from itertools import product
import json,subprocess
from audit import response,mod,need
D=Path(__file__).resolve().parent

def tt(n):return [sum(1<<j for j in range(i+1,n)) for i in range(n)]
def op(Q):return [((1<<len(Q))-1)^(1<<i)^r for i,r in enumerate(Q)]
def blow(Q,size):
 m=len(Q);n=m*size;T=[0]*n;blocks=[list(range(i*size,(i+1)*size)) for i in range(m)]
 for i in range(m):
  for u in blocks[i]:T[u]=sum(1<<v for v in blocks[i] if u<v)|sum(sum(1<<v for v in blocks[j]) for j in range(m) if Q[i]>>j&1)
 return T,blocks

def extend(T,a,c):
 n=len(T);return [r|(((a>>i)&1)<<n)|(((c>>i)&1)<<(n+1)) for i,r in enumerate(T)]+[(((1<<n)-1)^a)|(1<<(n+1)),((1<<n)-1)^c]
def subs(a):
 b=a
 while True:
  yield b
  if not b:break
  b=(b-1)&a

def main():
 configs=[dict(name='C3_pair_children',Q=[2,4,1],size=2,modulus=4,full_named=True),dict(name='edge_TT4_children',Q=[2,0],size=4,modulus=16,full_named=True),dict(name='C3_TT4_children',Q=[2,4,1],size=4,modulus=16,full_named=False)]
 cases=[];meta=[];details=[]
 for ci,conf in enumerate(configs):
  record=dict(conf,physical_old_rows=[],physical_child_sets=None,named_basis_queries=None,named_basis_covers=[None,None],coarse_queries=[],coarse_physical_queries=[],coarse_covers=[None,None],complete_named_hamilton_values=[None,None],complete_named_coefficients=[None,None]);details.append(record)
  for side in range(2):
   Q=conf['Q'] if side==0 else op(conf['Q']);T,blocks=blow(Q,conf['size']);n=len(T);record['physical_old_rows'].append(T);record['physical_child_sets']=blocks;record['physical_port_names']=[n,n+1]
   if conf['full_named']:
    half=[x for x in range(1<<n) if x.bit_count()<=2];qs=list(product(half,repeat=2));record['named_basis_queries']=qs
    for a,c in qs:cases.append(extend(T,a,c));meta.append((ci,side,'named',a,c))
   for a,c in product(range(1<<len(Q)),repeat=2):
    aa=sum(1<<v for i,b in enumerate(blocks) if a>>i&1 for v in b);cc=sum(1<<v for i,b in enumerate(blocks) if c>>i&1 for v in b);cases.append(extend(T,aa,cc));meta.append((ci,side,'coarse',a,c))
    if side==0:record['coarse_queries'].append([a,c]);record['coarse_physical_queries'].append([aa,cc])
 subprocess.run(['g++','-std=c++17','-O3','-Wall','-Wextra','-Werror',str(D/'physical.cpp'),'-o',str(D/'physical')],check=True);inp=str(len(cases))+'\n'+''.join(str(len(T))+'\n'+' '.join(map(str,T))+'\n' for T in cases);out=subprocess.check_output([str(D/'physical')],input=inp.encode()).decode();rows=[list(map(int,l.split())) for l in out.splitlines()];need(len(rows)==len(cases),'all literal physical cases');(D/'literal_cover_values.txt').write_text(out)
 named={};coarse={}
 for (ci,side,kind,a,c),values in zip(meta,rows):
  target=named if kind=='named' else coarse;target.setdefault((ci,side),{})[a,c]=values
 for ci,conf in enumerate(configs):
  rec=details[ci]
  for side in range(2):
   Q=conf['Q'] if side==0 else op(conf['Q']);observed=coarse[ci,side];M=conf['modulus']
   for (a,c),v in observed.items():need([x%M for x in v]==[mod(x)%M for x in response(Q,a,c)],'physical first-five covers versus exact geometric Gram')
   rec['coarse_covers'][side]=[observed[tuple(q)] for q in rec['coarse_queries']]
   if conf['full_named']:
    n=len(rec['physical_old_rows'][side]);basis=named[ci,side];qs=rec['named_basis_queries'];coefs=[];f=[0]*(1<<(2*n))
    for a,c in qs:
     coeff=[sum((-1)**((a^x).bit_count()+(c^y).bit_count())*basis[x,y][h] for x in subs(a) for y in subs(c)) for h in range(5)];coefs.append(coeff);f[a+(c<<n)]=coeff[0]
    for bit in range(2*n):
     b=1<<bit
     for mask in range(1<<(2*n)):
      if mask&b:f[mask]+=f[mask^b]
    need(all(v>0 and v%2 for v in f),'full physical Hamilton parity')
    for coarse_query,physical_query in zip(rec['coarse_queries'],rec['coarse_physical_queries']):
     a,c=physical_query;need(f[a+(c<<n)]==observed[tuple(coarse_query)][0],'full named cube cross-check at every coarse physical mask')
    rec['named_basis_covers'][side]=[basis[tuple(q)] for q in qs];rec['complete_named_coefficients'][side]=coefs;rec['complete_named_hamilton_values'][side]=f
  if not conf['full_named']:rec.pop('complete_named_hamilton_values');rec.pop('complete_named_coefficients');rec.pop('named_basis_covers');rec.pop('named_basis_queries')
 (D/'physical_banks.json').write_text(json.dumps(details,separators=(',',':'))+'\n');print(json.dumps(dict(status='PASS',physical_cases=len(cases),maximum_augmented_order=max(map(len,cases)),complete_named_hamilton_values=sum(len(v) for r in details if 'complete_named_hamilton_values' in r for v in r['complete_named_hamilton_values']),physical_coarse_queries=sum(len(r['coarse_queries'])*2 for r in details),cover_coordinates=5),sort_keys=True,indent=2))
if __name__=='__main__':main()
