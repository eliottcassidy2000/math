from pathlib import Path
from itertools import product
import json,subprocess,hashlib,math
ROOT=Path(__file__).resolve().parent

def need(c,m):
 if not c:raise RuntimeError(m)
def subsets(x):
 y=x
 while True:
  yield y
  if not y:break
  y=(y-1)&x

def prime(T):
 for m in range(1,127):
  if m.bit_count()<2:continue
  if all((T[v]&m) in (0,m) for v in range(7) if not m>>v&1):return False
 return True

def run_bank(rows,queries):
 data=f'{len(rows)} {len(queries)}\n'+''.join(f'{a} {c}\n' for a,c in queries)+''.join(' '.join(map(str,r))+'\n' for r in rows)
 result=subprocess.check_output([str(ROOT/'even_bank')],input=data.encode()).decode()
 ls=[list(map(int,line.split())) for line in result.splitlines()]
 need(len(ls)==len(rows)*len(queries)*2 and all(len(x)==110 for x in ls),'complete 110-channel bank')
 need(all(all(0<=x<16 for x in v[:5]) and all(x in (0,1) for x in v[5:]) for v in ls),'residue domains')
 return result,ls

def main():
 rows=json.loads((ROOT/'rows.json').read_text());need(len(rows)==4 and all(prime(r) for r in rows),'four prime quotients')
 for T in rows:
  need(all(not(T[i]>>i&1) and all(((T[i]>>j)&1)+((T[j]>>i)&1)==1 for j in range(7) if j!=i) for i in range(7)),'tournament rows')
 need(rows[1]==[sum(1<<((i+d)%7) for d in (1,2,3)) for i in range(7)],'carousel physical labels')
 half=[m for m in range(128) if m.bit_count()<=2];qs=list(product(half,repeat=2));qindex={q:i for i,q in enumerate(qs)};need(len(qs)==841,'complete bidegree basis')
 subprocess.run(['g++','-std=c++17','-O2','-Wall','-Wextra','-Werror',str(ROOT/'even_bank.cpp'),'-o',str(ROOT/'even_bank')],check=True)
 raw,bank=run_bank(rows,qs);(ROOT/'even_queries.json').write_text(json.dumps(qs));(ROOT/'even_bank.txt').write_text(raw)
 coefficients=[];tables=[];moduli=[16]*5+[2]*105
 for ri in range(4):
  sides=[]
  for side in range(2):
   coeff={}
   for a,c in qs:
    value=bank[2*(ri*841+qindex[a,c])+side]
    out=[]
    for channel,mod in enumerate(moduli):
     z=sum((-1)**((a^x).bit_count()+(c^y).bit_count())*bank[2*(ri*841+qindex[x,y])+side][channel] for x in subsets(a) for y in subsets(c))%mod
     out.append(z)
    coeff[a,c]=out
   for a,c in qs:
    rebuilt=[sum(coeff[x,y][h] for x in subsets(a) for y in subsets(c))%moduli[h] for h in range(110)]
    need(rebuilt==bank[2*(ri*841+qindex[a,c])+side],'complete basis reconstruction')
   sides.append(coeff)
  tables.append(sides)
  coefficients.append(dict(rows=rows[ri],bank=[[a,c,sides[0][a,c],sides[1][a,c]] for a,c in qs]))
 (ROOT/'coefficients.json').write_text(json.dumps(dict(channel_order=['baseline_q'+str(h) for h in range(1,6)]+[f'D_k{k}_child{i}_q{h}' for k in (2,3,5) for i in range(7) for h in range(1,6)],records=coefficients),separators=(',',':'))+'\n')
 high=[(127,127),(127,0),(0,127),(85,42),(42,85),(109,91),(73,54),(31,103)]
 _,controls=run_bank(rows,high)
 for ri in range(4):
  for qi,(a,c) in enumerate(high):
   for side in range(2):
    cf=tables[ri][side];parts=[v for (x,y),v in cf.items() if x&a==x and y&c==y];rebuilt=[sum(v[h] for v in parts)%moduli[h] for h in range(110)]
    need(rebuilt==controls[2*(ri*len(high)+qi)+side],'nonbasis full-query check, all channels')
 specifications=[(0,(1,0),8,[]),(1,(0,3),0,[3,5]),(1,(4,68),8,[3,5]),(2,(0,3),0,[0,1]),(2,(4,68),8,[0,1]),(3,(3,1),8,[])]
 certificates=[]
 for ri,q,b,active in specifications:
  at=2*(ri*841+qindex[q]);source,target=bank[at:at+2];difference=[(x-y)%m for x,y,m in zip(source,target,moduli)]
  expected=[b,0,0,0,0]+[int(k==2 and i in active and h==1) for k in (2,3,5) for i in range(7) for h in range(1,6)]
  need(difference==expected,'complete all21-bit all-cover identity')
  certificates.append(dict(quotient_index=ri,query=q,baseline_contrast=difference[:5],hamilton_q2_bit_indices=active,all_higher_coordinate_contrasts_zero=True))
 independent=subprocess.check_output(['python3',str(ROOT/'independent.py')]).decode();(ROOT/'independent.json').write_text(independent);ind=json.loads(independent);need(ind['status']=='PASS' and ind['physical_cases']==1524 and ind['physical_order_max']==13,'independent literal physical interpolation')
 data=['rows.json','even_queries.json','even_bank.txt','coefficients.json','independent.json'];sources=['PROOF.md','verify.py','even_bank.cpp','physical.cpp','independent.py']
 def receipt(files):return {f:dict(bytes=(ROOT/f).stat().st_size,sha256=hashlib.sha256((ROOT/f).read_bytes()).hexdigest()) for f in files}
 print(json.dumps(dict(status='PASS',scope='Four explicitly named seven-quotient converses, arbitrary even children of order at least four, complete all-cover modulo16 identities and every fixed ambient. No classification completeness claim in this packet.',basis_queries=841,nonbasis_complete_channel_checks=64,channels_per_side=110,certificates=certificates,independent_physical_cases=1524,source_pins=receipt(sources),data_pins=receipt(data)),sort_keys=True,indent=2))
if __name__=='__main__':main()
