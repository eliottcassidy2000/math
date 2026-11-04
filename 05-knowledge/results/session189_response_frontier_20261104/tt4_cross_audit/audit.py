#!/usr/bin/env python3
from pathlib import Path
from itertools import combinations,permutations
import hashlib,json,subprocess,sys,tempfile,shutil
from algebra import *
P=Path(__file__).resolve().parent

def fullbanks(graphs,name):
 text=str(len(graphs))+'\n'+''.join(str(len(r))+' '+' '.join(map(str,r))+'\n'for r in graphs);(P/(name+'.in')).write_text(text);out=subprocess.check_output([str(P/'literal_polynomial')],input=text,text=True);(P/(name+'.json')).write_text(out);return json.loads(out)
def coef(B,p,q,r,i):
 ix={a:j for j,a in enumerate(B['masks'])};return B['coefficients'][ix[(1<<p)|(1<<q)]*len(ix)+ix[(1<<r)|(1<<i)]]
def physical(A,R,into):
 m=len(A);out=A[:]+[r<<m for r in R]
 for i in range(m):
  for j in range(len(R)):
   if into>>j&1:out[m+j]|=1<<i
   else:out[i]|=1<<(m+j)
 return out

def main():
 src=Path(sys.argv[1]).resolve()if len(sys.argv)>1 else P.parent/'session189-tt4-cross';manifest=json.load(open(src/'RELEASE_FILES.json'));pins={f:hashlib.sha256((src/f).read_bytes()).hexdigest()for f in manifest};need(all(pins[f]==v['sha256']for f,v in manifest.items()),'producer frozen pins')
 for source,exe in [('literal_polynomial.cpp','literal_polynomial'),('queries.cpp','queries')]:subprocess.run(['g++','-O3','-std=c++17',str(P/source),'-o',str(P/exe)],check=True)
 C=json.load(open(src/'controls.json'));expected=json.load(open(src/'bank.json'));banks=fullbanks([b['rows']for b in expected],'literal_banks');degree=[0]*5
 for b,e in zip(banks,expected):
  need(all(b[k]==e[k]for k in ['rows','masks','coefficients','values']),'every literal named coordinate')
  for a in b['masks']:
   for c in b['masks']:degree[a.bit_count()+c.bit_count()]+=1
 # All four crossed targets and exact constant kernels.
 keys=[(p,q,r)for p,q in combinations(range(4),2)for r in range(4)if r not in(p,q)];A=tt([0,1,2,3]);inner=[]
 for z in json.load(open(src/'inner_kernels.json')):
  rows=tt(z['order']);out=[]
  for item in z['coefficients']:
   p,q,r=item['p'],item['q'],item['r'];L,R=kernel(rows,p,q,r);need((L,R)==(item['L'],item['R']),'complete inner kernels');out.append({'p':p,'q':q,'r':r,'L':L,'R':R})
  inner.append({'order':z['order'],'coefficients':out})
 cross=[v for v in permutations(range(4))if{v[0],v[-1]}=={1,2}];need(cross==[tuple(v)for v in C['crossed_targets']],'exact four residual orders');cert=[]
 for target,combo,gamma in zip(C['crossed_targets'],C['combinations'],C['constants']):
  B=tt(list(reversed(target)));Psum=[0]*3;Qsum=[0]*3;terms=[]
  for p,q,r,w in combo:
   L,R=kernel(A,p,q,r);l,rr=kernel(B,p,q,r);pp=[x-y for x,y in zip(L,l)];qq=[x-y for x,y in zip(R,rr)];Psum=[x+w*y for x,y in zip(Psum,pp)];Qsum=[x+w*y for x,y in zip(Qsum,qq)];terms.append({'p':p,'q':q,'r':r,'weight':w,'P':pp,'Q':qq})
  need(Psum==Qsum==[gamma,0,0],'four constant symmetric certificates');cert.append({'target':target,'local_B_order':list(reversed(target)),'gamma':gamma,'terms':terms})
 (P/'kernel_certificates.json').write_text(json.dumps({'complete_inner':inner,'certificates':cert},sort_keys=True,indent=2)+'\n')
 # Complete literal capacity matrices and nontrivial ambient excursion channel.
 capacities=[];individual=[];outcomes=[]
 for caseid,case in enumerate(C['cases']):
  offsets=[0,4]
  for size in case['others'][:-1]:offsets.append(offsets[-1]+size)
  sizes=[4]+case['others'];N=sum(sizes)
  def rebuild(order,flip):
   out=[0]*N
   blocks=[tt(order)]+[tt(list(range(z)))for z in case['others']]
   for bi,block in enumerate(blocks):
    for u,row in enumerate(block):
     out[offsets[bi]+u]=row<<offsets[bi]
     for bj,size in enumerate(sizes):
      if bi!=bj and bool(case['Q'][bi]>>bj&1)!=flip:out[offsets[bi]+u]|=((1<<size)-1)<<offsets[bj]
   return[outrow|(1<<N)for outrow in out]+[1<<(N+1),(1<<N)-1]if case['ambient']else out
  need(banks[case['source']]['rows']==rebuild([0,1,2,3],False),'source physical quotient provenance')
  for idx,order in zip(case['targets'],C['crossed_targets']):need(banks[idx]['rows']==rebuild(order,True),'target physical quotient provenance')
  source=banks[case['source']];rows=source['rows'];R=[r>>4 for r in rows[4:]];into=sum((r&1)<<j for j,r in enumerate(rows[4:]));V=case['virtual'];group=[i-3 for i in case['group']];Vall=V
  if case['ambient']:
   n=len(V);Vall=[r|(1<<n)for r in V]+[1<<(n+1),(1<<n)-1]
  cap=capacity(Vall);sym=sum(cap[0][i]+cap[i][0]for i in group);s=sympoly(V,0,set(group));need(s[0]>0 and all(x>=0 for x in s),'positive symmetrized polynomial')
  if case['ambient']:
   alpha,beta,tau,*_=outer([2,0],2);theta=[];power=[1,0,0,0]
   for _ in s:theta.append(mul(mul(alpha,beta),power)[-1]);power=mul(power,tau)
  else:theta=[1]+[0]*(len(s)-1)
  need(sym==sum(a*b for a,b in zip(s,theta)),'positive excursion exact transfer');capacities.append({'case':caseid,'rows':Vall,'matrix':cap,'group':group,'symmetric_capacity':sym,'core_sym_polynomial':s,'theta':theta})
  for h,(target,combo,gamma)in enumerate(zip(C['crossed_targets'],C['combinations'],C['constants'])):
   B=tt(list(reversed(target)));local=banks[case['local_replacements'][h]];actualtarget=banks[case['targets'][h]]
   # Source, target and local replacement are checked on all original physical arcs.
   need(local['rows']==physical(B,R,into),'literal local replacement rows')
   for p,q,r in keys:
    for i in range(len(R)):
     pred,_,_=replacement_formula(A,B,R,into,p,q,r,i);actual=coef(source,p,q,r,i+4)-coef(local,p,q,r,i+4);need(pred==actual,'general replacement all individual marks');individual.append([caseid,h,p,q,r,i+4,actual])
   loc=sum(w*sum(coef(source,p,q,r,i)-coef(local,p,q,r,i)for i in case['group'])for p,q,r,w in combo);act=sum(w*sum(coef(source,p,q,r,i)-coef(actualtarget,p,q,r,i)for i in case['group'])for p,q,r,w in combo);need(loc==act==gamma*sym,'uniform group quotient and ambient restoration');outcomes.append({'case':caseid,'target':h,'gamma':gamma,'symmetric_capacity':sym,'actual_grouped_contrast':act})
 need(individual==json.load(open(src/'individual_quartics.json')),'complete producer individual contrast bank')
 (P/'capacity_transfer.json').write_text(json.dumps(capacities,sort_keys=True,indent=2)+'\n');(P/'individual_quartics.json').write_text(json.dumps(individual,separators=(',',':'))+'\n');(P/'outcomes.json').write_text(json.dumps(outcomes,sort_keys=True,indent=2)+'\n')
 # All actual producer grouped queries, independently build ports on the physical rows.
 physicalq=json.load(open(src/'physical_queries.json'));requests=[]
 for q in physicalq['queries']:requests.append((banks[q['graph']]['rows'],q['a'],q['c']))
 text=str(len(requests))+'\n'+''.join(str(len(r))+'\n'+' '.join(map(str,r))+'\n'+str(a)+' '+str(c)+'\n'for r,a,c in requests);(P/'physical_queries.in').write_text(text);out=subprocess.check_output([str(P/'queries')],input=text,text=True);(P/'physical_query_values.txt').write_text(out);need(list(map(int,out.split()))==list(map(int,(src/'physical_values.txt').read_text().split())),'all physical query counts')
 # Nontransitive arbitrary replacements, without internal profile equality.
 examples=[];extra_graphs=[]
 for aa,bb in [([10,12,9,0],[2,4,9,3]),([6,12,24,17,3],[30,28,24,16,0])]:
  for erows,ein in [([2,4,1],2),([2,4,1],3),([6,4,0],2)]:
   examples.append((aa,bb,erows,ein));extra_graphs.extend([physical(aa,erows,ein),physical(bb,erows,ein)])
 extra=fullbanks(extra_graphs,'arbitrary_replacement_banks');extra_controls=[];different_channels=False;derivative_channel=False
 for j,(aa,bb,erows,ein)in enumerate(examples):
  m=len(aa);items=[]
  for p,q in combinations(range(m),2):
   for r in range(m):
    if r in[p,q]:continue
    for i in range(len(erows)):
     pred,pp,qq=replacement_formula(aa,bb,erows,ein,p,q,r,i);actual=coef(extra[2*j],p,q,r,m+i)-coef(extra[2*j+1],p,q,r,m+i);need(pred==actual,'nontransitive replacement identity');different_channels|=qq!=[-v for v in pp];derivative_channel|=any(derivative(pp))or any(derivative(qq));items.append([p,q,r,m+i,actual])
  extra_controls.append({'A':aa,'B':bb,'exterior':erows,'into':ein,'checks':items})
 need(different_channels and derivative_channel,'nonconverse and derivative channels genuinely exercised');(P/'arbitrary_replacement_controls.json').write_text(json.dumps(extra_controls,separators=(',',':'))+'\n')
 # Clean source verifier replay, without touching the producer.
 with tempfile.TemporaryDirectory(prefix='tt4-cross-audit-')as td:
  D=Path(td)
  for f in manifest:shutil.copyfile(src/f,D/f)
  for opt,name in [([], 'normal'),(['-O'],'optimized')]:
   out=subprocess.check_output([sys.executable,*opt,str(D/'verify.py')]);need(out==(src/'receipt.json').read_bytes(),'clean producer receipt replay');(P/('producer.'+name+'.json')).write_bytes(out)
 need(pins=={f:hashlib.sha256((src/f).read_bytes()).hexdigest()for f in manifest},'producer bytes unchanged')
 result={'status':'PASS','producer_named_graphs':len(banks),'producer_named_coefficients':sum(degree),'degree_counts':degree,'general_individual_quartics':len(individual),'actual_grouped_queries':len(requests),'arbitrary_replacement_graphs':len(extra),'arbitrary_replacement_coefficients':sum(len(b['coefficients'])for b in extra),'arbitrary_replacement_quartics':sum(len(x['checks'])for x in extra_controls),'producer_pins':pins,'sources':{f:hashlib.sha256((P/f).read_bytes()).hexdigest()for f in ['audit.py','algebra.py','literal_polynomial.cpp','queries.cpp']},'artifacts':{f:hashlib.sha256((P/f).read_bytes()).hexdigest()for f in ['literal_banks.json','kernel_certificates.json','capacity_transfer.json','individual_quartics.json','outcomes.json','physical_queries.in','physical_query_values.txt','arbitrary_replacement_banks.json','arbitrary_replacement_controls.json','producer.normal.json','producer.optimized.json']}};(P/'receipt.json').write_text(json.dumps(result,sort_keys=True,indent=2)+'\n');print(json.dumps(result,sort_keys=True))
if __name__=='__main__':main()
