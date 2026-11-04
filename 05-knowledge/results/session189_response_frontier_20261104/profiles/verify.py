"""Regenerate, independently cross-check, and byte-compare the complete packet."""
from pathlib import Path
import hashlib,json,shutil,subprocess,sys,tempfile
ROOT=Path(__file__).resolve().parent
SOURCES=['core.py','anchor_explore.py','prepare_models.py','prepare_extra.py','analyze_models.py','analyze_residual.py','analyze_ambient.py','check_run.py','count_named.cpp','cycle_named.cpp','run_words.cpp','anchor_check.cpp','census.cpp']
ARTIFACTS=['named_banks.json','extra_named_banks.json','affine_banks.json','anchor_exploration.json','cycle_banks.json','ambient_analysis.json','type_witnesses.json','residual_mod4.json','census_values.txt','prime_owners.u16le','run_values.txt','anchor_independent.txt','bare_full_H_values.txt','bare_full_cycle_values.txt']
def need(ok,msg):
 if not ok:raise RuntimeError(msg)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 with tempfile.TemporaryDirectory(prefix='prime7-profile-verification-') as td:
  D=Path(td)
  for f in SOURCES+['quotients.json']:shutil.copyfile(ROOT/f,D/f)
  def py(f):
   cmd=[sys.executable]+(['-O'] if sys.flags.optimize else [])+[str(D/f)]
   return subprocess.check_output(cmd,cwd=D,stderr=subprocess.PIPE).decode()
  def cc(name):subprocess.run(['g++','-std=c++17','-O3','-Wall','-Wextra','-Werror',str(D/(name+'.cpp')),'-o',str(D/name)],check=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
  def run(name,inp,out,*args):
   with open(D/inp,'rb') as a,open(D/out,'wb') as b:subprocess.run([str(D/name),*args],stdin=a,stdout=b,stderr=subprocess.PIPE,cwd=D,check=True)
  Q=json.loads((D/'quotients.json').read_text())
  need(len(Q)==197,'197 prime representatives')
  for q in Q:
   T=q['rows'];need(len(T)==7,'quotient order')
   need(all(isinstance(r,int) and 0<=r<128 and not(r>>i&1) for i,r in enumerate(T)),'valid outgoing rows')
   need(all(((T[i]>>j)&1)+((T[j]>>i)&1)==1 for i in range(7) for j in range(i+1,7)),'tournament edge complement')
   need(q['p']==[(a+b)%2 for a,b in zip(q['S'],q['E'])],'endpoint parity definition')
  for name in ['count_named','cycle_named','run_words','anchor_check','census']:cc(name)
  py('anchor_explore.py');py('prepare_models.py');run('count_named','models.in','named_values.txt');py('analyze_models.py')
  py('prepare_extra.py');run('count_named','extra_models.in','extra_named_values.txt');py('analyze_residual.py')
  M=[m for m in json.loads((D/'models.json').read_text()) if m['sizes']==[1]*7]
  (D/'cycle_models.json').write_text(json.dumps(M))
  (D/'cycle_models.in').write_text(str(len(M))+'\n'+''.join('7 '+' '.join(map(str,m['rows']))+'\n' for m in M))
  run('cycle_named','cycle_models.in','cycle_values.txt');py('analyze_ambient.py')
  S=json.loads((D/'affine_banks.json').read_text())
  (D/'run.in').write_text(str(len(S))+'\n'+''.join(' '.join(map(str,s['Q']))+'\n' for s in S));run('run_words','run.in','run_values.txt')
  (D/'census.in').write_text(str(len(Q))+'\n'+''.join(' '.join(map(str,q['rows']))+'\n' for q in Q));run('census','census.in','census_values.txt',str(D/'prime_owners.u16le'))
  zero=[q for q in Q if not any(q['p'])]
  (D/'anchor.in').write_text(str(len(zero))+'\n'+''.join(' '.join(map(str,q['rows']))+'\n' for q in zero));run('anchor_check','anchor.in','anchor_independent.txt')
  run('count_named','cycle_models.in','bare_full_H_values.txt','all');run('cycle_named','cycle_models.in','bare_full_cycle_values.txt','all')
  check_log=py('check_run.py')
  # Import only the just-copied known mathematical primitives, not environmental modules.
  sys.path.insert(0,str(D));import core
  for q in Q:
   T=q['rows'];U=core.converse(T);sig=[(s-e)%4 for s,e in zip(q['S'],q['E'])]
   for j in range(7):
    a=core.repeated_endpoint(T,j);b=core.repeated_endpoint(U,j);sig=[(v+2*(y-x))%4 for v,x,y in zip(sig,a,b)]
   need(sig==q['sigma'],'complete sigma, including endpoint-nonzero quotients')
  A=json.loads((D/'anchor_exploration.json').read_text());V=[list(map(int,l.split())) for l in (D/'anchor_independent.txt').read_text().splitlines()]
  need(len(A)==len(V)==38,'38 complete endpoint-zero records')
  for a,v in zip(A,V):need(v==sum(a['start_run_skews'],[])+a['one_port_base']+sum(a['one_port_double'],[]),'3990 independent exact anchor coordinates')
  need(sum(len(a['patterns']) for a in A)==4864,'all singleton-position patterns')
  need(sum(len(a['blind']) for a in A)==10,'ten anchor-blind patterns')
  R=json.loads((D/'residual_mod4.json').read_text());need(len(R)==6 and all(r['query'] is not None for r in R),'all six singleton residuals have nonzero complete mod4 bank')
  expected_safe=[[],[127],[0,127],[]];expected_local=[[],[127],[],[]]
  B=json.loads((D/'ambient_analysis.json').read_text())
  need([r['safe_rank'] for r in B]==[14,14,13,14],'complete safe-query ranks')
  need([r['safe_blind_masks'] for r in B]==expected_safe,'complete safe-query zero domain')
  need([r['local_blind_masks'] for r in B]==expected_local,'complete local zero domain')
  W=json.loads((D/'type_witnesses.json').read_text());need(sum(len(w['records']) for w in W)==65536,'complete four quotient type domains')
  need(sum(r[2]<0 for w in W for r in w['records'])==1,'sole omitted even carousel E7 branch')
  # Complete physical cube checks on all seven-name quotient controls.
  Hbanks=[m for m in json.loads((D/'named_banks.json').read_text()) if m['sizes']==[1]*7]
  Cbanks=json.loads((D/'cycle_banks.json').read_text());cube_checks=0
  for filename,banks in [('bare_full_H_values.txt',Hbanks),('bare_full_cycle_values.txt',Cbanks)]:
   it=iter(map(int,(D/filename).read_text().split()))
   for bank in banks:
    need(next(it)==7 and next(it)==128,'full physical cube dimensions')
    wanted=[next(it) for _ in range(16384)];co=[0]*16384
    for i,a in enumerate(bank['named_masks']):
     for j,c in enumerate(bank['named_masks']):co[(a<<7)|c]=bank['named_coefficients'][i][j]
    for bit in range(14):
     for mask in range(16384):
      if mask>>bit&1:co[mask]+=co[mask^(1<<bit)]
    need(co==wanted,'every full physical query reconstructed from complete named coefficients');cube_checks+=len(co)
   need(next(it,None) is None,'consume full cube counts')
  # Byte-level comparisons are full comparisons, not hash-only checks.
  artifact_receipts={}
  for name in ARTIFACTS:
   need((D/name).read_bytes()==(ROOT/name).read_bytes(),'retained artifact differs: '+name)
   artifact_receipts[name]={'bytes':(D/name).stat().st_size,'sha256':sha(D/name)}
  sources={f:sha(ROOT/f) for f in ['PROOF.md','verify.py','quotients.json']+SOURCES}
  return {'status':'PASS','prime_representatives':197,'labeled_tournaments_checked':2097152,'labeled_prime_tournaments':973680,'endpoint_zero_representatives':38,'exact_endpoint_coordinates':2758,'anchor_pattern_records':4864,'independent_anchor_coordinates':3990,'singleton_residuals_resolved':6,'independent_exact_run_coordinates':242208,'full_physical_cube_checks':cube_checks,'nontrivial_type_records':65536,'safe_query_ranks':[14,14,13,14],'sole_remaining_branch':'carousel (14,28,56,112,97,67,7), seven even children of order at least four','check_log':check_log.strip().splitlines(),'source_pins':sources,'artifacts':artifact_receipts}
if __name__=='__main__':print(json.dumps(main(),sort_keys=True,indent=2))
