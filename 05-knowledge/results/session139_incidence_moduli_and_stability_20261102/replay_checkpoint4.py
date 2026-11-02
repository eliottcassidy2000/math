#!/usr/bin/env python3
"""Fresh exact replay of final Session139 accepted inputs and active controls."""
from pathlib import Path
import argparse,hashlib,json,os,shutil,subprocess,sys,tempfile

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path);args=ap.parse_args();src=Path(__file__).resolve().parent
 def check(ok,label):
  if not ok:raise RuntimeError(label)
 kernels={
 'tensor_matching_moduli/verify_moduli.py':'tensor_matching_moduli/results.json',
 'tensor_matching_audit/verify_exact.py':'tensor_matching_audit/exact_results.json',
 'tensor_matching_audit/verify_mutations.py':'tensor_matching_audit/mutation_results.json',
 'tensor_matching_audit/verify_arithmetic_slices.py':'tensor_matching_audit/arithmetic_slices.json',
 'stability/verify_eisenstein_bernstein.py':'stability/eisenstein_bernstein.json',
 'stability/verify_n140.py':'stability/verification_n140.json',
 'independent_threshold140/audit_n140.py':'independent_threshold140/normal.json',
 'stability/verify_local_factor_hostiles.py':'stability/local_factor_hostiles.json',
 'independent/audit_local_factors.py':'independent/audit_local_factors.out',
 'all_moduli/verify_padding_orbits.py':'all_moduli/padding.normal.json'}
 provenance=['tensor_matching_audit/verify_curated_manifest.py','stability/verify_curated_checkpoint4.py']
 manifest=json.loads((src/'stability/CURATED_CHECKPOINT4.json').read_text())
 files={x['path'] for sec in ['producer_files','independent_evidence'] for x in manifest[sec]}
 files.update(manifest['generated_manifest_exclusions']);files.update(kernels);files.add('all_moduli/PADDING_SOURCE_PINS.json')
 mutations=[
 ('tensor_matching_moduli/verify_moduli.py','require(gap>=0,(key,str(gap)))','require(gap>0,(key,str(gap)))','RuntimeError'),
 ('tensor_matching_audit/verify_exact.py','require(gap>=0,(key,gap))','require(gap>0,(key,gap))','RuntimeError'),
 ('stability/verify_eisenstein_bernstein.py','quadratic([-1,1],[2])==1<2','quadratic([-1,1],[2])==2','negative_control_non_even_constant'),
 ('independent/audit_local_factors.py','F(-7,8192)','F(7,8192)','n5 conditional factor hostile'),
 ('all_moduli/verify_padding_orbits.py','fixed(11,2)==21','fixed(11,2)==22','odd fixed/zero control')]
 flags=[('stability/verify_n140.py',['--threshold','139'],'bootstrap_step_0'),('independent_threshold140/audit_n140.py',['--threshold','139'],'bootstrap_step_0'),('independent_threshold140/audit_n140.py',['--mutate-perfect-indicator'],'perfect_penalty_case_n140_r1')]
 runs=[];outputs={}
 with tempfile.TemporaryDirectory(prefix='session139-checkpoint4-relocated-') as name:
  dst=Path(name)
  for f in files:(dst/f).parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(src/f,dst/f)
  mutated=[]
  for i,(f,old,new,needle) in enumerate(mutations):
   text=(src/f).read_text();check(text.count(old)==1,'unique mutation '+f);target=str(Path(f).parent/f'mutation_{i}.py');(dst/target).write_text(text.replace(old,new,1));mutated.append((f,target,needle))
  def run(label,script,flags=(),expected=0,needle=None,opt=False,kind='scientific'):
   proc=subprocess.run([sys.executable,'-B']+(['-O'] if opt else [])+[script,*flags],cwd=dst,text=True,capture_output=True,env={**os.environ,'PYTHONDONTWRITEBYTECODE':'1'},timeout=900)
   check(proc.returncode==expected,f'{label}: exit {proc.returncode}: {proc.stderr[-1500:]}')
   if needle:check(needle in proc.stdout+proc.stderr,f'{label}: missing rejection '+needle)
   runs.append({'label':label,'kind':kind,'mode':'optimized' if opt else 'normal','exit_code':proc.returncode,'stdout_sha256':hashlib.sha256(proc.stdout.encode()).hexdigest(),'expected_rejection':expected!=0});return proc.stdout
  for opt in (False,True):
   for f in kernels:outputs[f,opt]=run(f,f,opt=opt)
   for f in provenance:outputs[f,opt]=run(f,f,opt=opt,kind='source-integrity')
   for f,flag,needle in flags:run(f+' '+' '.join(flag),f,flag,1,needle,opt)
   for f,target,needle in mutated:run(f+'-negative',target,expected=1,needle=needle,opt=opt)
   proof=dst/'tensor_matching/PROOF.md';original=proof.read_bytes();altered=original.replace(b'# Tensor',b'# Xensor',1);check(len(original)==len(altered) and altered!=original,'same-size proof corruption')
   proof.write_bytes(altered)
   run('tensor-manifest-corruption',provenance[0],expected=1,needle='SHA256 mismatch',opt=opt,kind='source-integrity')
   run('producer-manifest-corruption',provenance[1],expected=1,needle='pin mismatch',opt=opt,kind='source-integrity')
   proof.write_bytes(original)
  for f,receipt in kernels.items():
   check(outputs[f,False]==outputs[f,True],f+' normal/optimized identity');check(outputs[f,False]==(src/receipt).read_text(),f+' frozen output')
  for f in provenance:check(outputs[f,False]==outputs[f,True],f+' normal/optimized identity')
  check(outputs[provenance[1],False]==(src/'stability/checkpoint4_manifest_verification.json').read_text(),'producer manifest receipt')
  report={'status':'PASS','fresh_relocation':True,'total_executions':len(runs),'positive_scientific_runs':20,'positive_source_integrity_runs':4,'expected_rejections':20,'mathematical_false_strengthenings_rejected_inside_tensor_control':8,'appendix_internal_negative_controls':4,'normal_optimized_identical':[*kernels,*provenance],'frozen_outputs_matched':list(kernels.values()),'n140_independent_named_checks':len(json.loads(outputs['independent_threshold140/audit_n140.py',False])['active_checks']),'source_pins':{f:hashlib.sha256((src/f).read_bytes()).hexdigest() for f in sorted(files)},'runs':runs}
 body=json.dumps(report,sort_keys=True,indent=2)+'\n'
 if args.output:args.output.write_text(body)
 print(body,end='')
if __name__=='__main__':main()
