#!/usr/bin/env python3
"""Fresh relocated checkpoint-two replay; no inherited research imports.
Use python3 replay_checkpoint2.py [--output report.json].
"""
from pathlib import Path
import argparse,hashlib,json,os,shutil,subprocess,sys,tempfile

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path);args=ap.parse_args();src=Path(__file__).resolve().parent
 def check(ok,label):
  if not ok:raise RuntimeError(label)
 kernels={
  'stability/verify_n526.py':'stability/verification_n526.json',
  'independent_threshold526/audit_n526.py':'independent_threshold526/normal.json',
  'all_moduli/verify.py':'all_moduli/normal.json',
  'independent/audit_moduli.py':'independent/audit_moduli.out',
  'independent/audit_function_states.py':'independent/audit_function_states.out',
  'independent/audit_allmod_padding.py':'independent/audit_allmod_padding.out'}
 mutations=[
  ('stability/verify_n526.py','lam=d/64','lam=d/63','dual_polynomial_identity'),
  ('all_moduli/verify.py','clock(16)==(4,2)','clock(16)==(4,4)','dyadic special period'),
  ('independent/audit_moduli.py','predicted(16)==(4,2)','predicted(16)==(4,4)','dyadic period is not Euler period'),
  ('independent/audit_function_states.py','expected=2*p**((p-1)//2)-1','expected=2*p**((p-1)//2)-2','full unary orbit incomplete'),
  ('independent/audit_allmod_padding.py','response==172','response==173','actual E4 response')]
 runs=[];outputs={}
 with tempfile.TemporaryDirectory(prefix='session139-checkpoint2-relocated-') as name:
  dst=Path(name)
  for f in kernels:(dst/f).parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(src/f,dst/f)
  for i,(f,old,new,needle) in enumerate(mutations):
   text=(src/f).read_text();check(text.count(old)==1,'unique mutation '+f);(dst/f'mutation_{i}.py').write_text(text.replace(old,new,1))
  def run(label,script,flags=(),expected=0,needle=None,opt=False):
   command=[sys.executable,'-B']+(['-O'] if opt else [])+[script,*flags]
   proc=subprocess.run(command,cwd=dst,text=True,capture_output=True,env={**os.environ,'PYTHONDONTWRITEBYTECODE':'1'},timeout=600)
   check(proc.returncode==expected,f'{label}: exit {proc.returncode}: {proc.stderr[-1500:]}')
   if needle:check(needle in proc.stdout+proc.stderr,f'{label}: missing rejection '+needle)
   runs.append({'label':label,'mode':'optimized' if opt else 'normal','exit_code':proc.returncode,'stdout_sha256':hashlib.sha256(proc.stdout.encode()).hexdigest(),'expected_rejection':expected!=0})
   return proc.stdout
  for opt in (False,True):
   for f in kernels:outputs[f,opt]=run(f,f,opt=opt)
   run('threshold525-negative','independent_threshold526/audit_n526.py',('--threshold','525'),1,'final_strict_margin',opt)
   for i,(f,old,new,needle) in enumerate(mutations):run(f+'-negative',f'mutation_{i}.py',expected=1,needle=needle,opt=opt)
  for f,receipt in kernels.items():
   check(outputs[f,False]==outputs[f,True],f+' normal/optimized identity')
   check(outputs[f,False]==(src/receipt).read_text(),f+' frozen receipt identity')
  report={'status':'PASS','fresh_relocation':True,'normal_optimized_identical':list(kernels),'scientific_executions':len(runs),'positive_runs':12,'expected_rejections':12,'frozen_outputs_matched':list(kernels.values()),'n526_independent_named_checks':len(json.loads(outputs['independent_threshold526/audit_n526.py',False])['active_checks']),'source_pins':{f:hashlib.sha256((src/f).read_bytes()).hexdigest() for f in kernels},'runs':runs}
 body=json.dumps(report,sort_keys=True,indent=2)+'\n'
 if args.output:args.output.write_text(body)
 print(body,end='')
if __name__=='__main__':main()
