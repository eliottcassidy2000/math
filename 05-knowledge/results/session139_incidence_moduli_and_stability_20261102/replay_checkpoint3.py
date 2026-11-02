#!/usr/bin/env python3
"""Fresh relocated checkpoint-three exact replay with active failure controls."""
from pathlib import Path
import argparse,hashlib,json,os,shutil,subprocess,sys,tempfile

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path);args=ap.parse_args();src=Path(__file__).resolve().parent
 def check(ok,label):
  if not ok:raise RuntimeError(label)
 kernels={
 'independent/audit_integral_cocycle.py':'independent/audit_integral_cocycle.out',
 'independent/audit_weighted_cocycle.py':'independent/audit_weighted_cocycle.out',
 'regular_spectral/verify_regular_critical_line.py':'regular_spectral/verify_regular_critical_line.out',
 'independent/audit_regular_spectrum.py':'independent/audit_regular_spectrum.out',
 'independent/audit_regular_circle.py':'independent/audit_regular_circle.out',
 'all_moduli/verify_content.py':'all_moduli/content.normal.json',
 'all_moduli/verify_weighted_simulation.py':'all_moduli/weighted.normal.json',
 'independent/audit_content_simulation.py':'independent/audit_content_simulation.out',
 'independent/audit_regular_content_bridge.py':'independent/audit_regular_content_bridge.out'}
 dependencies=['independent/audit_bridges.py','all_moduli/content.normal.json']
 flag_controls=[
 ('independent/audit_integral_cocycle.py','--mutate-five-sign','oriented five-face integral conservation'),
 ('independent/audit_weighted_cocycle.py','--mutate-second-constant','RuntimeError'),
 ('regular_spectral/verify_regular_critical_line.py','--mutate-conjugation','conjugacy'),
 ('regular_spectral/verify_regular_critical_line.py','--mutate-difference','H difference'),
 ('regular_spectral/verify_regular_critical_line.py','--mutate-functional-sign','reduced polynomial functional equation'),
 ('independent/audit_regular_spectrum.py','--mutate-h-coefficient','RuntimeError')]
 mutations=[
 ('independent/audit_regular_circle.py','*product_num*d*F(1,4)','*product_num*d*F(1,5)','exact circle H identity'),
 ('all_moduli/verify_content.py','check(g==want,label)','check(g==want+1,label)','RuntimeError'),
 ('all_moduli/verify_weighted_simulation.py','FA==[0,Q(3,2),0,5,0,1]','FA==[0,Q(5,2),0,5,0,1]','weighted A F'),
 ('independent/audit_content_simulation.py',"wants={'W5':1,","wants={'W5':2,",'content'),
 ('independent/audit_regular_content_bridge.py',"==[880,660]","==[880,880]",'noncospectral hostile')]
 runs=[];outputs={}
 with tempfile.TemporaryDirectory(prefix='session139-checkpoint3-relocated-') as name:
  dst=Path(name)
  for f in [*kernels,*dependencies]:(dst/f).parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(src/f,dst/f)
  mutated=[]
  for i,(f,old,new,needle) in enumerate(mutations):
   text=(src/f).read_text();check(text.count(old)==1,'unique mutation '+f);target=str(Path(f).parent/f'mutation_{i}.py');(dst/target).write_text(text.replace(old,new,1));mutated.append((f,target,needle))
  def run(label,script,flags=(),expected=0,needle=None,opt=False):
   proc=subprocess.run([sys.executable,'-B']+(['-O'] if opt else [])+[script,*flags],cwd=dst,text=True,capture_output=True,env={**os.environ,'PYTHONDONTWRITEBYTECODE':'1'},timeout=900)
   check(proc.returncode==expected,f'{label}: exit {proc.returncode}: {proc.stderr[-1500:]}')
   if needle:check(needle in proc.stdout+proc.stderr,f'{label}: missing rejection '+needle)
   runs.append({'label':label,'mode':'optimized' if opt else 'normal','exit_code':proc.returncode,'stdout_sha256':hashlib.sha256(proc.stdout.encode()).hexdigest(),'expected_rejection':expected!=0});return proc.stdout
  for opt in (False,True):
   for f in kernels:outputs[f,opt]=run(f,f,opt=opt)
   for f,flag,needle in flag_controls:run(f+flag,f,(flag,),1,needle,opt)
   for f,target,needle in mutated:run(f+'-negative',target,expected=1,needle=needle,opt=opt)
  for f,receipt in kernels.items():
   check(outputs[f,False]==outputs[f,True],f+' normal/optimized identity')
   check(outputs[f,False]==(src/receipt).read_text(),f+' frozen receipt identity')
  report={'status':'PASS','fresh_relocation':True,'scientific_executions':len(runs),'positive_runs':18,'expected_rejections':22,'normal_optimized_identical':list(kernels),'frozen_outputs_matched':list(kernels.values()),'imported_independent_helpers':['independent/audit_bridges.py','independent/audit_integral_cocycle.py'],'pinned_input_receipt':'all_moduli/content.normal.json','source_pins':{f:hashlib.sha256((src/f).read_bytes()).hexdigest() for f in [*kernels,*dependencies]},'runs':runs}
 body=json.dumps(report,sort_keys=True,indent=2)+'\n'
 if args.output:args.output.write_text(body)
 print(body,end='')
if __name__=='__main__':main()
