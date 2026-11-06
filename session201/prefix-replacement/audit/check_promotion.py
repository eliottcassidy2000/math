"""Verify exact approved promotion deltas and replay final/preserved artifacts."""
from pathlib import Path
import hashlib,json,subprocess,sys,tempfile,shutil
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent
FINAL=HERE/'promoted-inputs'
OLD=HERE/'frozen-inputs'
EXPECT={
 'PROOF.md':'fac1169ad772839fd80be856c67d05d01d01f703404456f24a139b8d738c54e0',
 'DIGIT-ANSATZ.md':'cd7b3fcd2a11e9542e2312d02faade85c508c35f503452fd9c85c780d67d33c7',
 'DYADIC-SPAN.md':'9f37c511c6dcd1cda8fb292826470d9218e55fd9b02bb0e93232c872e1c4d859',
 'README.md':'06c316a94534acd652e5e28255f915b91822659bde05ed417e2d8bf969b562fd',
 'verify.py':'2868922f90f5ba2cc8d6122408d87df692bfd64d1afa5747b814ca4f070af8dc',
 'check_digit_ansatz.py':'c4bf47af37fd6839309d9d4def09b7ec7dee01e147c64bdd6e0979a864db56b8',
 'check_dyadic_span.py':'a121e03109243afe3218a1ac56fd2cb87a2b9b90c951b8a36d2179042533108d',
 'VALIDATION-normal.json':'bb10aecd4f968a2868737cd88525d0c5b374ff1119e7e0e21e3fe4197dec7da1',
 'VALIDATION-optimized.json':'bb10aecd4f968a2868737cd88525d0c5b374ff1119e7e0e21e3fe4197dec7da1',
 'DIGIT-normal.json':'caeda63b8c707e6d8b4fcbe15fac125f8cca1df5e9acb2bf200938348a46ee1c',
 'DIGIT-optimized.json':'caeda63b8c707e6d8b4fcbe15fac125f8cca1df5e9acb2bf200938348a46ee1c',
 'DYADIC-normal.json':'dd30a0decff4e51eef6c13f74a31b4c984c8475ceca47e990f1797babd172149',
 'DYADIC-optimized.json':'dd30a0decff4e51eef6c13f74a31b4c984c8475ceca47e990f1797babd172149',
}
def need(x,label):
 if not x:raise RuntimeError(label)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def exact_replacement(name,pairs):
 text=(OLD/name).read_text()
 for before,after in pairs:
  need(text.count(before)==1,('unique delta anchor',name,before))
  text=text.replace(before,after)
 need(text==(FINAL/name).read_text(),('exact approved delta',name))
def main():
 for name,h in EXPECT.items():need(sha(FINAL/name)==h,('final input hash',name))
 exact_replacement('PROOF.md',[
  ('**PROVED, pending independent audit.**','**PROVED, independently audited in the stated scope.**'),
  ('The search has55 profiles. Its completeness follows from Section2, so the',
   'The U-only inventory has55 profiles, of which34 also have D<=33. Its\ncompleteness follows from Section2, so the')])
 exact_replacement('DIGIT-ANSATZ.md',[
  ('**PROVED, pending independent audit; separate stopping-rule sidecar.**',
   '**PROVED, independently audited; separate stopping-rule sidecar.**')])
 exact_replacement('DYADIC-SPAN.md',[
  ('**PROVED, pending independent audit; separate addendum to DIGIT-ANSATZ.md.**',
   '**PROVED, independently audited; separate addendum to DIGIT-ANSATZ.md.**'),
  ('## 1. A three-dimensional affine span in every tail',
   '## 1. Three linearly independent homogenized vectors in every tail')])
 exact_replacement('check_dyadic_span.py',[
  ('Every dyadic class has actual affine span3;',
   'Every dyadic class has actual homogenized linear span3;')])
 for name in ['verify.py','check_digit_ansatz.py','VALIDATION-normal.json','VALIDATION-optimized.json','DIGIT-normal.json','DIGIT-optimized.json']:
  need((FINAL/name).read_bytes()==(OLD/name).read_bytes(),('unchanged main/digit',name))
 for p in (OLD/'dependencies').glob('*.md'):
  need(p.read_bytes()==(FINAL/'dependencies'/p.name).read_bytes(),('unchanged dependency',p.name))
 old=json.loads((OLD/'DYADIC-normal.json').read_text());new=json.loads((FINAL/'DYADIC-normal.json').read_text())
 need({k:v for k,v in old.items()if k not in ['scope','verifier_sha256']}=={k:v for k,v in new.items()if k not in ['scope','verifier_sha256']},'unchanged dyadic numbers')
 runs=[]
 # A fresh tree checks that the documented paths do not rely on the scratch
 # packet name or its neighboring files. Copy only the needed audit sources.
 with tempfile.TemporaryDirectory(prefix='prefix-repair-release-check-') as tmp:
  root=Path(tmp)/'packet';shutil.copytree(FINAL,root)
  audit=root/'audit';audit.mkdir();shutil.copytree(OLD,audit/'frozen-inputs')
  for script,receipt in [('check_independent.py','independent-normal.json'),('check_digit_independent.py','digit-independent-normal.json'),('check_dyadic_independent.py','dyadic-independent-normal.json')]:
   shutil.copy2(HERE/script,audit/script)
   for optimized in [False,True]:
    cmd=[sys.executable,'-B']+(['-O']if optimized else[])+['audit/'+script]
    result=subprocess.run(cmd,cwd=root,check=True,capture_output=True).stdout
    need(result==(HERE/receipt).read_bytes(),('portable independent replay',script,optimized))
    runs.append(dict(script='audit/'+script,optimized=optimized,status='MATCH'))
  for script,receipt in [('verify.py','VALIDATION-normal.json'),('check_digit_ansatz.py','DIGIT-normal.json'),('check_dyadic_span.py','DYADIC-normal.json')]:
   for optimized in [False,True]:
    cmd=[sys.executable,'-B']+(['-O']if optimized else[])+[script]
    result=subprocess.run(cmd,cwd=root,check=True,capture_output=True).stdout
    need(result==(FINAL/receipt).read_bytes(),('portable author replay',script,optimized))
    runs.append(dict(script=script,optimized=optimized,status='MATCH'))
  demos=[]
  for args in [['--repair-word','ZOOZ'],['--repair-word','OZJJ','--target','5']]:
   result=json.loads(subprocess.run([sys.executable,'-B','verify.py']+args,cwd=root,check=True,capture_output=True).stdout)
   need(result['status']=='REPAIRED','README demo success')
   demos.append(result)
  need(demos[0]['retained_prefix']=='Z'and demos[0]['removed_letters']==3 and len(demos[0]['extension'])==9,'README target33 scope')
 return dict(status='PASS',final_input_hashes=EXPECT,approved_deltas_only=True,portable_normal_and_optimized_replays=runs,README_demos=demos,verifier_sha256=sha(Path(__file__)),scope='Final promotion and portable packaging verification; all three previously audited mathematical scopes remain unchanged.')
if __name__=='__main__':print(json.dumps(main(),sort_keys=True,indent=2))
