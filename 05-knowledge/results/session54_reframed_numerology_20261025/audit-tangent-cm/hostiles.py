import json, pathlib, subprocess, tempfile
base=pathlib.Path(__file__).resolve().parent; src=(base/'inputs'/'producer-results.json').read_text(); controls=[]
for name,script,mutate in [
 ('alter-j8000','independent.py',lambda d:d['candidate_conductors'][1].update(polynomial='J - 8001')),
 ('omit-n168','independent.py',lambda d:d['candidate_conductors'].pop()),
 ('remove-nonfundamental-minus12','audit_small_discriminants.py',lambda d:d['reduced_forms_abs_discriminant_below32'].update({'-12':[]})),
 ('wrong-trace-bound','audit_small_discriminants.py',lambda d:d['discriminant_bounds'][0].update(trace=10**100))]:
 with tempfile.TemporaryDirectory(dir=base) as tmp:
  temp=pathlib.Path(tmp); d=json.loads(src);mutate(d); p=temp/'data.json';p.write_text(json.dumps(d))
  code=(base/script).read_text().replace("(Path(__file__).resolve().parent/'inputs'/'producer-results.json')", "Path("+repr(str(p))+")")
  test=temp/'test.py';test.write_text(code)
  for opt in [[],['-O']]:
   run=subprocess.run(['python',*opt,str(test)],capture_output=True,text=True)
   if run.returncode==0:raise RuntimeError('hostile accepted '+name)
   controls.append({'name':name,'optimized':bool(opt),'exit':run.returncode,'last_line':run.stderr.splitlines()[-1]})
print(json.dumps(controls,indent=2))
