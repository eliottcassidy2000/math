import pathlib,tempfile,shutil,importlib.util,json
P=pathlib.Path(__file__).resolve().parent/'producer'
s=importlib.util.spec_from_file_location('packet_replay',P/'replay.py');m=importlib.util.module_from_spec(s);s.loader.exec_module(m)
out=[]
for kind in ['missing','nested-extra','modified']:
 with tempfile.TemporaryDirectory() as td:
  q=pathlib.Path(td)/'packet';shutil.copytree(P,q)
  if kind=='missing':(q/'independent.json').unlink()
  if kind=='nested-extra':(q/'nested').mkdir();(q/'nested'/'extra').write_text('x')
  if kind=='modified':(q/'fixed-points.json').write_bytes((q/'fixed-points.json').read_bytes()+b' ')
  try:m.validate(q)
  except RuntimeError as e:out.append({'hostile':kind,'rejected':True,'reason':str(e)})
  else:raise RuntimeError('hostile accepted: '+kind)
print(json.dumps(out,indent=2))
