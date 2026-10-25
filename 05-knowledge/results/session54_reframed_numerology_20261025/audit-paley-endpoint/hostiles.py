import pathlib,tempfile,shutil,json,subprocess,sys,hashlib
ROOT=pathlib.Path(__file__).resolve().parent
P=ROOT/'inputs/session54-paley-endpoint-observer'
def need(c,m):
 if not c:raise RuntimeError(m)
def validate(p):
 entries=json.loads((p/'MANIFEST.json').read_text())
 need({str(x.relative_to(p)) for x in p.rglob('*') if x.is_file()}=={x['path'] for x in entries}|{'MANIFEST.json'},'inventory')
 for x in entries:need(hashlib.sha256((p/x['path']).read_bytes()).hexdigest()==x['sha256'],'pin')
checks=[]
for kind in ['count mutation','source mutation','missing source','extra file','nested extra file','false integral product','false q19 congruence']:
 with tempfile.TemporaryDirectory() as td:
  p=pathlib.Path(td)/'packet';shutil.copytree(P,p)
  if kind=='count mutation':
   f=p/'counts.txt';f.write_text(f.read_text().replace('3040229428','3040229429'))
  if kind=='source mutation':
   f=next((p/'sources').iterdir());f.write_text(f.read_text()+'\nMUTATION\n')
  if kind=='missing source':next((p/'sources').iterdir()).unlink()
  if kind=='extra file':(p/'extra.txt').write_text('unexpected')
  if kind=='nested extra file':(p/'extra').mkdir();(p/'extra/a.txt').write_text('unexpected')
  try:
   validate(p)
   if kind=='false integral product':need((1-0-0*1)%3==0,'(1,0) fails q3 fiber modulus 3')
   if kind=='false q19 congruence':need(3040229428%19==1,'a mod19 is3 not1')
  except (RuntimeError,FileNotFoundError):checks.append(dict(hostile=kind,status='REJECTED'))
  else:raise RuntimeError('hostile accepted '+kind)
print(json.dumps(dict(status='PASS',checks=checks),indent=2))
