import pathlib,tempfile,shutil,json,subprocess,sys
D=pathlib.Path(__file__).resolve().parent
out=[]
for name,target,change in [
 ('wrong_exponent','converse.json',lambda d:d['map'][0].__setitem__('exponent',(d['map'][0]['exponent']+1)%63)),
 ('missing_blocker','result.json',lambda d:d['bijection'].pop()),
 ('wrong_subset','result.json',lambda d:d['bijection'][0].__setitem__('even_subset',[0,1])),
 ('wrong_fourier','result.json',lambda d:d['bijection'][0].__setitem__('fourier',[0,0]))]:
 with tempfile.TemporaryDirectory() as t:
  p=pathlib.Path(t);shutil.copy(D/'audit.py',p/'audit.py');shutil.copytree(D/'frozen',p/'frozen')
  f=p/'frozen'/target;data=json.loads(f.read_text());change(data);f.write_text(json.dumps(data))
  for opt in ([],['-O']):
   r=subprocess.run([sys.executable,*opt,str(p/'audit.py')],capture_output=True,text=True)
   if r.returncode==0:raise RuntimeError('accepted mutation '+name)
   out.append({'mutation':name,'optimized':bool(opt),'returncode':r.returncode,'reason':r.stderr.strip().splitlines()[-1]})
print(json.dumps({'status':'PASS','rejections':out},indent=2))
