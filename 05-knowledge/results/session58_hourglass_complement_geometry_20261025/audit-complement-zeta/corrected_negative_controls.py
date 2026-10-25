from pathlib import Path
import subprocess,json,sys
root=Path(__file__).resolve().parent
source=(root/'replays/absolute/corrected.py').read_text()
old='if not (N - f) % 2 == 0:'
if source.count(old)!=1:raise RuntimeError('unique negative mutation site missing')
source=source.replace(old,'if not (N - f) % 2 == 1:')
path=root/'replays/absolute/corrected_negative.py';path.write_text(source)
rows=[]
for optimized in (False,True):
 r=subprocess.run([sys.executable]+(['-O'] if optimized else [])+[str(path)],capture_output=True,text=True)
 if r.returncode==0 or 'Verification failed:' not in r.stderr:raise RuntimeError('corrected negative did not fail')
 rows.append({'optimized':optimized,'exit':r.returncode,'detected':True,'mutation':'require odd involution nonfixed count'})
path.unlink()
(root/'corrected_negative_results.json').write_text(json.dumps(rows,indent=2)+'\n')
print(json.dumps(rows,indent=2))
