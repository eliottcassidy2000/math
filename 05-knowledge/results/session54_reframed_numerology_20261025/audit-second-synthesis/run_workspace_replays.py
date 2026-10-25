from pathlib import Path
import subprocess,concurrent.futures,json,sys,os
A=Path(__file__).resolve().parent;R=A.parent
names=['sixtythree','p11','budget-transfer','boundary-polarity','cycle-packing','cycle-potential']
def run(name):
 p=R/('session54-audit-'+name)/'replay.py';s=subprocess.run([sys.executable,'-B',str(p),'--self-test'],cwd='/tmp',capture_output=True,text=True,env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1'),timeout=1200)
 if s.returncode:raise RuntimeError(name+' '+s.stdout+s.stderr)
 (A/(name+'-replay.txt')).write_text(s.stdout)
 return {'packet':name,'status':'PASS','command':['python','-B','replay.py','--self-test']}
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as ex:rows=list(ex.map(run,names))
(A/'workspace-replays.json').write_text(json.dumps(rows,indent=2)+'\n')
print(json.dumps(rows))
