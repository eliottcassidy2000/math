from pathlib import Path
import json,subprocess,sys,time,resource
from dyadic_completion import integer_interval,recover_from_interval
sys.set_int_max_str_digits(0)
P=Path(__file__).parent;R=P.parents[1]
values=list(map(int,json.loads((R/'session270_chanpak_enumeration_20261111/enumeration/prefix.json').read_text())['values']))
binary='/tmp/s277_modular'
subprocess.run(['g++','-O3','-std=c++17',str(P/'modular_counter.cpp'),'-lgmpxx','-lgmp','-o',binary],check=True)
def limits():resource.setrlimit(resource.RLIMIT_AS,(800*1024**2,800*1024**2))
controls=0
for n in range(31):
 for b in (0,1,5,20,99,300):
  q=subprocess.run([binary,str(n),str(b)],capture_output=True,text=True,check=True,preexec_fn=limits,timeout=15)
  assert int(q.stdout)==values[n]%(1<<b);controls+=1
rows=[]
for n in (100,150,200,250):
 lo,hi=integer_interval(n,3);b=(hi-lo).bit_length();start=time.perf_counter()
 try:
  q=subprocess.run([binary,str(n),str(b)],capture_output=True,text=True,preexec_fn=limits,timeout=90)
  elapsed=time.perf_counter()-start
  row=dict(n=n,bits=b,elapsed=elapsed,returncode=q.returncode)
  if q.returncode==0:
   assert recover_from_interval(lo,hi,int(q.stdout),b)==values[n]
   row['status']='PASS';row.update(json.loads(q.stderr))
  else:row.update(status='RESOURCE_BLOCKED',stderr=q.stderr)
 except subprocess.TimeoutExpired:row=dict(n=n,bits=b,status='TIMEOUT_90_SECONDS')
 rows.append(row);print(row,flush=True)
 (P/'cpp_results.json').write_text(json.dumps(dict(controls=controls,rows=rows,memory_limit_bytes=800*1024**2,per_large_case_timeout_seconds=90),indent=2)+'\n')
