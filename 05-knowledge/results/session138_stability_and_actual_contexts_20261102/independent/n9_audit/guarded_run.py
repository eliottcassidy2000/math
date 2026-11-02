"""Reproduce normal/UBSan independent graph audits under recorded caps."""
import argparse,json,pathlib,resource,subprocess,time
p=argparse.ArgumentParser();p.add_argument('--compile',action='store_true');p.add_argument('--sanitize',action='store_true');a=p.parse_args()
root=pathlib.Path(__file__).resolve().parent
exe=root/('graph9_independent_ubsan' if a.sanitize else 'graph9_independent_optimized')
if a.compile:
 command=['gcc','-std=c11','-O1','-Wall','-Wextra','-Werror','-Wno-misleading-indentation']
 if a.sanitize:command+=['-fsanitize=undefined','-fno-sanitize-recover=undefined']
 command+=[str(root/'graph9_independent.c'),'-o',str(exe)]
 def cap():
  data=(64 if a.sanitize else 24)*1024*1024
  rss=(128 if a.sanitize else 40)*1024*1024
  resource.setrlimit(resource.RLIMIT_DATA,(data,data));resource.setrlimit(resource.RLIMIT_RSS,(rss,rss))
else:
 command=[str(exe)]
 def cap():
  address=(128 if a.sanitize else 40)*1024*1024
  resource.setrlimit(resource.RLIMIT_AS,(address,address))
  resource.setrlimit(resource.RLIMIT_DATA,(24*1024*1024,24*1024*1024))
  resource.setrlimit(resource.RLIMIT_RSS,(40*1024*1024,40*1024*1024))
start=time.monotonic();result=subprocess.run(command,preexec_fn=cap,capture_output=True,text=True)
print(result.stdout,end='')
receipt={'command':command,'returncode':result.returncode,'seconds':time.monotonic()-start,'max_rss_kib':resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss,'stderr':result.stderr}
name='guarded_'+('compile' if a.compile else 'replay')+('_ubsan' if a.sanitize else '')+'.json'
(root/name).write_text(json.dumps(receipt,indent=2)+'\n')
if result.returncode:raise SystemExit(result.returncode)
if not a.compile and receipt['max_rss_kib']>=40*1024:raise RuntimeError('runtime RAM budget exceeded')
