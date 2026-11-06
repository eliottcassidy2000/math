"""Run one independently compiled checker and retain exact output and resources."""
import json, platform, resource, subprocess, sys, time
if len(sys.argv) != 3:
    raise SystemExit('usage: python3 run_check.py EXECUTABLE RECEIPT_PREFIX')
exe, prefix = sys.argv[1:]
t = time.perf_counter()
with open(prefix+'.json', 'w') as out, open(prefix+'_PROGRESS.txt', 'w') as err:
    proc = subprocess.run([exe], stdout=out, stderr=err)
r = resource.getrusage(resource.RUSAGE_CHILDREN)
with open(prefix+'_RESOURCES.json', 'w') as out:
    json.dump(dict(exit_code=proc.returncode, wall_seconds=time.perf_counter()-t,
                   child_user_seconds=r.ru_utime, child_system_seconds=r.ru_stime,
                   max_rss_KiB=r.ru_maxrss, machine=platform.machine(),
                   platform=platform.platform()), out, indent=2)
    out.write('\n')
if proc.returncode:
    raise SystemExit(proc.returncode)
