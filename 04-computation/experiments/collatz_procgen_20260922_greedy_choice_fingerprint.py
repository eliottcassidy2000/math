import subprocess,re
import os, tempfile
HERE=os.path.dirname(os.path.abspath(__file__)); _t=tempfile.mkdtemp(); EXE=os.path.join(_t,'excgen'); subprocess.run(['clang','-O3','-o',EXE,os.path.join(HERE,'collatz_procgen_20260922_exceptional_general.c'),'-lm'],check=True); M=20; J=6
def count(S):
    if not S: args=[EXE,'3','1',str(M),'0']
    else: args=[EXE,'3','1',str(M),'2',str(J)]+[str(s) for s in sorted(S)]
    out=subprocess.run(args,capture_output=True,text=True).stdout.strip().splitlines()[-1]
    return int(re.search(r'exceptional=(\d+)',out).group(1))
evens=list(range(0,64,2))
S=set(); base=count(S); print("S=empty:",base, flush=True)
full=count(set(evens)); print("S=all evens mod 64:",full, flush=True)
for rnd in range(8):
    best=None
    for r in evens:
        if r in S: continue
        c=count(S|{r})
        if best is None or c<best[0]: best=(c,r)
    S.add(best[1]); print(f"round {rnd+1}: add {best[1]:2d} (mod 64, = {best[1]%8} mod 8) -> exceptional {best[0]}  S={sorted(S)}", flush=True)
