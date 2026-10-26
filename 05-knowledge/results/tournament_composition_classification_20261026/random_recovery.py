"""Replay the frozen random repository sample, independent of later additions."""
import hashlib
import json
from pathlib import Path
import random
import subprocess
import sys

sys.stdout.reconfigure(encoding='utf-8', newline='\n')
REPO=Path(__file__).resolve().parents[3]
COMMIT='fc3768ac5fcd'
ROOTS=['01-canon/theorems','05-knowledge/results','07-reflections']
raw=subprocess.check_output(['git','ls-tree','-r','--name-only',COMMIT,'--']+ROOTS,cwd=REPO,text=True)
paths=[s for s in raw.splitlines() if s.endswith('.md') and 'normseven_arithmetic_f1_20261026' not in s]
# Match rg --files' Windows spelling used for the actual first sample.
windows=[]
for p in paths:
    prefix=next(r for r in ROOTS if p.startswith(r+'/'))
    windows.append(prefix+'\\'+p[len(prefix)+1:].replace('/','\\'))
population=sorted(windows)
chosen=random.Random(812720261026).sample(population,24)
records=[]
for path in chosen:
    posix=path.replace('\\','/')
    blob=subprocess.check_output(['git','show',COMMIT+':'+posix],cwd=REPO)
    lines=blob.decode('utf-8',errors='replace').splitlines()
    title=next((line for line in lines if line.startswith(('title:','# '))), '(no title)')
    records.append({'path':posix,'title':title,'sha256':hashlib.sha256(blob).hexdigest()})
print(json.dumps({'commit':COMMIT,'seed':812720261026,'population_size':len(population),
                  'population_sha256':hashlib.sha256('\n'.join(population).encode()).hexdigest(),
                  'sample':records},indent=2,ensure_ascii=False))
