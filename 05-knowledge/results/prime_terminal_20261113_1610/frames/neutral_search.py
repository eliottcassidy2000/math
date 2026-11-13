"""All at-most-four arc flips of one specified order-ten target."""
from pathlib import Path
import json
import subprocess
HERE=Path(__file__).resolve().parent
exe=HERE/'neutral10.exe'
subprocess.run(['g++','-O3','-std=c++17',str(HERE/'neutral10.cpp'),'-o',str(exe)],check=True)
out=json.loads(subprocess.check_output([str(exe)],text=True))
out['scope']='All distance<=4 arc reversals from the named H601 Bezout witness; no all-order negative inference'
(HERE/'neutral_receipt.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8')
print({k:v for k,v in out.items() if k!='hits'})
print('hits',len(out['hits']),out['hits'][:2])
