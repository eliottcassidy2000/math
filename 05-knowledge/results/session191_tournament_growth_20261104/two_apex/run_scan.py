#!/usr/bin/env python3
import hashlib,json,subprocess
from pathlib import Path
from verify import dp,graph,require
HERE=Path(__file__).resolve().parent
subprocess.run(['g++','-O3','-std=c++17',str(HERE/'scan.cpp'),'-o',str(HERE/'scan.bin')],check=True)
r=subprocess.run([str(HERE/'scan.bin'),'12'],cwd=HERE,check=True,capture_output=True,text=True)
out=json.loads(r.stdout);require(out=={'chain_length_max':12,'presentations':22369621,'distinct_values':200235,'maximum_H':744405,'odd_prefix_except_7_21':230519,'first_additional_finite_gap':230521},'scan output')
raw=(HERE/'scan_witnesses.tsv').read_bytes();out['witnesses_sha256']=hashlib.sha256(raw).hexdigest()
small=[];total=0;larger=[]
for line in raw.decode().splitlines():
 h,w=line.split('\t');h=int(h);word=tuple(map(int,w));total+=1
 if h<=1001:
  require(dp(graph(word))==h,'independent witness '+str(h));small.append({'H':h,'word':w})
 # Additional deterministic witnesses across the full prefix and upper band.
 elif h%10009==1:
  require(dp(graph(word))==h,'larger independent witness '+str(h));larger.append(h)
require(total==200235,'witness support size')
require([x['H'] for x in small]==[h for h in range(1,1002,2) if h not in (7,21)],'complete small prefix')
out['independently_checked_prefix_through']=1001;out['independently_checked_witnesses']=len(small);out['larger_independent_witnesses']=larger
(HERE/'prefix_witnesses.json').write_text(json.dumps(small,indent=2)+'\n')
(HERE/'scan_receipt.json').write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
print(json.dumps(out,sort_keys=True))
