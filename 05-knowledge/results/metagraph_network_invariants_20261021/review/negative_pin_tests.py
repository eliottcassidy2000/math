#!/usr/bin/env python3
"""Pin enforcement fails closed, including under python -O. No original writes."""
import hashlib,json,shutil,subprocess,sys,tempfile
from pathlib import Path
H=Path(__file__).resolve().parent;repo=H.parent/'math';producer=H.parent/'session20-metagraph-networks'
pins=json.loads((H/'frozen_inputs.json').read_text());results=[]
with tempfile.TemporaryDirectory(prefix='metagraph-audit-') as tmp:
    t=Path(tmp)
    fake_atlas=t/'math/05-knowledge/results/merged_metagraph_black_quotient_bridges_20261020';fake_atlas.mkdir(parents=True)
    for name in pins['atlas']:shutil.copy2(repo/'05-knowledge/results/merged_metagraph_black_quotient_bridges_20261020'/name,fake_atlas/name)
    with (fake_atlas/'census_n3.json').open('ab') as f:f.write(b'\n')
    fake_producer=t/'producer';fake_producer.mkdir()
    for name in pins['producer']:shutil.copy2(producer/name,fake_producer/name)
    with (fake_producer/'PROOF.md').open('ab') as f:f.write(b'\n')
    for optimized in (False,True):
        for mode,argv,message in [('atlas',['--repo',str(t/'math'),'--producer',str(producer)],'Changed frozen atlas census_n3.json'),('producer',['--repo',str(repo),'--producer',str(fake_producer)],'Changed producer PROOF.md')]:
            p=subprocess.run([sys.executable]+(['-O'] if optimized else [])+[str(H/'independent_audit.py')]+argv,capture_output=True,text=True)
            if p.returncode==0 or message not in p.stderr:raise RuntimeError('Negative pin control failed: '+mode)
            results.append({'mutation':mode,'optimized':optimized,'returncode':p.returncode,'reason':message,'status':'correctly rejected'})
print(json.dumps({'status':'PASS','cases':results},indent=2))
