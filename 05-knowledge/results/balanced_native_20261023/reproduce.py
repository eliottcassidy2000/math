#!/usr/bin/env python3
"""Normal/-O exact replay and frozen local-file manifest. No search output is trusted."""
from pathlib import Path
import argparse,hashlib,json,subprocess,sys
HERE=Path(__file__).resolve().parent
RUNS=['one_bit','charts','flow','zero_sum','search28','native18','audit','route18','route18_audit','multi_target','product_flow','four_move','six_trade']
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--check',action='store_true');args=ap.parse_args()
    for stem in RUNS:
        results=[]
        for mode,flags in [('normal',[]),('optimized',['-O'])]:
            raw=subprocess.check_output([sys.executable,*flags,'-B',str(HERE/(stem+'.py'))])
            data=json.loads(raw);results.append(data)
            path=HERE/(stem+'.'+mode+'.json')
            if args.check:
                if json.loads(path.read_text())!=data:raise RuntimeError('receipt mismatch: '+path.name)
            else:path.write_text(json.dumps(data,indent=2,sort_keys=True)+'\n')
        if results[0]!=results[1]:raise RuntimeError('optimized drift: '+stem)
        print(stem,'PASS',results[0].get('checks','bounded search'),flush=True)
    files={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(HERE.iterdir()) if p.is_file() and p.name!='MANIFEST.json'}
    manifest={'status':'normal and optimized exact replay PASS','input_revision':'0f9dec9c3e1d',
      'incoming_domain_revision':'66db753a59','incoming_compensator_revision':'a793424e20','files_sha256':files}
    path=HERE/'MANIFEST.json'
    if args.check:
        if json.loads(path.read_text())!=manifest:raise RuntimeError('local-file manifest drift')
    else:path.write_text(json.dumps(manifest,indent=2,sort_keys=True)+'\n')
if __name__=='__main__':main()
