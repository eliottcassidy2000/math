#!/usr/bin/env python3
"""Hostile mutation controls against this auditor, not producer assertions."""
from copy import deepcopy
from pathlib import Path
import json
import sys
from check_rows import audit

def main():
    path=Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).resolve().parent.parent/'ramsey'/'rows.normal.json'
    original=json.loads(path.read_text())
    cases=[]
    z=deepcopy(original);z['certificates'].pop();cases.append(('missing_row',z))
    z=deepcopy(original);z['certificates'][0]['sharp_patterns'].pop();cases.append(('missing_clique_pattern',z))
    z=deepcopy(original);z['certificates'][0]['sharp_patterns'][0]['dual'][0]='0';cases.append(('false_dual',z))
    z=deepcopy(original);z['certificates'][0]['cutoff']-=1;cases.append(('false_cutoff',z))
    z=deepcopy(original);z['certificates'][0]['finite_upper'][0][1]+=1;cases.append(('false_finite_root',z))
    z=deepcopy(original);z['certificates'][0]['seeds'][0]['n']+=1;cases.append(('false_seed_order',z))
    rejections=[]
    for name,receipt in cases:
        try:audit(receipt)
        except RuntimeError as err:rejections.append({'mutation':name,'rejected':str(err)})
        else:raise RuntimeError('mutation unexpectedly accepted: '+name)
    print(json.dumps({'status':'PASS','rejections':rejections},sort_keys=True,indent=2))

if __name__=='__main__':main()
