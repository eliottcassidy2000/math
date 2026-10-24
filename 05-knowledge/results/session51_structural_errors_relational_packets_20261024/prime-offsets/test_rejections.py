#!/usr/bin/env python3
"""Hostile mutations must be rejected by the standalone verifier."""
import contextlib, copy, io, json, tempfile
from pathlib import Path
import check
original_root=check.ROOT
original=json.loads((original_root/'certificates.json').read_text())
mutations={}
a=copy.deepcopy(original); a['rows'].pop(); mutations['missing cutoff row']=a
a=copy.deepcopy(original); a['rows'][5]['p']=17; mutations['wrong primorial index']=a
a=copy.deepcopy(original); a['rows'][5]['plus']['witnesses'].pop(0); mutations['omitted smaller offset']=a
a=copy.deepcopy(original); a['rows'][5]['plus']['witnesses'][0]={'prime':True}; mutations['composite endpoint called prime']=a
a=copy.deepcopy(original); a['primality_certificates']['30047']['base']=1; mutations['invalid Lucas base']=a
a=copy.deepcopy(original); a['primality_certificates']['30047']['factors'][0][1]+=1; mutations['false n minus one factorization']=a
a=copy.deepcopy(original); a['rows'][0]['plus']['witnesses'][0]={'divisor':4}; mutations['divisor equals composite target']=a
for name,obj in mutations.items():
    with tempfile.TemporaryDirectory() as d:
        check.ROOT=Path(d)
        (check.ROOT/'certificates.json').write_text(json.dumps(obj))
        try:
            with contextlib.redirect_stdout(io.StringIO()): check.verify()
        except ValueError:
            print('PASS rejected:',name)
        else:
            raise RuntimeError('FALSE CERTIFICATE ACCEPTED: '+name)
check.ROOT=original_root
print('PASS: all 7 hostile certificate mutations rejected')
