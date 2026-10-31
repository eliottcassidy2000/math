"""Bind already independently reconstructed records to a producer receipt.

Run rebuild.py first. For a single combined audit use check.py instead.
This split permits independent reconstruction to overlap producer execution.
"""
import sys,json,hashlib
from pathlib import Path
from check import validate_header,need,verify_family
raw=Path(sys.argv[1]).read_bytes();d=json.loads(raw);validate_header(d)
records=json.loads(Path(sys.argv[2]).read_bytes())
need(len(records)==1287,'independent rebuild complete')
totals=[0]*6;cuts=[0]*6;inherited=tree_bytes=0
for rec,expected in zip(d['families'],records):
    verify_family(rec,expected)
    totals=[a+b for a,b in zip(totals,expected['node_counts'])]
    cuts=[max(a,b) for a,b in zip(cuts,expected['max_cutoffs'])]
    inherited+=expected['inherited_branches'];tree_bytes+=expected['tree_bytes']
need(totals==d['node_counts'] and cuts==d['max_cutoffs'] and inherited==d['inherited_branches'],'aggregate counters')
print(json.dumps({'status':'PASS','core_count':1287,'node_counts':totals,'max_cutoffs':cuts,'inherited_branches':inherited,'tree_commitments_checked':1287,'full_tree_bytes_reconstructed':tree_bytes,'receipt_sha256':hashlib.sha256(raw).hexdigest()},sort_keys=True,indent=2))
