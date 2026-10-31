"""Retain all family definitions and bind the regenerable complete finite trees."""
from pathlib import Path
import argparse,hashlib,json

def canonical(obj):
    return (json.dumps(obj,sort_keys=True,separators=(',',':'))+'\n').encode()
def compact_data(data,raw=None):
    if raw is None: raw=canonical(data)
    summary={k:v for k,v in data.items() if k!='families'}
    summary['complete_ledger_sha256']=hashlib.sha256(raw).hexdigest()
    summary['complete_ledger_bytes']=len(raw)
    summary['complete_second_joint_branches']=sum(len(b['finite_branches']) for f in data['families'] for b in f['branches'])
    summary['families']=[]
    for f in data['families']:
        branches=f['branches']
        summary['families'].append({
            'core':f['core'],'safe_intervals':f['safe_intervals'],
            'width':f['width'],'triple_cutoff':f['triple_cutoff'],
            'first_joint_count':len(branches),
            'second_joint_count':sum(len(b['finite_branches']) for b in branches),
            'finite_query_count':sum(max(0,row[2]-row[0]) for b in branches for row in b['finite_branches']),
            'all_rejection_lists_empty':all(not row[3] for b in branches for row in b['finite_branches']),
            'complete_branch_tree_sha256':hashlib.sha256(canonical(branches)).hexdigest()})
    summary['scope']='Hashes bind reproducible complete trees; they are not substitutes for replay or the analytic stopping proof.'
    return summary
def compact(source):
    raw=source.read_bytes()
    return compact_data(json.loads(raw),raw)

if __name__=='__main__':
    ap=argparse.ArgumentParser(); ap.add_argument('source',type=Path); ap.add_argument('output',type=Path); args=ap.parse_args()
    data=compact(args.source); args.output.write_bytes(canonical(data))
    print(json.dumps({'compact_bytes':args.output.stat().st_size,'complete_ledger_sha256':data['complete_ledger_sha256']}))
