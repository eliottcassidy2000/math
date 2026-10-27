#!/usr/bin/env python3
"""Check completed replay receipts and their source/scope bindings, not rerun jobs."""
import argparse,hashlib,json
from pathlib import Path
from c3_transport import need,read_regular,sha

def check(packet):
 q='verification.checkpoint3/'
 def load(n):return json.loads(read_regular(packet,q+n))
 driver=sha(read_regular(packet,'c3_replay.py'));modes={};receipts=[]
 observed=load('MATERIALIZED_INPUTS.final.json');transport=json.loads(read_regular(packet,'TRANSPORT.checkpoint3.json'))
 actual={r['path']:(r['sha256'],r['bytes'],r['mode']) for r in observed}
 expected={r['logical_path']:(r['logical_sha256'],r['logical_bytes'],int(r['logical_mode'],8)) for r in transport['files']}
 need(len(observed)==len(actual)==600 and actual==expected,'executed materialized inputs differ from transported original bytes/modes')
 for mode in ('normal','optimized'):
  results={}
  for part in ('basic','joint','joint-audit'):
   name='final-'+part+'-'+mode+'.json';j=load(name)
   need(j['status']=='PASS' and j['mode']==mode and j['adapter_sha256']==driver,'replay status/mode/driver mismatch')
   for k,v in j['results'].items():
    need(k not in results and v['status']=='PASS','duplicate or failed replay component');results[k]=v
   need(j['jobs'] and all(r['returncode']==0 and r['stderr_sha256']==sha(b'') for r in j['jobs']),'failed/warning replay command')
   need(j['negative_adapter_controls']['status']=='PASS','adapter hostile controls failed')
   receipts.append({'path':q+name,'sha256':sha(read_regular(packet,q+name))})
  expected={'boundary','boundary-audit','twins','twins-audit','root-twins','joint','joint-audit','tc','tc-window'}
  need(set(results)==expected,'incomplete replay components')
  need(results['joint']['scales']==19 and results['joint']['producer_rectangles']==245826 and results['joint']['source_records']==1000908,'producer scope count')
  need(results['joint-audit']['rectangles_independently_recomputed']==42425 and results['joint-audit']['direct_controls']==43 and results['joint-audit']['source_records_rebuilt']==1000908,'independent audit coverage')
  need(results['tc-window']['window_jobs']==83 and results['tc-window']['window_strong']==185642 and results['tc-window']['alternate_n12_jobs']==45 and results['tc-window']['alternate_n12_strong']==46968,'window coverage')
  need(results['tc-window']['full_window_replayed'] and not results['tc-window']['original_large_n12_scan_replayed'],'window scope')
  modes[mode]=results
 need(modes['normal']==modes['optimized'],'normal/optimized semantic results differ')
 u=load('final-ubsan.json');need(u['status']=='PASS' and u['mode']=='ubsan' and u['adapter_sha256']==driver,'UBSan receipt/source')
 need(set(u['results'])=={'boundary','boundary-audit','tc','tc-window'},'UBSan scope')
 need(all(r['returncode']==0 and r['stderr_sha256']==sha(b'') for r in u['jobs']),'UBSan warning/failure')
 for stem in ('transport','provenance','corruption','provenance-corruption','git-origins'):
  a=read_regular(packet,q+stem+'.normal.json');b=read_regular(packet,q+stem+'.optimized.json')
  need(a==b and json.loads(a)['status']=='PASS','normal/optimized verification receipt mismatch '+stem)
 cpp=read_regular(packet,q+'CPP_REJECTIONS.normal.json')
 need(cpp==read_regular(packet,q+'CPP_REJECTIONS.optimized.json') and json.loads(cpp)['status']=='PASS' and json.loads(cpp)['rejections']==30,'C++ rejection controls')
 inherited=load('inherited-replay.normal.json');need(inherited['status']=='PASS' and len(inherited['replays'])==16,'inherited C1/C2 replay')
 binding=load('inherited-execution-binding.json');snapshot=read_regular(packet,q+'inherited-execution-manifest.json')
 need(sha(snapshot)==binding['execution_manifest_sha256'] and sha(read_regular(packet,q+'inherited-replay.normal.json'))==binding['receipt_sha256'],'inherited execution binding')
 oldrows={r['path']:r for r in json.loads(snapshot)['files']}
 need(len(oldrows)==inherited['files_verified']==binding['files_verified_at_execution'],'inherited execution inventory count')
 prior=json.loads(read_regular(packet,'FILE_MANIFEST.checkpoint2.json'))
 oldnames=[r['path'] for r in prior['files']]+['FILE_MANIFEST.checkpoint2.json']
 need(len(oldnames)==206,'C1/C2 preservation count')
 for n in oldnames:
  b=read_regular(packet,n);r=oldrows[n]
  need(len(b)==r['bytes'] and sha(b)==r['sha256'],'C1/C2 replay input changed since execution')
 return dict(status='PASS',driver_sha256=driver,normal_and_optimized_components=9,producer_rectangles=245826,independently_recomputed_rectangles=42425,independent_direct_controls=43,independently_rebuilt_source_records=1000908,tc_window_jobs=83,tc_alternate_order12_jobs=45,original_large_order12_scan='NOT_RERUN',ubsan='Full n4..9 universe authentication and declared Cayley/Aut controls; selected n22 module/n16 nonmodule scan jobs, not full window',inherited_C1_C2_replay_entries=16,receipts=receipts)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--packet',type=Path,default=Path(__file__).resolve().parent);a=p.parse_args();print(json.dumps(check(a.packet),sort_keys=True,indent=2))
