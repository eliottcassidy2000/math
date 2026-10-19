from pathlib import Path
import json,hashlib,importlib.util,time
ROOT=Path(__file__).resolve().parent
sp=importlib.util.spec_from_file_location('S',ROOT/'scan_c2.py');S=importlib.util.module_from_spec(sp);sp.loader.exec_module(S)
data=S.prepare();ids=[i for i,item in enumerate(data['cores']) if len(item)==1]
(ROOT/'residual_certificates.jsonl.gz').write_bytes(b'')
S.write(ROOT/'progress.json',{'status':'PAUSED_FOR_SHARDS','completed_residual_cores':0,'last_core_index':-1,'committed_bytes':0})
plan={'prefix_completed':0,'prefix_last_core_index':-1,'prefix_record_sha256':hashlib.sha256(b'').hexdigest(),'shards':[[i for pos,i in enumerate(ids) if pos%4==j] for j in range(4)],'cover_sha256':S.sha(ROOT/'analytic_cover.json.gz')}
(ROOT/'shard_plan.json').write_text(json.dumps(plan,separators=(',',':'))+'\n');print('SHARDS',list(map(len,plan['shards'])),flush=True)
