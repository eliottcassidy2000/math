#!/usr/bin/env python3
"""Reconcile frozen inventory counts and blob selection, not all lexical rereading."""
import json,gzip,pathlib,collections,subprocess,re,hashlib
base=pathlib.Path(__file__).resolve().parents[1];p=base/'history'
def get(n):return json.loads(gzip.decompress((p/n).read_bytes()) if n.endswith('.gz') else (p/n).read_text())
def req(v,msg):
 if not v:raise RuntimeError(msg)
S=get('summary.json');I=get('inventory.json.gz');E=get('expansion.json.gz');M=get('searched_manifest.json.gz');U=get('unavailable.json')
req(len(M)==S['read_paths']==34864,'read count')
req(len(I)==S['seed_hit_paths']==7509,'seed count')
req(len({x['blob'] for x in I})==S['seed_hit_distinct_blobs']==6877,'distinct seeds')
req(len(E)==S['expansion_only_paths_with_at_least_two_topics']==10069,'expansion count')
req(len(U)==S['unavailable_paths']==23,'omissions')
req(set(x['path'] for x in I).isdisjoint(x['path'] for x in E),'seed/expansion partition')
counts=collections.Counter(k for x in I for k in x['queries'])
req(all(counts[k]==v for k,v in S['query_path_counts'].items()),'query counts')
ms={x['path']:x['blob'] for x in M};req(all(ms[x['path']]==x['blob'] for x in I+E),'manifest membership')
root=subprocess.check_output(['git','rev-parse','--show-toplevel'],cwd=p,text=True).strip()
raw=subprocess.check_output(['git','ls-tree','-rz',S['ref']],cwd=root);candidate={}
for x in raw.split(b'\0'):
 if not x:continue
 a,b=x.split(b'\t',1); mode,kind,sha=a.decode().split();name=b.decode()
 if kind=='blob' and re.match(r'0[0-7]-',name) and pathlib.Path(name).suffix.lower() in {'.md','.txt','.tex','.rst'}:candidate[name]=sha
req(len(candidate)==S['candidate_paths']==34887,'candidate count')
req(all(candidate[k]==v for k,v in ms.items()),'pinned objects')
req(set(ms)|{x['path'] for x in U}==set(candidate),'complete scope partition')
for name,h in get('catalog_hashes.json').items():
 data=(p/(name+'.gz')).read_bytes();raw=gzip.decompress(data)
 req(hashlib.sha256(data).hexdigest()==h['gzip_sha256'],'gzip hash')
 req(hashlib.sha256(raw).hexdigest()==h['raw_sha256'],'raw hash')
 req(len(data)==h['gzip_bytes'] and len(raw)==h['raw_bytes'],'byte counts')
print(json.dumps({'status':'PASS','candidate_paths':len(candidate),'read_paths':len(M),'seed_paths':len(I),'expansion_paths':len(E),'omitted_paths':len(U),'scope':'manifest and aggregate reconciliation; no full lexical replay'},indent=2,sort_keys=True))
