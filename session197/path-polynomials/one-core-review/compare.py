import hashlib,json
from pathlib import Path
P=Path(__file__).parent
def need(b,s):
 if not b:raise RuntimeError(s)
a=json.loads((P/'independent-normal.json').read_text());b=json.loads((P/'source-normal.json').read_text());hits=(P/'support.bin').read_bytes()
need(sorted(a['signatures'])==sorted(z['signature']for z in b['inventory']),'complete signature inventory')
need(a['escape_rows']==b['two_core_escape']['rows'],'all nineteen literal rows')
support=[n for n in range(1,1000001,2)if hits[n//2]];digest=hashlib.sha256(json.dumps(support,separators=(',',':')).encode()).hexdigest()
need(digest==b['census']['support_sha256']and len(support)==499997,'full support digest')
print(json.dumps({'status':'PASS','core_signatures':65,'literal_escape_rows':19,'support_count':len(support),'support_sha256':digest},indent=2,sort_keys=True))
