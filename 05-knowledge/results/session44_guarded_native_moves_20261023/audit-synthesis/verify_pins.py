"""Check the accepted synthesis/source snapshot in research or published layout."""
from pathlib import Path
import argparse,hashlib,json
p=argparse.ArgumentParser();p.add_argument('--workspace',type=Path,default=Path(__file__).resolve().parent.parent);a=p.parse_args()
W=a.workspace;report=json.loads((Path(__file__).resolve().parent/'SOURCE_PINS.json').read_text());count=0
for item in report['files']:
 path=W/item['research_path']
 if not path.exists():path=W/item['logical_path']
 actual=hashlib.sha256(path.read_bytes()).hexdigest()
 if actual!=item['sha256']:raise RuntimeError('Source snapshot changed: '+item['logical_path'])
 count+=1
print(json.dumps({'status':'PASS','frozen_source_files':count,'synthesis_sha256':report['files'][0]['sha256']},indent=2))
