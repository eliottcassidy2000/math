"""Pinned inherited exact engines; keep the old common module unshadowed."""
from pathlib import Path
import sys,gzip,json,hashlib
HERE=Path(__file__).resolve().parent
BIRTH=HERE.parent/'birth_grammar_repairs_20261112'
JOINT=HERE.parent/'joint_stability_forbidden_20261112'
sys.path.insert(0,str(BIRTH))
from common import need,tables,ib,strong,defect,automorphisms,move,np,value,v2,expcoeff
from hamilton_clocks import skeleton,trim

def load(name,bank=BIRTH):
    path=bank/name
    return json.loads(gzip.decompress(path.read_bytes()) if path.suffix=='.gz' else path.read_text(encoding='utf8'))

def save(name,data):
    if not __debug__:name=name.replace('.json.gz','.optimized.json.gz')
    (HERE/name).write_bytes(gzip.compress((json.dumps(data,separators=(',',':'))+'\n').encode(),mtime=0))

