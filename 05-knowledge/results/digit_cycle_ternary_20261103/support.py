"""Inherited marked tournament interfaces, plus deterministic local receipts."""
from pathlib import Path
import sys,json
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
PREV=HERE.parent/'return701_algebra_20261103'
sys.path.insert(0,str(PREV))
from common import junction,need,observer,values
from return_native import paths,flip,phases
def save(name,obj):(HERE/name).write_text(json.dumps(obj,indent=2)+'\n',encoding='utf8',newline='\n')
def source():
    d=json.loads((PREV/'native_return.json').read_text())
    T=tuple(d['next_rows']);names=d['next_names'];r=d['next_root'];b=names.index(6)
    need(junction.H(T)==2107 and names[r]==21,'actual inherited 2107 root and source')
    return T,r,b,names
