"""Pinned inherited interfaces; this packet starts from the actual 6325 child."""
from pathlib import Path
import sys,json
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
PREV=HERE.parent/'digit_cycle_ternary_20261103'
sys.path.insert(0,str(PREV));sys.path.insert(0,str(HERE.parent/'return701_algebra_20261103'))
from common import junction,need,observer,values
from return_native import paths,flip,phases
def save(name,obj):(HERE/name).write_text(json.dumps(obj,indent=2)+'\n',encoding='utf8',newline='\n')
def source():
    d=json.loads((PREV/'return2107.json').read_text());T=tuple(d['next_rows']);names=d['next_names']
    r=d['next_root'];b=names.index(6)
    need(junction.profile(T,r)==(6325,2107,2107,6325) and names[r]==31,'actual 6325 source')
    return T,r,b,names
