"""Actual saved sources and pinned inherited interfaces for the 18979 session."""
from pathlib import Path
import json,sys
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2];PREV=HERE.parent/'return6325_observer_20261103'
sys.path.insert(0,str(PREV));sys.path.insert(0,str(HERE.parent/'return701_algebra_20261103'))
from common import junction,need,observer,values
from return_native import paths,flip,phases
def save(name,d):(HERE/name).write_text(json.dumps(d,indent=2)+'\n',encoding='utf8',newline='\n')
def source():
    d=json.loads((PREV/'return6325.json').read_text());T=tuple(d['next_rows']);r=d['next_root'];names=d['next_names'];b=names.index(6)
    need(junction.profile(T,r)==(18979,6325,6325,18979) and names[r]==43,'actual18979 source and old root')
    return T,r,b,names
def barrier_source():
    d=json.loads((PREV/'return6325.json').read_text());return tuple(d['final_rows']),d['root'],d['spare'],d['names']
