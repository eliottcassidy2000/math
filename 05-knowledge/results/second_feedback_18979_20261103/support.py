"""Source-preserving imports for the second-feedback research session."""
from pathlib import Path
import json,sys
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
PREV=HERE.parent/'return18979_ports_20261103'
INCOMING=HERE.parent/'return18979_palindrome_20261103'
sys.path.insert(0,str(HERE.parent/'return701_algebra_20261103'))
from common import junction,need,observer,values
from return_native import phases,flip,paths
sys.path.insert(0,str(HERE.parent/'return6325_observer_20261103'))
def save(name,d):(HERE/name).write_text(json.dumps(d,indent=2)+'\n',encoding='utf8',newline='\n')
def source():
    d=json.loads((HERE.parent/'return6325_observer_20261103/return6325.json').read_text())
    return tuple(d['next_rows']),d['next_root'],d['next_names'].index(6),d['next_names']
from matching_carry import clone,matching
