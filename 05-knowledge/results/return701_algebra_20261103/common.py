"""Current-session imports; inherited sources are pinned by the replay driver."""
from pathlib import Path
import sys,json
import sympy as s
from sympy.polys.matrices import DomainMatrix
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
PREV=HERE.parent/'joint_port_geometry_20261103';sys.path.insert(0,str(PREV))
from return_probe import old,junction,need,observer,values
import effects as e

def save(name,obj):(HERE/name).write_text(json.dumps(obj,indent=2)+'\n',encoding='utf8',newline='\n')
def data(M):return e.data(M)
def dm(M):return DomainMatrix.from_Matrix(M).to_field()
def rref(M):
    R,P=dm(M).rref();return R.to_Matrix()[:len(P),:],P
def nullrows(M):return dm(M).nullspace().to_Matrix()
def read(name):return json.loads((PREV/name).read_text())

def source701():
    d=read('return_native.json');T=tuple(d['final_rows']);names=d['final_names'];r=names.index(7);b=names.index(6)
    U=junction.ear(junction.ear(T,1<<r),1<<len(T));names=names+[13,14]
    need(junction.profile(U,len(T))==(701,233,233,701),'exact saved source701')
    return U,len(T),b,names
