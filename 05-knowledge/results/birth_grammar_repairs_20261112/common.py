from pathlib import Path
import sys,gzip,json,hashlib
HERE=Path(__file__).resolve().parent
OLD=HERE.parent/'joint_stability_forbidden_20261112'
sys.path.insert(0,str(OLD))
from census_nine import need,tables,ib,strong,defect,automorphisms,move
from stable_census import selected_defects,certify
from joint_nine import quotient,key
import numpy as np

def read(name,old=False):
 return json.loads(gzip.decompress(((OLD if old else HERE)/name).read_bytes()))
def save(name,data):
 if not __debug__:name=name.replace('.json.gz','.optimized.json.gz')
 (HERE/name).write_bytes(gzip.compress((json.dumps(data,separators=(',',':'))+'\n').encode(),mtime=0))
def expcoeff(u):
 from math import comb
 return [sum(u[m-1]*(-1)**(m-j)*comb(m,j) for m in range(j,len(u)+1)) for j in range(1,len(u)+1)]
def value(e,a):return sum(c*(i+1)**a for i,c in enumerate(e))
def v2(n):
 need(n!=0,'finite valuation');return (abs(n)&-abs(n)).bit_length()-1
