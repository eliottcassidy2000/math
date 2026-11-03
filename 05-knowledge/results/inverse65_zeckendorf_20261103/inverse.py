"""Exact guarded inverse; terminates even when a requested state is absent."""
from functools import lru_cache
from pathlib import Path
import argparse,json
import sympy as sp
from bank import MATS,COSTS,need,step,vector,closure,recover,actual_target
HERE=Path(__file__).resolve().parent
ROOT=(1,0,0)
INVERSES={}
for label,rows in MATS.items():
    matrix=sp.Matrix(rows);det=int(matrix.det())
    need(det!=0,('singular atom',label))
    INVERSES[label]=(det,tuple(tuple(int(v) for v in row) for row in matrix.adjugate().tolist()))

def norm(v):return 6*v[0]+9*v[1]+5*v[2]

def physical_cone(v):
    if v==ROOT:return True
    h,x,y=v
    return h>=1 and h%2==1 and 0<=x<=y<=h+x and h<=x+y and x<=16*h and y<=16*h

def predecessors(v):
    for label,(det,adj) in INVERSES.items():
        nums=step(v,adj)
        if any(z%det for z in nums):continue
        parent=tuple(z//det for z in nums)
        if not physical_cone(parent) or parent[0]>v[0] or 9*norm(parent)>5*norm(v):continue
        need(step(parent,MATS[label])==v,('inverse arithmetic',label,v))
        yield parent,label

@lru_cache(None)
def decode_state(v):
    if v==ROOT:return ()
    if not physical_cone(v):return None
    options=sorted(predecessors(v),key=lambda z:(norm(z[0]),COSTS[z[1]],z[1]))
    for parent,label in options:
        prefix=decode_state(parent)
        if prefix is not None:return prefix+(label,)
    return None

def solve_h(h):
    """Uniform complete scalar decision procedure; no surjectivity assumption."""
    if h<1 or h%2==0:return None
    if h==1:return ('a',)  # scalar witnesses are nonempty actual tournaments
    states,parents,_=closure(h)
    for v in states:
        if v[0]!=h:continue
        word=decode_state(v)
        need(word is not None and vector(word)==v,('forward/backward agreement',h,v))
        return word
    return None

def source_certificate(word):
    v=ROOT;n=0;clock=[]
    for label in word:
        v=step(v,MATS[label]);n+=COSTS[label]
        clock.append({'letter':label,'vertices':n,'state':v,'port_in_processing_frame':n-1})
    return {'source_type':'actual_tournament' if word else 'formal_empty_unit',
            'processing_word':word,'source_order_word':tuple(reversed(word)),
            'clock':clock,'vertices':n,'H':v[0],'state':v,'target_rows':actual_target(word)}

def main():
    parser=argparse.ArgumentParser();parser.add_argument('H',type=int);parser.add_argument('--x',type=int);parser.add_argument('--y',type=int)
    args=parser.parse_args()
    if (args.x is None)!=(args.y is None):parser.error('supply both --x and --y')
    word=solve_h(args.H) if args.x is None else decode_state((args.H,args.x,args.y))
    print(json.dumps({'status':'ABSENT_IN_DECLARED_BANK'} if word is None else {'status':'CONSTRUCTED',**source_certificate(word)}))
if __name__=='__main__':main()
