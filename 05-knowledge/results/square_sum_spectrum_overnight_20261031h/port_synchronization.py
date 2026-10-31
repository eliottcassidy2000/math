#!/usr/bin/env python3
"""Equal square class is enough after an explicitly charged duplication."""
from math import isqrt
import hashlib,json
from pythagorean_glue import Atlas,Closed,from_cube,glue,scale,edge,need,digest

def square_scale(A,z):
    f=z*z
    out=Closed({f*v for v in A.vertices},scale(A.base,f),
               [(scale(r,f),scale(a,f)) for r,a in A.moves],
               edge(f*A.port[0],f*A.port[1]),A.square_class,A.norm_seed)
    out.check_port();return out

def synchronize(A,B):
    A.check_port();B.check_port();need(A.square_class==B.square_class,'shared square class')
    if max(A.vertices)>=max(B.vertices):
        AA,ra=glue(A,A);BB,rb=glue(B,B,ra)
    else:
        BB,rb=glue(B,B);AA,ra=glue(A,A,rb)
    d=A.square_class
    a,b=isqrt(A.port[0]//d),isqrt(B.port[0]//d)
    left,right=square_scale(AA,b),square_scale(BB,a)
    need(left.port==right.port,'common new shape and exact integer port alignment')
    result,outer=glue(left,right)
    need(len(result.moves)==2*(len(A.moves)+len(B.moves)),'duplication pays both source histories')
    return result,{'source_ports':[A.port,B.port],'square_class':d,'first_frames':[ra,rb],
                  'alignment_square_roots':[b,a],'aligned_port':left.port,'outer':outer,
                  'source_multiplicities':[2,2]}

def main():
    a=Atlas();A=from_cube(a.cycle(2812,'A'));B=from_cube(a.cycle(2861,'A'),(196,2304))
    rejected=False
    try:glue(A,B)
    except RuntimeError as exc:rejected=str(exc)=='same normalized port type'
    need(rejected,'equal square class alone does not align the original two ports')
    out,record=synchronize(A,B);hashes=[]
    for salt in range(4):
        bits=''.join(str(hashlib.sha256(f'{salt}:{b}'.encode()).digest()[0]%2) for b in range(len(out.moves)))
        hashes.append(digest(out.walk(bits)))
    need(len(out.vertices)==11346 and len(out.moves)==100,'four exact source copies')
    print(json.dumps({'status':'passed','rejected_unaligned_direct_join':rejected,'synchronization':record,
          'full_cycle_sha256':hashes,'scope':'Different primitive port shapes in one square class become compatible after two copies of each source; counts square before multiplying. No unpaid one-copy heterogeneous join is asserted.'},indent=2,sort_keys=True))

if __name__=='__main__':main()
