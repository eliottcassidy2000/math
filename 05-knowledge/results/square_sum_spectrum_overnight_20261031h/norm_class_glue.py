#!/usr/bin/env python3
"""Nontrivial square-class ports and an explicitly paid variable port."""
import hashlib,json
from pythagorean_glue import Atlas,from_cube,condition_cube,glue,need,digest

def main():
    atlas=Atlas();records=[]
    for N,port,d in [(2812,(20,80),5),(2828,(578,1922),2),(2812,(160,1440),10)]:
        c=atlas.cycle(N,'A');A,sidecar=condition_cube(c,port,d)
        C,rec=glue(A,A)
        samples=[]
        for salt in range(4):
            bits=''.join(str(hashlib.sha256(f'{salt}:{b}'.encode()).digest()[0]%2) for b in range(len(C.moves)))
            samples.append(digest(C.walk(bits)))
        D,rec2=glue(C,C)
        bits='01'*(len(D.moves)//2)+'0'*(len(D.moves)%2)
        second=digest(D.walk(bits))
        need(D.norm_seed==A.norm_seed and D.square_class==d,'fixed norm seed under recursion')
        records.append({'N':N,'port':port,'class':d,'source_context':sidecar,'first':rec,
                        'first_cycle_sha256':samples,'second':rec2,'second_cycle_sha256':second})
    need(records[-1]['source_context']['pins']=={22:1} and records[-1]['first']['output_bits']==48,'pay the variable port in both copies')
    # A tempting forgotten terminal-state guard fails at the first tested hostile.
    c=atlas.cycle(2815,'A')
    need(c.compile([((9,16),True)]) is None,'9--16 is not universal across terminal states')
    print(json.dumps({'status':'passed','records':records,'missing_port_hostile_N':2815,
          'scope':'Exact selected-carrier products in square classes2,5,10, retaining the source-bit pin for the variable class10 port. No full induced-graph count is claimed.'},indent=2,sort_keys=True))

if __name__=='__main__':main()
