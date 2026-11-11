#!/usr/bin/env python3
"""Construct or exclude an exact multiplier q, at every base order.

python3 construct_seed.py --q 20 --tail-length 1
The complete finite search may be enormous for large q. It does not decide
whether an odd slope c succeeds at some unspecified exponent.
"""
from pathlib import Path
import argparse,json,sys
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'seed_decision'))
from check import need,tournament,strong,endpoint_tables,block_data,observables,join,bound

def ispow2(x):return x>0 and x&(x-1)==0

def construct(q,t=1):
    if q<1 or t<1:raise ValueError('q and tail length must be positive')
    R=bound(q//2) if q>1 else 1
    found=None;tested=0
    if ispow2(q):found=((0,),1,1)
    else:
        for n in range(3,R+1):
            for code in range(1<<(n*(n-1)//2)):
                T=tournament(n,code)
                if not strong(T):continue
                h=sum(endpoint_tables(T)[-1])
                if h>q//2:continue
                for O in range(1<<n):
                    if O.bit_count()<2 or O.bit_count()*h>q:continue
                    tested+=1;hh,w,u=observables(block_data(T,O))
                    if hh-w==1 and u and q%u==0 and ispow2(q//u):
                        found=T,O,u;break
                if found:break
            if found:break
    if not found:return {'status':'NO_SEED','exact_multiplier':q,'scope':'All base orders within the selected-set tail mechanism; not all exponents of the odd slope','strong_core_order_bound':R,'selected_core_presentations_tested':tested}
    S,O,u=found;ratio=q//u;b=ratio.bit_length()-1;B=S
    for _ in range(b):B=join(B,(0,))
    h,w,U=observables(block_data(B,O));need(h-w==1 and U==q,'literal base seed')
    T=B
    for _ in range(t):T=join(T,(0,))
    rows=list(T);v=len(T)
    for a in range(v):
        if not(O>>a&1):rows[a]|=1<<v
    rows.append(O);actual=sum(endpoint_tables(tuple(rows))[-1]);expected=1+(q<<t)
    need(actual==expected,'literal emitted tournament')
    c=q;j=0
    while c%2==0:c//=2;j+=1
    return {'status':'SEED','exact_multiplier':q,'odd_slope':c,'multiplier_exponent':j,'ray_start_exponent':j+1,'core_rows':list(S),'selected_mask':O,'sink_padding':b,'base_rows':list(B),'H':h,'W':w,'U':U,'tail_length':t,'tournament_rows':rows,'expected_H':expected,'counted_H':actual,'proof':'seed_decision/CORE_CLASSIFICATION.md','selected_core_presentations_tested':tested}

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--q',type=int);p.add_argument('--tail-length',type=int,default=1);p.add_argument('--self-test',action='store_true');a=p.parse_args()
    if a.self_test:
        cases=[1,3,5,6,7,9,10,11,13,20,29];yes={1,6,11,20,29};rows=[]
        for q in cases:
            r=construct(q);need((r['status']=='SEED')==(q in yes),'complete constructor control');rows.append(r)
        out={'status':'PASS','scope':'Exact-q CLI controls, positive literal adjacency and negative all-order core exhaustion','cases':rows}
        (ROOT/'CONSTRUCTOR_RECEIPT.json').write_text(json.dumps(out,indent=2)+'\n');print('PASS: fixed-multiplier constructor, 11 exact-q controls')
    else:
        if a.q is None:p.error('--q is required unless --self-test is used')
        print(json.dumps(construct(a.q,a.tail_length),indent=2))
if __name__=='__main__':main()
