"""Coordinatewise mod4 stabilization under arbitrary transitive blowups."""
from pathlib import Path
import json
from frame_probe import ordinary_bank,hp,valid,need

HERE=Path(__file__).resolve().parent

def blowup(T,lengths):
    labels=[v for v,m in enumerate(lengths) for _ in range(m)]
    return tuple(sum(int(i!=j and ((a==b and i<j) or (a!=b and T[a]>>b&1)))<<j for j,b in enumerate(labels)) for i,a in enumerate(labels))

def run():
    controls=0;signatures=[]
    for (n,code),T in sorted(ordinary_bank().items()):
        if n>5:continue
        signature=[]
        for S in range(1<<n):
            lengths=[2 if S>>i&1 else 1 for i in range(n)]
            base=blowup(T,lengths);need(valid(base),'literal transitive blowup')
            h=hp(base);signature.append(h%4)
            for i in range(n):
                if not(S>>i&1):continue
                longer=lengths.copy();longer[i]=3
                B=blowup(T,longer);need(valid(B) and hp(B)%4==h%4,'coordinatewise mod4 stabilization')
                controls+=1
        signatures.append(dict(n=n,ordinary_code=code,rows=T,signature=signature))
    cyclic=[]
    for m in range(1,5):
        T=blowup((2,4,1),(m,m,m));h=hp(T)
        cyclic.append([m,h,h%4])
        if m>=2:need(h%4==cyclic[1][2],'cyclic port pattern also stabilizes')
    halvings=[(9*3**m+9*2**m-6)%4 for m in range(2,12)]
    need(halvings==[3,1]*5,'halved curve does not stabilize modulo4')
    result=dict(status='FINITE-EXACT controls; all-order stabilization proved in CLONE_PARITY.md',scope='All ordinary frames through5; every subset of doubled ports; each doubled port extended to3 while all other lengths retained. Also cyclic three-port blowups at widths1..4.',coordinate_controls=controls,signatures=signatures,cyclic_three_port=cyclic,halved_curve_residues=halvings)
    suffix='' if __debug__ else '.optimized'
    (HERE/('clone_parity'+suffix+'.json')).write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print('PARITY_CONTROLS',controls,'FRAMES',len(signatures),'CYCLIC_PORTS',cyclic,flush=True)

if __name__=='__main__':run()
