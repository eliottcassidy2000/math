"""Literal and affine controls for ENDCAP_NEXT.md, including the module carrier."""
import importlib.util
import json
from pathlib import Path


spec = importlib.util.spec_from_file_location('main_ports', Path(__file__).with_name('check.py'))
main = importlib.util.module_from_spec(spec)
spec.loader.exec_module(main)
base, need = main.base, main.need


def replacement(n):
    w, r = tuple(range(4,n-3)), tuple(range(n-3,n))
    return [(0,1)+w+r+(2,3), (1,)+w+(0,)+r+(2,3), (1,2,3)+w+(0,)+r]


def check_carrier(n):
    src = base.source(n,{3,n-3,n-2},[(3,1)])
    p = tuple(range(n))
    bflips = set(src.edges)^base.pathset(p)
    d3 = main.words(n,3,False)[2]
    src.anti(d3,bflips)
    need(src.ok,'two-C3 carrier exists')
    def bforce(u,v,value):
        need(src.forced(u,v,value^(base.pair(u,v) in bflips)),('carrier arc',n,u,v,value))
    left, right = (0,1,2), (n-3,n-2,n-1)
    for a in left:
        for b in right:
            bforce(a,b,1)
        for v in range(3,n-3):
            bforce(a,v,int(v==3))
    for a in right:
        for v in range(3,n-3):
            bforce(a,v,int(v!=n-4))
    for triple in (left,right):
        for a,b in zip(triple,triple[1:]+triple[:1]):
            bforce(a,b,1)
    need(src.forced(1,4,1),'carrier forces repair chord')
    hostile=src.clone();hostile.force(4,1)
    need(not hostile.ok,'opposite chord violates carrier')
    qs=replacement(n)
    main.old.certify(src,qs)


def main_check():
    orders=list(range(12,61))+[80,100,140]
    cells=0
    for n in orders:
        for chord in ((1,4),(4,1)):
            qs=replacement(n) if chord==(1,4) else main.words(n,3,False)
            src,_=base.symbolic(n,{3,n-3,n-2},qs,[(3,1),chord],True)
            for seed in (0,373):
                base.literal(src.model(seed),qs)
            cells+=1
        check_carrier(n)
    print(json.dumps({'status':'PASS','complete_affine_chord_cells':cells,
        'orders':orders,'carrier_cells':len(orders),'literal_models_per_cell':2,
        'hostile':'D3-blue carrier forbids4->1',
        'scope':'endpoint-four h3,3->1; literal exact target; no quotient of the middle'},indent=2))


if __name__=='__main__':
    main_check()
