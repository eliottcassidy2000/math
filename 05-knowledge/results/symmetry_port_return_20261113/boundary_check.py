"""Exact controls for changing a uniform module boundary."""
from pathlib import Path
from math import gcd
import json
ROOT=Path(__file__).resolve().parent
def need(ok,why):
    if not ok:raise RuntimeError(why)
def wrapper(n):
    h,p,w,u=95**n,1276**n,8709**n,40760**n
    return h+p,2*w+4*p-2*h,3*u+9*w-3*p
def main():
    out=[]
    for n,expected in [(1,3),(2,1),(191,3)]:
        h,p,w=wrapper(n);g2=gcd(h,p);g3=gcd(g2,w)
        need(g3==expected,'three-observation content')
        if n==191:need(g2==1371 and w%457==312,'higher observation hostile')
        out.append(dict(n=n,H=h if n<3 else '95^191+1276^191',P2=p if n<3 else '2*8709^191+4*1276^191-2*95^191',gcd_H_P2=g2,gcd_H_P2_P3=g3,P3_mod457=w%457))
    h1,p1,_=wrapper(1);h2,p2,_=wrapper(2)
    G1=60*3*p1+h1;G2=60*9*p2+h2
    need((G1,G2)==(4021131,85423165841) and gcd(G1,G2)==1,'outer release')
    d=json.loads((ROOT/'pair_growth2.json').read_text())
    need(all(all(mask not in (0,255) for mask in h[side][:2]) for h in d['hits'] for side in ('left_attachment','right_attachment')),'old module boundaries are actually broken')
    receipt=dict(status='PASS',controls=out,outer_endpoints=[G1,G2],outer_gcd=1,all_determinant_extension_masks_nonconstant=True)
    (ROOT/'boundary_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print('PASS rigid context3→content1 release; higher-coordinate hostile; actual mixed incidences')
if __name__=='__main__':main()
