#!/usr/bin/env python3
"""Exact formal-signature control, not a tournament-realization census."""
from fractions import Fraction as Q
from functools import lru_cache
import json

@lru_cache(None)
def g(a):
    if a < 2:
        return 1
    q = a // 2
    return g(q) if a % 2 else g(q) + q + q % 2

def response(s,a,b):
    d,kp,km,ip,im,op,om,j=s
    return (d*(g(a)+g(b)-1)+(d+kp)*a*b+(d+km)*(a%2)*(b%2)
            +ip*a-im*(a%2)+op*b-om*(b%2)+j)

PROBES=((0,0),(1,0),(2,0),(3,0),(0,1),(0,2),(0,3),(1,1),(2,2))
def recover(f):
    ip=Q(f[3,0]-f[1,0],2)
    d=Q(f[2,0]-f[0,0],2)-ip
    im=ip-f[1,0]+f[0,0]
    op=Q(f[0,3]-f[0,1],2)
    om=op-f[0,1]+f[0,0]
    def delta(a,b):
        return f[a,b]-f[a,0]-f[0,b]+f[0,0]
    # The eight-probe decoder reconstructs the optional (0,2) value.
    f=dict(f)
    f.setdefault((0,2),f[0,0]+2*d+2*op)
    kp=Q(delta(2,2),4)-d
    km=delta(1,1)-2*d-kp
    return (d,kp,km,ip,im,op,om,f[0,0]-d)

def main():
    count=0
    for i in range(8):
        s=tuple(int(i==j) for j in range(8))
        f={p:response(s,*p) for p in PROBES}
        assert recover(f)==s
        assert f[2,0]-f[3,0]+f[1,0]==f[0,2]-f[0,3]+f[0,1]
        f.pop((0,2))
        assert recover(f)==s
        for a in range(16):
            for b in range(16):
                assert response(recover(f),a,b)==response(s,a,b)
                count+=1
    # The consumer sees j_plus+j_minus; this direction is identically zero.
    def expanded_j(jp,jm,a,b):
        return response((1,2,3,4,5,6,7,jp+jm),a,b)
    for a,b in PROBES:
        assert expanded_j(11,19,a,b)==expanded_j(12,18,a,b)
    print(json.dumps({"status":"PASS","formal_signature_dimension":8,
          "independent_probe_count":8,"consistency_probe_count":1,
          "basis_grid_equalities":count,"invisible_j_difference":True,
          "scope":"formal response formula; no new graph coverage claim"},indent=2))

if __name__=='__main__':
    main()
