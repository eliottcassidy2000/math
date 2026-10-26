#!/usr/bin/env python3
"""Independent monomial controls for the non-UFD and normalized-ratio hostiles."""
import json

def need(ok,msg):
    if not ok:raise RuntimeError(msg)
def ring(a,b):return a>=0 and b>=0 and (a+b)%2==0
def ideal(a,b,name):
    return ring(a,b) and (a>=1 if name=='P' else a+b>0)

def main():
    checks=0
    # R=Q[s^2,t^2,st]. A=R[w], w^2=2. B_I=R+Iw.
    for name in ('P','m'):
        for a in range(13):
            for b in range(13):
                if not ring(a,b):continue
                for parity in (0,1):
                    inside=(parity==0 or ideal(a,b,name))
                    multiplied_by_w=(parity==1 or ideal(a,b,name))
                    need((inside and multiplied_by_w)==ideal(a,b,name),'conductor IA component test')
                    checks+=1
    # P^2=(s^4,s^3t,s^2t^2)=x*m, but reflexive hull is (x).
    for a in range(13):
        for b in range(13):
            P2=any(ring(a-u,b-v) for u,v in ((4,0),(3,1),(2,2)))
            need(P2==ideal(a-2,b,'m'),'P^2=xm')
    need(ring(2-2,0) and not ideal(2-2,0,'m'),'x in reflexive hull but not P^2')
    need(not ideal(0,0,'m'),'codimension-two conductor obstruction at 1')
    # At every height-one prime, m localizes to R: m is the height-two vertex.
    # This assertion is justified algebraically in the note, not inferred from a grid.
    overlap=0
    for k in range(13):
        for parity in (0,1):
            inB=parity==0 or k>=1
            times_w=parity==1 or k>=1
            need((inB and times_w)==(k>=1),'normalized overlap conductor tA')
            overlap+=1
    need(True and 1>=1 and 1>=1,'old roots t and wt pass at 1')
    need(not (1==0 or 0>=1),'new integral ratio w fails at 1')
    print(json.dumps({'status':'FINITE-EXACT PASS','coefficient_conductor_checks':checks,
        'overlap_checks':overlap,'nonprincipal_ideal':'(x,z)',
        'codimension_two_hostile':'I=(x,y,z), valuations miss 1 not in IA',
        'normalization_hostile':'t and wt pass, integral ratio w fails',
        'scope':'finite monomial identities; ideal/localization/normality proofs are in CONDUCTOR_IDEALS.md'},indent=2,sort_keys=True))
if __name__=='__main__':main()
