#!/usr/bin/env python3
"""Reject-path controls for independent exact verifier; survives Python -O."""
import json
from audit_checks import *
from audit_families import build

def reject(name,fun):
    try:fun()
    except RuntimeError as e:return {'control':name,'rejection':str(e)}
    raise RuntimeError('false acceptance: '+name)

def main():
    t,_=build(1,False);p=(2,1,0,5,6,7,8,9,4,3);q=(9,8,7,6,5,0,1,3,4,2)
    tests=[reject('nonpermutation',lambda:reverse(t,p[:-1]+(p[0],))),reject('reversed nondirected word',lambda:reverse(t,p[::-1])),reject('blue C3 move',lambda:replay(from_rows([2,4,1]),[(0,1,2)])),reject('nonstrong target is not a strong core',lambda:replay(t,[p,q],require_core=True)),reject('malformed tournament encoding',lambda:from_rows([0,0,0]))]
    print(json.dumps({'status':'PASS','reject_controls':tests},sort_keys=True,indent=2))
if __name__=='__main__':main()
