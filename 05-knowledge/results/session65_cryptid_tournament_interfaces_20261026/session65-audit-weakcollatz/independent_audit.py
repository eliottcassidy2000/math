#!/usr/bin/env python3
"""Independent exact audit. No asserts, imports of target code, or external programs.

The native simulator is compiled and run separately. This script checks its JSON
as data. All source files are read-only; its only output is stdout.
"""
import ast
from fractions import Fraction
from hashlib import sha256, sha1
from itertools import product
import json
from pathlib import Path
import re

HERE=Path(__file__).resolve().parent
TARGET=HERE.parent/'session65-weakcollatz'
def gate(ok, why):
    if not ok: raise RuntimeError(why)

def advance(n):
    r=n%3
    return None if r==2 else (5*n+18+4*r)//3

def itinerary(n,depth):
    digits=[]
    for _ in range(depth):
        r=n%3
        if r==2: return None
        digits.append(r)
        n=(5*n+18+4*r)//3
    return tuple(digits)

def offset(w):
    b=0
    for j,r in enumerate(w): b=5*b+(18+4*r)*3**j
    return b

def valuation_unit(n):
    gate(n!=0,'unit requested at zero')
    v=0
    while n%3==0: n//=3;v+=1
    return v,n%3

def divisor_sum(n):
    # Prime-power product, independent of producer's divisor enumeration.
    remaining=n;answer=1;p=2
    while p*p<=remaining:
        power=1;subtotal=1
        while remaining%p==0:
            remaining//=p;power*=p;subtotal+=power
        answer*=subtotal;p+=1
    if remaining>1:answer*=remaining+1
    return answer

sequence=[6,16,34,64,114,196,334,564,946,1584,2646,4416,7366,12284]
orbit=[0];costs=[]
for n in [0]+sequence:
    gate(n==orbit[-1],'exact orbit mismatch')
    q,r=divmod(n,3)
    costs.append(5*q*q+(19+6*r)*q+15+12*r if r<2 else 6*q+12)
    nxt=advance(n)
    if nxt is not None:orbit.append(nxt)
gate(orbit==[0]+sequence and nxt is None,'orbit end')
gate(sum(costs)==47176870,'native-cost arithmetic')
ratios={str(n):str(Fraction(divisor_sum(n),n)) for n in sequence}
gate([n for n in sequence if divisor_sum(n)%n==0]==[6],'multiply-perfect hostile')
gate([n for n in sequence[:-1] if divisor_sum(n)==advance(n)]==[946],'sigma successor coincidence')
gate(advance(-9)==-9 and advance(-11)==-11,'negative fixed-point controls')

censuses=[];pairs=0;periodic=0
for m in range(1,11):
    mod=3**m
    classes={}
    for w in product(range(2),repeat=m):
        b=offset(w)
        representative=(-b*pow(5**m,-1,mod))%mod
        gate(representative not in classes,'inverse collision')
        classes[representative]=w
        gate(itinerary(representative,m)==w,'forward affine address')
        # The periodic extension is a negative rational fixed by g^m.
        point=Fraction(-b,5**m-3**m)
        gate(-11<=point<=-9,'negative periodic extension range')
        x=point
        for digit in w:
            r=x.numerator*pow(x.denominator,-1,3)%3
            gate(r==digit,'rational periodic itinerary')
            x=(5*x+18+4*r)/3
        gate(x==point,'rational periodic closure');periodic+=1
    direct={n:itinerary(n,m) for n in range(mod) if itinerary(n,m) is not None}
    gate(classes==direct and len(classes)==2**m,'direct/inverse census')
    for n,w in classes.items():
        x,y=n,n+mod
        for _ in range(m):x,y=advance(x),advance(y)
        gate(y-x==5**m,'guarded affine difference')
    if m<=9:
        ordered=sorted(classes,key=lambda n:tuple(d^(j%2) for j,d in enumerate(classes[n])))
        for a_index,a in enumerate(ordered):
            for b in ordered[a_index+1:]:
                w,v=classes[a],classes[b]
                j=next(j for j in range(m) if w[j]!=v[j])
                valuation,unit=valuation_unit(b-a)
                gate(valuation==j,'first differing digit valuation')
                gate(unit==((-1)**j*(v[j]-w[j]))%3==1,'alternating lex orientation')
                pairs+=1
    censuses.append({'m':m,'modulus':mod,'safe_classes':len(classes)})

producer=json.loads((TARGET/'bridge_results.json').read_text())
gate(producer['orbit']==sequence,'producer orbit')
gate(producer['native_costs']==costs,'producer native cost list')
gate(producer['sigma_ratios']==ratios,'producer sigma ratios')
gate(producer['alternating_lex_orientation_pairs']==sum((2**m)*(2**m-1)//2 for m in range(1,9)),'producer pair total')
gate((TARGET/'bridge_results.json').read_bytes()==(TARGET/'bridge_results_optimized.json').read_bytes(),'producer normal/-O bytes')
gate(not any(isinstance(n,ast.Assert) for n in ast.walk(ast.parse((TARGET/'verify_bridge.py').read_text()))),'producer contains optimization-sensitive assert')
receipt=json.loads((TARGET/'verification_receipt.json').read_text())
for path,digest in receipt['sha256'].items():
    gate(sha256((TARGET/path).read_bytes()).hexdigest()==digest,'producer sealed hash '+path)

native=json.loads((HERE/'native_mb.json').read_text())
gate(native['steps']==sum(costs) and native['ones']==4098,'independent native vs macro')
gate([row['cost'] for row in native['macro_checkpoints']]==costs,'native macro checkpoints')
gate(native['fresh_native_macro_cases']==121 and native['status']=='PASS','fresh native cases')

manifest=json.loads((TARGET/'tinycryptids/source-manifest.json').read_text())
verified_blobs=[];not_downloaded=[]
for row in manifest['files']:
    path=TARGET/'tinycryptids'/row['path']
    if not path.exists():not_downloaded.append(row['path']);continue
    data=path.read_bytes()
    digest=sha1(b'blob '+str(len(data)).encode()+b'\0'+data).hexdigest()
    gate(digest==row['git_blob_sha'],'external source blob integrity '+row['path'])
    verified_blobs.append(row['path'])

leanroot=TARGET/'tinycryptids/lean'
roots=['TinyCryptids.WeakCollatz.Main','TinyCryptids.WeakCollatz3.Main']
closure=set();todo=list(roots)
while todo:
    module=todo.pop()
    if module in closure:continue
    file=leanroot/(module.replace('.','/')+'.lean')
    gate(file.exists(),'missing targeted local Lean import '+module)
    closure.add(module)
    src=file.read_text()
    gate(not re.search(r'\b(sorry|admit|axiom|unsafe|native_decide)\b',src),'forbidden Lean source token '+module)
    for imp in re.findall(r'^import\s+(\S+)',src,re.M):
        if imp.startswith('TinyCryptids'):todo.append(imp)
gate(len(closure)==15,'targeted import closure size')
root_imports=re.findall(r'^import\s+(\S+)',(leanroot/'TinyCryptids.lean').read_text(),re.M)
gate(all(r in root_imports for r in roots),'root source imports targeted theorems')
missing_root_imports=[imp for imp in root_imports if not (leanroot/(imp.replace('.','/')+'.lean')).exists()]

print(json.dumps({
    'status':'PASS',
    'orbit_with_seed':orbit,'native_costs':costs,'native_total':sum(costs),
    'independent_native':native,'sigma_ratios':ratios,
    'cylinders':censuses,'orientation_pairs':pairs,'negative_periodic_extensions':periodic,
    'source_blobs_verified':len(verified_blobs),'manifest_docs_not_downloaded':not_downloaded,
    'targeted_lean_import_closure':len(closure),'unavailable_root_imports':missing_root_imports,
    'external_software_executed_by_auditor':False,
    'lean_or_rust_build_performed_by_auditor':False,
    'scope':'Independent native and finite arithmetic verification; proofs of infinite statements are in AUDIT.md. No all-seed halting conclusion.'
},indent=2,sort_keys=True))
