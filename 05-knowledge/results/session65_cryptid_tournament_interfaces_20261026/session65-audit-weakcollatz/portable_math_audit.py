#!/usr/bin/env python3
"""NEW PORTABLE CHECKER, not a replacement for the sealed independent_audit.py.

Preserves every authored-math, producer-output, native and receipt-integrity
check. Independently recompiles/runs only our native_mb.c. External source blob
hashes, Lean token scan and import closure are NOT rerun: those checks remain
explicitly historical observed receipts. See PORTABILITY.md for the distinction.
Requires Python standard library and gcc, no network or downloaded source code.
"""
import ast
from fractions import Fraction
from hashlib import sha256, sha1
from itertools import product
import json
from pathlib import Path
import subprocess
import tempfile
import shutil

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

compiler=shutil.which('gcc')
gate(compiler is not None,'portable native replay requires gcc; no install is attempted')
with tempfile.TemporaryDirectory(prefix='s65-authored-native-') as tmp:
    binary=Path(tmp)/'native_mb'
    build=subprocess.run([compiler,'-std=c11','-O2','-Wall','-Wextra','-Werror',
                          str(HERE/'native_mb.c'),'-o',str(binary)],
                         capture_output=True,text=True,timeout=60)
    gate(build.returncode==0,'authored native compile failed: '+build.stderr)
    run=subprocess.run([str(binary)],capture_output=True,text=True,timeout=60)
    gate(run.returncode==0,'authored native execution failed: '+run.stderr)
    gate(run.stdout.encode()==(HERE/'native_mb.json').read_bytes(),'native replay differs from archived native JSON')
    native=json.loads(run.stdout)
gate(native['steps']==sum(costs) and native['ones']==4098,'independent native vs macro')
gate([row['cost'] for row in native['macro_checkpoints']]==costs,'native macro checkpoints')
gate(native['fresh_native_macro_cases']==121 and native['status']=='PASS','fresh native cases')

# These are historical source-audit receipts, not newly repeated source checks.
# The old checker remains unchanged, but is not a portable execution target.
historical_manifest=json.loads((HERE/'manifest.json').read_text())
for path,digest in historical_manifest['audit_sha256'].items():
    gate(sha256((HERE/path).read_bytes()).hexdigest()==digest,'historical audit seal '+path)
for path,digest in historical_manifest['target_sha256'].items():
    gate(sha256((TARGET/path).read_bytes()).hexdigest()==digest,'historical producer seal '+path)
historical=json.loads((HERE/'normal.json').read_text())
gate((HERE/'normal.json').read_bytes()==(HERE/'optimized.json').read_bytes(),'historical source audit normal/-O bytes')
for key,value in {'orbit_with_seed':orbit,'native_costs':costs,'native_total':sum(costs),
                  'independent_native':native,'sigma_ratios':ratios,'cylinders':censuses,
                  'orientation_pairs':pairs,'negative_periodic_extensions':periodic}.items():
    gate(historical[key]==value,'portable math differs from historical audit '+key)
remote=json.loads((HERE/'remote_git_pin.json').read_text())
gate(remote['commit']=='3c43f97768d4ea06d70e0eb154ba88060bc36dec' and not remote['truncated'],
     'historical remote tree pin')
gate(len(remote['checks'])==33 and all(row['match'] and row['expected']==row['remote'] for row in remote['checks']),
     'historical remote source-pin consistency')
print(json.dumps({
    'status':'PASS_PORTABLE_AUTHORED_MATH_AND_NATIVE_REPLAY',
    'orbit_with_seed':orbit,'native_costs':costs,'native_total':sum(costs),
    'independent_native':native,'native_execution_performed':True,'sigma_ratios':ratios,
    'cylinders':censuses,'orientation_pairs':pairs,'negative_periodic_extensions':periodic,
    'historical_source_audit':{
        'status':'HISTORICAL_RECEIPTS_HASH_VERIFIED_NOT_SOURCE_CHECKS_REPERFORMED',
        'source_blobs_previously_verified':historical['source_blobs_verified'],
        'targeted_lean_import_closure_previously_inspected':historical['targeted_lean_import_closure'],
        'unavailable_root_imports_at_original_audit':historical['unavailable_root_imports'],
        'remote_manifest_entries_previously_verified':len(remote['checks'])},
    'external_tinycryptids_sources_available_or_required':False,
    'external_tinycryptids_software_executed':False,
    'lean_or_rust_build_performed':False,
    'scope':'Portable authored mathematics/native replay only. Original source-integrity and Lean-token/import audits remain historical observed receipts, not offline-reperformed checks. No all-seed halting conclusion.'
},indent=2,sort_keys=True))
