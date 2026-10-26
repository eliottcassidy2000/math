#!/usr/bin/env python3
"""Portable derivative of the frozen Session67 core-triangle audit.

Default: offline snapshot SHA256 and locally computed Git-blob SHA1 checks.
This DOES NOT verify membership at the manifest's named commit.
Optional --repo adds direct git-show comparison against that commit.
All numerical checks are unchanged from audit_core_triangle.py.
"""
from collections import Counter
from fractions import Fraction as Q
from itertools import combinations, permutations
from math import comb, factorial
from pathlib import Path
import argparse, hashlib, json, subprocess
from verify import parts, exponent, adjacency, hpairs, core_fixed, z

base=Path(__file__).resolve().parent
ap=argparse.ArgumentParser(description=__doc__)
ap.add_argument('--metagraph-dir', type=Path, default=base.parent/'metagraph',
                help='Peer lane containing receipt.json, source_manifest.json and frozen sources')
ap.add_argument('--repo', type=Path,
                help='Optional existing Git checkout containing the declared base commit; no fetch is performed')
args=ap.parse_args()
peer=args.metagraph_dir.resolve()
counts=Counter()
def need(ok,label):
    counts[label]+=1
    if not ok:raise RuntimeError(label)
def E(n):
    return sum((Q(factorial(n)//factorial(n-m),1<<((m*(2*n-m-1)+2)//3))
                for m in range(3,n+1)),Q(0))
def encode(rows,p):
    return sum(((rows[p[a]]>>p[b])&1)<<k for k,(a,b) in enumerate(combinations(range(len(rows)),2)))

r=json.loads((peer/'receipt.json').read_text());N=r['arithmetic_order'];a=r['ordinary_counts'];c=r['pair_free_counts'];T=r['core_size_triangle']
# Offline mode authenticates bytes against the supplied manifest, not commit membership.
manifest=json.loads((peer/'source_manifest.json').read_text())
for item in manifest['files']:
    data=(peer/item['local_path']).read_bytes()
    need(hashlib.sha256(data).hexdigest()==item['sha256'],'source SHA256')
    header=b'blob '+str(len(data)).encode('ascii')+b'\0'
    need(hashlib.sha1(header+data).hexdigest()==item['git_blob'],'source git blob')
    if args.repo is not None:
        original=subprocess.check_output(['git','show',manifest['commit']+':'+item['repository_path']],
                                         cwd=args.repo.resolve())
        need(data==original,'source snapshot equals named Git commit')
need((peer/'receipt.json').read_bytes()==(peer/'receipt-optimized.json').read_bytes(),'peer normal optimized identical')
# Independent integer Burnside pair-orbit loss for every odd type.
for n in range(N+1):
    L=Q(1<<comb(n,2),factorial(n));eps=Q(a[n],1)/L-1
    need(0<=eps<=E(n),'exact epsilon support estimate')
    for lam in parts(n):
        moved=sum(d for d in lam if d>1)
        if moved:
            loss=comb(n,2)-exponent(lam)
            need(loss>=(moved*(2*n-moved-1)+2)//3,'odd type support loss')
    if n>=16:need(E(n)<=Q(2*n**3,1<<(n-1)),'rational geometric majorant n>=16')
# All-order inequalities tested on the whole exact finite triangle.
for n in range(1,N+1):
    for d in range(1,n):
        need(sum(T[n][:n-d+1])<=comb(n-1,d)*a[n-d],'core deficit tail inequality')
    for q in range(1,n+1):
        binomial=comb(n-1,q-1)
        delta=1-Q(T[n][q],a[q]*binomial)
        bound=Q(q*(q-1),1<<(q-1))*(1+E(q-1))+Q(3,2)*E(q)
        need(0<=delta<=bound,'explicit uniform relative error inequality')
    if n>=2:
        need(1-Q(c[n],a[n])<=Q(n*(n-1),1<<(n-1))*(1+E(n-1)),
             'ordinary pairfree mass deficit')
# Strong hostile: fixed q=3 does not satisfy the growing-q conclusion.
for n in range(3,N+1):
    exact=(comb(n-1,2)+2*(n%3==0))//3
    need(T[n][3]==exact,'fixed q3 exact hostile')
# Independent full labelled orbit exhaustion, actual rigid core counts.
rigid=[]
for n in range(1,7):
    allperms=list(permutations(range(n)));seen=set();classes=nonrigid=rigid_core=0
    for code in range(1<<comb(n,2)):
        if code in seen:continue
        rows=adjacency(n,code);orbit={encode(rows,p) for p in allperms};seen.update(orbit)
        classes+=1;aut=factorial(n)//len(orbit);nonrigid+=aut>1
        rigid_core+=aut==1 and not hpairs(rows)
    L=Q(1<<comb(n,2),factorial(n))
    need(classes==a[n],'independent literal class count')
    need(nonrigid<=Q(3,2)*(a[n]-L),'nonrigid class mass bound')
    need(rigid_core>=a[n]-(n-1)*a[n-1]-Q(3,2)*(a[n]-L),'rigid pairfree lower bound')
    for total in range(n,N+1):
        need(rigid_core*comb(total-1,n-1)<=T[total][n],'actual rigid core uniform lower bound')
    rigid.append(dict(n=n,all_classes=classes,nonrigid_classes=nonrigid,rigid_pairfree_classes=rigid_core))
selected=[]
for n in (12,16,20,30,40):
    L=Q(1<<comb(n,2),factorial(n))
    selected.append(dict(n=n,epsilon=str(Q(a[n])/L-1),E=str(E(n)),
                         core_fraction=str(Q(c[n],a[n])),
                         coarse_uniform_error=str(Q(n*(n-1),1<<(n-1))*(1+E(n-1))+Q(3,2)*E(n))))
out=dict(status='PASS',
         provenance_mode='named-commit' if args.repo is not None else 'offline-snapshot',
         declared_source_commit=manifest['commit'],
         named_commit_verified=args.repo is not None,
         snapshot_sha256_verified=True,
         snapshot_git_blob_sha1_verified=True,
         frozen_parent_script_sha256='676b704e028f97a05fb1641a4550564206b9cffa37af783a2b3cfc84be52550d',
         checks=dict(sorted(counts.items())),total_checks=sum(counts.values()),
         rigid_census=rigid,selected_exact_bounds=selected,
         peer_report_sha256=hashlib.sha256((peer/'REPORT.md').read_bytes()).hexdigest(),
         peer_script_sha256=hashlib.sha256((peer/'verify_core_families.py').read_bytes()).hexdigest(),
         peer_receipt_sha256=hashlib.sha256((peer/'receipt.json').read_bytes()).hexdigest())
print(json.dumps(out,indent=2,sort_keys=True))
