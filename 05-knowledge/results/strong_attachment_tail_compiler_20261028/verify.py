"""Exact finite audits and pinned-receipt validation; python -O keeps checks."""
from collections import Counter
from fractions import Fraction
from hashlib import sha256
from pathlib import Path
from math import factorial
import argparse,gzip,json,sys
import counting,native
from seal import digest

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
CHECKS=Counter()
def check(ok,label):
    CHECKS[label]+=1
    if not ok:raise RuntimeError(label)
def read(name):return json.loads((HERE/name).read_text(encoding='utf-8'))
def series(path):
    result={}
    for line in path.read_text(encoding='utf-8').splitlines():
        line=line.split('#')[0].strip()
        if line:
            n,v=line.split();result[int(n)]=int(v)
    return result
def strong(a):
    s=[0]*len(a)
    for n in range(1,len(a)):s[n]=a[n]-sum(s[k]*a[n-k] for k in range(1,n))
    s[0]=1 # A051337's explicitly prepended empty-object convention.
    return s

def main():
    if hasattr(sys,'set_int_max_str_digits'):sys.set_int_max_str_digits(0)
    p=argparse.ArgumentParser();p.add_argument('--record',action='store_true')
    p.add_argument('--replay-large',action='store_true');p.add_argument('--replay-roots',action='store_true');args=p.parse_args()
    saved=read('native-normal.json')
    check(saved==read('native-optimized.json'),'normal and optimized native receipts')
    banks,closed=native.closure(7)
    live=dict(bounded_unnamed=closed,actual_atoms=native.atom_experiments(banks),scope_max_order=7)
    # JSON roundtrip normalizes tuples without weakening values.
    check(json.loads(json.dumps(live))==saved,'full native replay')
    residual=closed[5]['order_eight_residual_probes']
    check(len(residual)==5 and all(w['separating_context'] is not None for w in residual),'full reconstruction through six with horizon eight')
    p7=tuple(sum(1<<((i+d)%7) for d in (1,2,4)) for i in range(7))
    for ports,expected in [((0,1,3),4281),((0,1,5),4315)]:
        t=native.marker_attachment(p7,ports,(0,1,2))
        check(native.core.is_strong(t),'marker result is strong')
        check(sum(1 for _ in native.core.paths(t))==expected,'independent marker DFS Hamilton count')
    a=read('ordinary-300.json')['values'];base=read('baseline-300.json')['values']
    check(len(a)==301 and a==base,'complete baseline agreement through 300')
    check(a[:251]==read('ordinary-250.json')['values'],'bound extension preserves prefix')
    old=json.loads((ROOT/'05-knowledge/results/tournament_future_boundary_20261020/receipt.json').read_text())
    check(a[:201]==list(map(int,old['large']['coefficients'])),'independent Fraction profile receipt through 200')
    ref=series(ROOT/'05-knowledge/results/tournament_future_boundary_20261020/a000568-reference.txt')
    check(len(ref)==77 and all(a[n]==v for n,v in ref.items()),'OEIS A000568 table through 76')
    direct=native.source('direct_june','04-computation/a000568_exact_burnside.py')
    for n in range(41):check(a[n]==direct.a000568(n),'independent direct odd partitions through 40')
    core50=read('core-50.json')['values']
    check(core50==list(map(int,old['core_classes'][:51])),'pair-free stored Fraction control through 50')
    fractions=native.source('fraction_october','05-knowledge/results/tournament_future_boundary_20261020/verify.py')
    for core in (False,True):
        truth=fractions.direct_series(20,core=core)
        for bound in range(21):
            vals,_=counting.contracted(bound,core=core)
            check(vals==truth[:bound+1],'all small horizons versus direct polynomial substitution')
        vals,_=counting.contracted(50 if core else 100,core=core)
        check(vals==(core50 if core else a[:101]),'production replay ordinary/core')
    for bound in range(21):
        vals,_=counting.contracted13(bound)
        check(vals==a[:bound+1],'two-stage negative-control boundary audit')
    vals,_=counting.contracted13(100)
    check(vals==a[:101],'two-stage contraction exact through 100')
    s=strong(a);sref=series(HERE/'a051337-reference.txt')
    check(len(sref)==51 and all(s[n]==v for n,v in sref.items()),'OEIS A051337 table through 50')
    for n in range(1,len(a)):
        check(sum(s[k]*a[n-k] for k in range(1,n+1))==a[n],'strong ordered-component inversion')
    derived=dict(A051337=s,pair_free=core50,conventions=dict(A051337_0=1,pair_free_0=1))
    if args.record:(HERE/'derived-series.json').write_text(json.dumps(derived,indent=2)+'\n',encoding='utf-8')
    else:check(derived==read('derived-series.json'),'derived series receipt')
    previous=None
    for n in range(3,51):
        term=Fraction(1<<((n-3)*(n-2)//2+1),3*factorial(n-3))
        if n>=4:check(term>1 and term>previous,'absolute tail historical correction')
        previous=term
    if args.replay_large:
        vals,_=counting.contracted(300);check(vals==a,'large contracted replay')
    roots=read('root-switch-skip.json')
    oldroots=json.loads((ROOT/'05-knowledge/results/relation_tower_frobenius_ports_20261028/root_switch.json').read_text())
    for new,oldr in zip(roots['results'],oldroots['results']):
        for key in ('order','pointed_states','classes','trace','automorphism_order_histogram'):
            check(new[key]==oldr[key],'root-switch full-row receipt agreement')
        check(new['computed_rows']==new['classes'] and new['parent_classification_calls']==new['order']*new['classes'],
              'root-switch exact work reduction')
    bank=json.loads(gzip.decompress((ROOT/roots['order_nine_bank']['source_bank']).read_bytes()))
    encoded=json.dumps(bank,separators=(',',':'))+'\n'
    check(sha256(encoded.encode()).hexdigest()==roots['order_nine_bank']['canonical_json_sha256'],'complete ordered nine-bank digest')
    from root_switch_skip import generate
    parents=[()]
    for n in range(1,10 if args.replay_roots else 8):
        before=parents;parents,rr=generate(before,n)
        check(len(parents)==a[n],'constructive root-switch replay count')
        if n<=7:
            control,_=generate(before,n,skip_blocks=False)
            check(parents==control,'matched whole-row constructive control')
    if args.replay_roots:check([list(t) for t in parents]==bank['rows'],'complete order-nine regeneration')
    import incoming
    incoming_live=incoming.run(100);incoming_saved=read('incoming.json')
    check(incoming_live['selector_class_bijections']==incoming_saved['selector_class_bijections'],
          'incoming selector class bijections through seven')
    check(incoming_saved['degree_bucket_bound']==300 and incoming_saved['degree_bucket_complete_agreement'],
          'incoming degree-bucket large receipt')
    pins=HERE/'sources.json'
    check(pins.exists(),'dependency pins exist')
    for path,expected_hash in json.loads(pins.read_text()).items():
        check(digest(ROOT/path)==expected_hash,'pinned dependency')
    if not args.record:
        check((HERE/'manifest.json').exists(),'packet manifest exists')
        for path,expected_hash in read('manifest.json').items():
            check(digest(HERE/path)==expected_hash,'sealed packet file')
    report=dict(checks=dict(sorted(CHECKS.items())),scope=dict(ordinary=300,external_ordinary=76,
        direct_partitions=40,strong=300,external_strong=50,pair_free=50,native_max_full=7,
        reconstructed_all_presentations_through=6,reconstruction_horizon=8,constructive_bank=9),
        a300_digits=len(str(a[300])),a300_sha256=sha256(str(a[300]).encode()).hexdigest(),
        normal_optimized_equal=True)
    if args.record:(HERE/'validation.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(report,sort_keys=True))

if __name__=='__main__':main()
