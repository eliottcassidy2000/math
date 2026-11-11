"""Seal independent nine-vertex results and extend only certified registry rows."""
from pathlib import Path
import csv,gzip,json
HERE=Path(__file__).resolve().parent
def read(name):return json.loads((HERE/name).read_text(encoding='utf-8-sig'))
def need(q,s):
    if not q:raise RuntimeError(s)
s=read('screen9.json');a=read('audit9.json');reg=read('registry.json')
for k in ('unit_sources','unit_selections_sum','unit_orbits_sum','unit_odd_coefficients'):
    need(s[k]==a[k],'independent nine-vertex statistic '+k)
cert=[{k:int(v) for k,v in r.items()} for r in csv.DictReader((HERE/'seed_witnesses9.tsv').open(),delimiter='\t')]
need(sorted(r['c'] for r in cert)==s['unit_coefficients'],'coefficient support independent equality')
old=json.loads((HERE.parent/'odd_c_seed_symmetry_20261110/verification.json').read_text())['new_seeds']
need(sorted(r['c'] for r in old)==a['missing_through_4097'],'same 33 lower-bound and upper-bound coefficients')
need(all(len(r['rows'])==10 for r in old),'eleven-vertex upper roots')
need(a['ordinary_sources_checked']==191536,'complete ordinary nine universe')
need(s['extension_support_distinct']==s['merged'] and not s['support_collisions'],'nine converse separation')
for k,v in s.items():
    if k in reg and isinstance(v,int):reg[k]['terms']=reg[k]['terms'][:8]+[v]
for prefix in ('extension_support','extension_distribution'):
    # Exact support separation leaves only the known converse pairs.
    reg[prefix+'_ambiguous_pairs']['terms']=reg[prefix+'_ambiguous_pairs']['terms'][:8]+[s['ordinary']-s['merged']]
    reg[prefix+'_largest_fibre']['terms']=reg[prefix+'_largest_fibre']['terms'][:8]+[2]
    reg[prefix+'_ambiguous_pairs']['order9_method']='derived from exact support separation and converse pairing'
    reg[prefix+'_largest_fibre']['order9_method']='derived from exact support separation and existence of a non-self-converse class'
for v in reg.values():v['last_order']=len(v['terms'])
(HERE/'registry.json').write_text(json.dumps(reg,indent=2,sort_keys=True)+'\n',encoding='utf-8')
with (HERE/'sequences.tsv').open('w',newline='',encoding='utf-8') as f:
    w=csv.writer(f,delimiter='\t');w.writerow(['id','kind','universe','definition']+[f'n={n}' for n in range(1,10)])
    for k,v in reg.items():w.writerow([k,v['kind'],v['universe'],v['definition']]+v['terms']+['']*(9-len(v['terms'])))
# Import the small independent functions without executing an experiment.
scope={'__file__':str(HERE/'verify.py')};src=(HERE/'verify.py').read_text();exec(src[:src.index('rs=[json.loads')],scope)
H=scope['H'];tail=scope['tail'];child=scope['child']
def decode9(code):
    out=[0]*9;b=35
    for i in range(9):
        for j in range(i+1,9):
            if code>>b&1:out[i]|=1<<j
            else:out[j]|=1<<i
            b-=1
    return out
tailchecks=0
for r in cert[::61]+cert[-1:]:
    rows=decode9(r['canonical_code'])
    if r['converse']:rows=[511^v^(1<<i) for i,v in enumerate(rows)]
    need(H(rows)==r['H'],'nine witness H')
    for t in (1,2):
        need(H(tail(rows,r['selected_mask'],t))==1+r['c']*2**(r['first_exponent']+t-1),'nine literal tail')
        tailchecks+=1
records=[json.loads(x) for x in gzip.open(HERE/'records.jsonl.gz','rt')]
need(all(min(r['extension'])==r['H'] for r in records),'support recovers parent H')
source_obstructions=read('source_obstructions.json');obstructionchecks=0
for n,ws in source_obstructions.items():
    if not ws:continue
    w=ws[0];n=int(n)
    for rep,mask,target in zip(w['sources'],w['attachments'],w['pair']):
        r=next(r for r in records if r['n']==n and r['rep']==rep)
        need(H(child(r['rows'],mask))==target,'incompatible pair actual witness');obstructionchecks+=1
    need(not any(set(w['pair'])<=set(r['extension']) for r in records if r['n']==n and r['H']==w['H']),'no common source')
# Check the bivariate tc11 law against every coefficient visible in the complete n<=8 corpus.
states={(0,0):1}
for n in range(1,9):
    for k in range(4):states[n,k]=states.get((n-1,k),0)+states.get((n-3,k-1),0)
def coeff(n,k):
    ans=0
    for order,power,mult in [(5,0,1),(6,1,2),(7,2,1)]:
        for l in range(n-order+1):
            for j in range(k-power+1):ans+=mult*states.get((l,j),0)*states.get((n-order-l,k-power-j),0)
    return ans
for n in range(1,9):
    for k in range(4):need(sum(r['n']==n and r['tc']==11 and r['aut']==3**k for r in records)==coeff(n,k),'bivariate q11 grammar')
receipt={'status':'PASS','independent_full_order9_seed_scans':2,'ordinary_order9_classes':191536,
         'certified_odd_coefficients':len(cert),'minimum_root_order_11_coefficients':a['missing_through_4097'],
         'extra_literal_nine_base_tails':tailchecks,'incompatible_pair_witness_checks':obstructionchecks,
         'bivariate_tc11_checks':32,'registry_sequences':len(reg),'registry_rows_extended_to_nine':sum(len(v['terms'])==9 for v in reg.values())}
(HERE/'final_checks.json').write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf-8')
print(json.dumps(receipt))
