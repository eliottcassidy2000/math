"""Signed absolute and relative error for every order in the inherited exact cache."""
import csv,json,sys,hashlib
from pathlib import Path
from decimal import Decimal,localcontext
from certified_interval import lower_and_tail,sector,odd_nonunit_partitions
from fractions import Fraction
sys.set_int_max_str_digits(0)
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[3]
CACHE=ROOT/'05-knowledge/results/session270_chanpak_enumeration_20261111/enumeration/prefix.json'
if not CACHE.exists():
    raise SystemExit('Required inherited exact prefix missing: '+str(CACHE))
values=[int(x) for x in json.loads(CACHE.read_text())['values']]

def log10q(q):
    if q==0:return '-Infinity'
    with localcontext() as ctx:
        ctx.prec=32
        return str(Decimal(q.numerator).log10()-Decimal(q.denominator).log10())

rows=[]
for n,a in enumerate(values):
    for s in (0,3,6,9):
        L=sum((sector(n,mu) for m in range(min(n,s)+1) for mu in odd_nonunit_partitions(m)),Fraction())
        error=Fraction(a)-L
        assert error>=0
        rows.append(dict(n=n,support_cut=s,signed_error_direction='exact' if not error else 'underestimate',
          log10_absolute_error=log10q(error),log10_relative_error=log10q(error/a),
          exact_error_numerator=str(error.numerator) if n<=24 else 'regenerate: cache[n] - support_sum(n,s)',
          exact_error_denominator=str(error.denominator),
          exact_error_sha256=hashlib.sha256((str(error.numerator)+'/'+str(error.denominator)).encode()).hexdigest()))
with (HERE/'ERRORS.csv').open('w') as f:
    writer=csv.DictWriter(f,fieldnames=rows[0]);writer.writeheader();writer.writerows(rows)
checks=[]
for n in (0,1,2,3,6,14,20,32,50,64,100,150,200,250):
    for s in (0,3,6,9):
        L,B=lower_and_tail(n,s)
        E=Fraction(values[n])-L
        assert 0<=E<=B
        checks.append({'n':n,'cut':s,'tail_bound_log10':log10q(B),'actual_error_log10':log10q(E)})
(HERE/'CONTROLS.json').write_text(json.dumps({'cache_sha256':hashlib.sha256(CACHE.read_bytes()).hexdigest(),
  'signed_error_rows':len(rows),'display_precision_digits':32,'display_is_not_certificate':True,
  'exact_rational_interval_checks':checks},indent=2)+'\n')
print('PASS:',len(rows),'signed errors;',len(checks),'exact rational intervals')
