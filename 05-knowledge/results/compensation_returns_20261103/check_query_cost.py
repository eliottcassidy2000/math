"""Literal controls for the named query that reads the lone-zero charge."""
from pathlib import Path
import hashlib
import json
import subprocess
import tempfile

HERE = Path(__file__).resolve().parent


def need(ok, why):
    if not ok:
        raise ValueError(why)


def birth(rows, out):
    return [r | (0 if out >> i & 1 else 1 << len(rows)) for i, r in enumerate(rows)] + [out]


def payload(graphs):
    return '\n'.join(' '.join(map(str, [i, len(r)] + r)) for i, r in enumerate(graphs)) + '\n'


def main():
    source = HERE / 'recovery/one_apex.json'
    counter = HERE / 'audit/exact_counter.cpp'
    children = [x for x in json.loads(source.read_text())['actual_source_controls'] if x['kind'] == 'child']
    records, extended = [], []
    for c in children:
        rows, names = c['graph']['rows'], c['graph']['names']
        w, n = names.index(0), len(rows)
        full = (1 << n) - 1
        rev = [full ^ r ^ (1 << i) for i, r in enumerate(rows)]
        query = birth(birth(rev, full ^ (1 << w)), full)
        preparation = birth(rows, 1 << w)
        extended += [query, preparation]
        records.append(dict(source_H=c['H'], original_names=names, mark=0, root=99,
                            source_rows=rows, query_names=names+[10000,10001], query_rows=query,
                            prepared_names=names+[10002], prepared_rows=preparation))
    with tempfile.TemporaryDirectory(prefix='named_query_charge_') as tmp:
        exe = Path(tmp) / 'counter.exe'
        subprocess.run(['g++','-O3','-std=c++20',str(counter),'-o',str(exe)],check=True)
        counts = [json.loads(s)['H'] for s in subprocess.check_output(
            [str(exe),'h'],input=payload(extended),text=True).splitlines()]
        profiles = [json.loads(s) for s in subprocess.check_output(
            [str(exe)],input=payload([r['source_rows'] for r in records]),text=True).splitlines()]
    for i, (r, p) in enumerate(zip(records, profiles)):
        w = r['original_names'].index(0)
        h, e, s2 = p['H'], p['E1'][w], p['S2'][w]
        query_h, prepared_h = counts[2*i:2*i+2]
        need(h == r['source_H'] and e == 1 and p['deletion_H'][w] == 3, 'source guards')
        need(p['E2'][r['original_names'].index(99)] == h, 'root99 end zero')
        need(query_h == h + 2*s2, 'converse singleton-query type translation')
        need(prepared_h == h + s2 - e, 'actual one-birth charge')
        r.update(query_H=query_h, endpoint=e, start_cover=s2, preparation_charge=s2-e, prepared_H=prepared_h)
    need([r['query_H'] for r in records] == [8675,8501], 'literal queried costs')
    result = dict(status='FINITE-EXACT controls for proved cost observer', records=records,
                  source_sha256={str(p.relative_to(HERE)):hashlib.sha256(p.read_bytes().replace(b'\r\n',b'\n')).hexdigest()
                                 for p in [source,counter,Path(__file__)]})
    (HERE/'query_cost.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8',newline='\n')
    print('PASS: named queries8675/8501 read exact preparation costs3632/3536')


if __name__ == '__main__':
    main()
