#!/usr/bin/env python3
"""Reproducible lexical census, not a census of mathematical meanings."""
import argparse,collections,gzip,json,re,subprocess
from pathlib import Path
HERE=Path(__file__).resolve().parent
parser=argparse.ArgumentParser()
parser.add_argument('--revision',default='2af5e7549da1f7a6cba03e42e15500eb060774c2')
REV=parser.parse_args().revision
patterns={str(n):re.compile(r'(?<![\w])'+str(n)+r'(?![\w])') for n in (1,23,25,63,65,67)}
themes={
 'micro_macro':r'microcosm|macrocosm|micro.?macro',
 'hidden_memory':r'hidden.bit|extra.bit|one.bit|two.bit|sidecar|lost.coordinate',
 'zeckendorf':r'zeckendorf|golden.mean|base.phi|base.φ|zero.slot',
 'dualities':r'rational.{0,12}irrational|algebraic.{0,12}transcend|positive.{0,12}negative|even.{0,12}odd|dualit',
 'square_sum':r'square.sum|square.total|square.sigma',
 'fields_codes':r'golay|mathieu|fano|finite.field|hamming|quadrat|heegner|class.group',
 'curves':r'elliptic|genus|jacobian|selmer|torsion',
 'dynamics':r'collatz|syracuse|logistic|dynamical|cubic.root',
 'spectra':r'hamilton|tournament|metagraph|eigenvalue',
 'units':r'identity.element|multiplicative.identity|field.with.one|unit.object|singleton|empty.product'}
regexes={k:re.compile(v,re.I) for k,v in themes.items()}
raw=subprocess.check_output(['git','ls-tree','-r','-z',REV]).split(b'\0')
by_blob=collections.defaultdict(list)
for row in raw:
 if not row:continue
 meta,path=row.split(b'\t',1); path=path.decode(); blob=meta.split()[2].decode()
 if not path.endswith(('.md','.tex','.org')) or path.startswith(('agents/','.','05-knowledge/results/six_roots_hidden_bits_20261103/')):continue
 by_blob[blob].append(path)
p=subprocess.Popen(['git','cat-file','--batch'],stdin=subprocess.PIPE,stdout=subprocess.PIPE)
records=[];files_by_number=collections.Counter();hits_by_number=collections.Counter();total_bytes=0
for blob,paths in sorted(by_blob.items()):
 p.stdin.write((blob+'\n').encode());p.stdin.flush()
 header=p.stdout.readline().split();size=int(header[2]);data=p.stdout.read(size);p.stdout.read(1);total_bytes+=size
 text=data.decode('utf-8','replace');counts={n:len(r.findall(text)) for n,r in patterns.items()};counts={n:c for n,c in counts.items() if c}
 tags=[k for k,r in regexes.items() if r.search(text)]
 if not counts:continue
 locations={}
 for n in counts:
  found=[i for i,line in enumerate(text.splitlines(),1) if patterns[n].search(line)]
  locations[n]={'line_count':len(found),'first_lines':found[:12]}
 for path in paths:
  records.append({'path':path,'blob':blob,'counts':counts,'locations':locations,'themes':tags})
  files_by_number.update(counts.keys());hits_by_number.update(counts)
p.stdin.close();p.wait()
result={'revision':REV,'scope':'all tracked .md/.tex/.org; excludes dot paths, agents and this packet; exact standalone decimal numerals; lexical appearances include historical/unproved/citation/index mentions, not validated meanings','files_scanned':sum(map(len,by_blob.values())),'unique_blobs_scanned':len(by_blob),'unique_blob_bytes':total_bytes,'files_with_number':dict(files_by_number),'raw_hits':dict(hits_by_number),'records':records}
with gzip.GzipFile(filename=str(HERE/'occurrences.json.gz'),mode='wb',mtime=0) as f:f.write(json.dumps(result,separators=(',',':')).encode())
summary={k:v for k,v in result.items() if k!='records'}
(HERE/'survey_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps(summary,indent=2))
