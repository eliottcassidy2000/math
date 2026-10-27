#!/usr/bin/env python3
import csv,gzip,pathlib,hashlib,json
H=pathlib.Path(__file__).resolve().parent
for item in json.loads((H/'SOURCE_TABLE_BINDINGS.json').read_text())['source_tables_exact_git_bytes']:
 p=H/item['local']
 if hashlib.sha256(p.read_bytes()).hexdigest()!=item['sha256']:raise RuntimeError('source hash mismatch')
for n in range(4,10):
 p=H/'inputs'/(f'classes_n{n}.tsv' if n<9 else 'n9_classes.tsv.gz')
 with (gzip.open if n==9 else open)(p,'rt') as f:
  for row in csv.DictReader(f,delimiter='\t'):print(n,row['canon' if n<9 else 'code'],row['Aut'])
