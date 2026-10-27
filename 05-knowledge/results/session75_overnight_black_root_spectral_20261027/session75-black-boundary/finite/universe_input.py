#!/usr/bin/env python3
import csv,gzip,pathlib
H=pathlib.Path(__file__).resolve().parent.parent/'inputs'
for n in range(4,10):
 p=H/f'classes_n{n}.tsv'if n<9 else H/'n9_classes.tsv.gz'
 with (gzip.open if n==9 else open)(p,'rt')as f:
  for r in csv.DictReader(f,delimiter='\t'):print(n,r['canon'if n<9 else'code'],r['Aut'])
