"""Decode the inherited complete order-nine table into an explicit adjacency input."""
from pathlib import Path
import csv,gzip,hashlib,json,sys
HERE=Path(__file__).resolve().parent
p=HERE.parent/'blue_quotient_20261020b/n9_classes.tsv.gz'
count=0
with gzip.open(p,'rt') as f,Path(sys.argv[1]).open('w') as out:
    for r in csv.DictReader(f,delimiter='\t'):
        code=int(r['code']);dual=int(r['converse_code'])
        if code>dual:continue
        rows=[0]*9;b=35
        for i in range(9):
            for j in range(i+1,9):
                if code>>b&1:rows[i]|=1<<j
                else:rows[j]|=1<<i
                b-=1
        out.write(' '.join(map(str,[code,dual,int(r['Aut']),int(r['SC_H'])]+rows))+'\n');count+=1
(HERE/'input9.json').write_text(json.dumps({'source':str(p.relative_to(HERE.parents[2])).replace('\\','/'),'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'merged_rows':count,'decoder':'upper triangle lexicographic, most significant bit first; bit one means i -> j'},indent=2)+'\n')
print(count,'complete converse-merged input rows')
