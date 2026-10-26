import csv,gzip,json,pathlib,hashlib
p=pathlib.Path(__file__).resolve().parent
census=json.loads((p/'output/census.json').read_text())
def need(x,msg):
 if not x: raise RuntimeError(msg)
checks=0
for c in census:
 n=c['n']; nodes=list(csv.DictReader((p/f'output/nodes_{n}.tsv').open(),delimiter='\t'))
 need(len(nodes)==c['ordinary_classes'],'node count')
 reps=set();merged={};hist={}
 for row in nodes:
  r,t,m,s=map(int,row.values());need(r not in reps,'duplicate representative');reps.add(r);merged[m]=merged.get(m,0)+t;hist[t]=hist.get(t,0)+1
 need(sum(merged.values())==c['tilings'],'cube mass')
 need(len(merged)==c['merged_classes'],'merged nodes')
 need(hist=={int(k):v for k,v in c['tc_histogram'].items()},'tc distribution')
 need(sum(int(row['self_converse']) for row in nodes)==c['self_converse'],'SC')
 file=p/f'output/lines_{n}.tsv';op=gzip.open if n==8 else open
 if n==8:file=file.with_suffix('.tsv.gz')
 counts=[0,0];loop=0;degrees=dict.fromkeys(merged,0);seen=set()
 with op(file,'rt') as f:
  for row in csv.DictReader(f,delimiter='\t'):
   u,v,b,m=map(int,row.values());need(u<=v and b in (0,1) and m>0,'line format');need((u,v,b) not in seen,'duplicate line row');seen.add((u,v,b));need(u in merged and v in merged,'line endpoint');degrees[u]+=m;degrees[v]+=m;counts[b]+=m;loop+=m*(u==v);checks+=1
 need(counts==[c['black_lines'],c['blue_lines']],'colours');need(loop==c['loop_lines'],'loops');need(len(seen)==c['distinct_coloured_endpoint_pairs'],'pairs')
 if n>=3:need(degrees==merged,'half-edge degrees');need(sum(counts)==c['tilings']//2,'pairing')
print(json.dumps({'status':'PASS','orders':len(census),'line_rows':checks,'ordinary_nodes':sum(c['ordinary_classes'] for c in census),'n8_line_gzip_sha256':hashlib.sha256((p/'output/lines_8.tsv.gz').read_bytes()).hexdigest()},sort_keys=True))
