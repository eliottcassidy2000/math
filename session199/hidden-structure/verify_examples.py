"""Replay the explicit large-target recipes without a large source census."""
import hashlib,json
from pathlib import Path
from verify import T,replay,read,graph,require
here=Path(__file__).resolve().parent
examples=json.loads((here/'EXAMPLES.json').read_text())
for r in examples:
 w=r['source_word_from_transitive_edge'];s=replay(T,w);n=r['target']
 require(list(s)==r['source_profile'] and read(s)==n,'source profile and target')
 require(w.count('J')==r['root_changes']<=3,'bounded root changes')
 require(r['balanced_word_from_transitive_edge']==w+'OJ','balanced wrapper')
 require(list(replay(T,w+'OJ'))==r['balanced_profile']==[n,s[0],n,s[0]],'balanced profile')
 out,root=graph(w+'OJ')
 require(len(out)==r['balanced_graph_order'],'literal graph order')
 # Each birth follows the exact physical source rule by graph(). No costly
 # unrestricted subset path count is claimed for these order22/23 examples.
 require(out[root].bit_count()==1,'balanced root has one out-neighbor')
print(json.dumps({'status':'PASS','targets':[r['target']for r in examples],
 'check':'exact physical recipe replay and graph construction; not independent subset DP',
 'examples_sha256':hashlib.sha256((here/'EXAMPLES.json').read_bytes()).hexdigest()},indent=2))
