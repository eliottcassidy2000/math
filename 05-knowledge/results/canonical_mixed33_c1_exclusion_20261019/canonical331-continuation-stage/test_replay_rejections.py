"""Negative checks for explicit split-tree completeness and exact margins."""
from copy import deepcopy
from pathlib import Path
import json
import exact_replay_engine as R
ROOT=Path(__file__).resolve().parent
record=json.loads((ROOT/'local_repair.json').read_text());p=(0,28,0,0,12,4,13,93,125)
def first_tree(r):return next(c['tree'] for c in r['certificates'] if c['kind']=='slack_one_tree')
def first_leaf(t):return t if t.get('leaf') else first_leaf(t['A_below'])
cases=[]
for name,mutation in [
 ('omitted_child',lambda t:t.pop('B_above')),
 ('wrong_slack_vertex',lambda t:t.update(split_vertex=0)),
 ('dropped_generator',lambda t:first_leaf(t)['generator_indices'].pop()),
 ('wrong_forced_closure',lambda t:first_leaf(t)['forced_A'].clear()),
 ('false_margin',lambda t:first_leaf(t).update(max_margin=1)),
 ('negative_weight',lambda t:first_leaf(t)['terms'][0].__setitem__(2,-1)),
]:
 r=deepcopy(record);mutation(first_tree(r))
 try:R.replay(p,r)
 except (RuntimeError,KeyError):cases.append(name)
 else:raise RuntimeError('Accepted corrupted certificate: '+name)
R.replay(p,record)
print(json.dumps({'status':'ALL_NEGATIVE_REPLAY_CHECKS_PASSED','rejected':cases,'positive_control':'passed'}))
