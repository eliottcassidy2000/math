"""Finite shadow of the official RH machine in the inherited marked tournament tape."""
import hashlib
from pathlib import Path
from sys import argv

raw=Path(argv[1]).read_bytes()
assert hashlib.sha256(raw).hexdigest()=='7e9731ed1b7c63d4e67c0f8d646f8498d9645f26e237a8604d4b7f91a6ab398d'
table=tuple(tuple(s.split()) for s in raw.decode('ascii').splitlines())
assert len(table)==95 and all(len(row)==2 for row in table)
assert [(q,s) for q,row in enumerate(table) for s,t in enumerate(row) if t=='---']==[(61,1)]
duration=int(argv[2]) if len(argv)>2 else 3634
head=state=0; ones=set(); trace=[]; lo=hi=0
for step in range(duration):
    t=table[state][int(head in ones)]
    if t=='---': raise AssertionError(('halt',step,state,head))
    glyph,motion,dest=t[0],t[1],int(t[2:])
    tgt=head+(-1 if motion=='L' else 1)
    trace.append((state,head,int(head in ones),dest,tgt,int(glyph)))
    if glyph=='1': ones.add(head)
    else: ones.discard(head)
    state,head=dest,tgt
    lo=min(lo,head); hi=max(hi,head)
if duration==3634: assert (state,head)==(15,-4)

# Each bit is one five-vertex tile block, differing on only the {0,2} arc.
# This replay uses sparse active blocks, not a full adjacency matrix.
def bit(cell,field):return (cell-lo)*96+field
live={bit(0,1)}
arc_cost=0
for q,x,s,r,y,t in trace:
    assert bit(x,q+1) in live and (bit(x,0) in live)==bool(s)
    changes=[bit(y,r+1),bit(x,q+1)]
    if t!=s: changes.insert(0,bit(x,0))
    assert len(changes)==len(set(changes))<=3
    for c in changes:
        if c in live:live.remove(c)
        else:live.add(c)
    arc_cost+=len(changes)
assert live=={bit(head,state+1)}|{bit(x,0) for x in ones}
L=hi-lo+1
print('native steps',duration,'state',state,'head',head,'visited interval',lo,hi)
print('tournament order',480*L,'arc flips',arc_cost,'fixed height 9^'+str(96*L))
