#!/usr/bin/env python3
from pathlib import Path
import sys,json
sys.path.insert(0,str(Path(__file__).resolve().parent.parent/'square_sum_spectrum_overnight_20261031j'))
from audit_tree import cycles,need
A={(0,1),(1,2),(0,2),(3,4),(4,5),(3,5),(0,3),(1,4),(2,5)}
B={(6,7),(7,8),(8,9),(6,9)}
portsA={(0,1),(3,4)};portsB={(6,7),(8,9)};seams={(0,6),(1,7),(3,8),(4,9)}
G=(A-portsA)|(B-portsB)|seams
actual=cycles(G);local=[c for c in cycles(A) if portsA<=c]
products={frozenset((a-portsA)|(b-portsB)|seams) for a in local for b in cycles(B)}
need(len(actual)==3 and len(local)==2 and len(products)==2 and products<actual,'strict graph versus product gap')
extra=next(iter(actual-products))
reclosed={e for e in extra if max(e)<6}|portsA
need(reclosed=={(0,1),(1,2),(0,2),(3,4),(4,5),(3,5)},'two local cycles after restoring ports')
print(json.dumps({'status':'passed','global_cycles':3,'source_A_port_cycles':2,'source_B_cycles':1,
 'selected_cycle_products':2,'extra_global_cycle':sorted(extra),'reclosed_local_edges':sorted(reclosed),
 'scope':'Exhaustive graph-only scope correction; no square-sum realization asserted.'},indent=2,sort_keys=True))
