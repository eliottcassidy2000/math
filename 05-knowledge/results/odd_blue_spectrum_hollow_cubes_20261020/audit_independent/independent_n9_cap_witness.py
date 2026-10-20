#!/usr/bin/env python3
"""Direct full-S9 witness against the historical 3^floor((n-2)/2) cap."""
from independent_pairing_audit import *
A,B=make_e(4,194);G=graph(A,B,0);exact=direct(G)
assert exact=={'H':2659,'aut':1,'fibre':2659,'blue':39,'anti_involutions':1,'centralizer':1}
row={'n':9,'m':4,'e':194,'b':0,'h_r':len(sym_paths(G)),'hollow_group':1,'adjacency_row_masks':[sum(x<<i for i,x in enumerate(r)) for r in G],'c3':c3(G),'exact_full_S9':exact,'historical_cap':27,'cap_refuted':True}
assert row['h_r']==39 and row['c3']==28
print(json.dumps(row));(OUT/'independent_n9_cap_witness.json').write_text(json.dumps(row,indent=2)+'\n')
