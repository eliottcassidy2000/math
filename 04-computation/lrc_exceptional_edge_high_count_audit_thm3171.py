#!/usr/bin/env python3
"""Exact ratio-graph audit for the sole weak high-channel edge.

The two-star pivots carry labels 1 and 2.  The only high channel whose
proved exact overlap floor is below 1/105 is the label-{1,6} edge on the
primitive channel 3:5.  Normalize the label-1 pivot level to one; its
label-6 leaf is therefore x=3/5 or 5/3.  This referee proves that fixing
that leaf raises the high-edge count from three to at least four.
"""
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations

LOW=frozenset((F(1,2),F(2,3),F(3,4),F(4,3),F(3,2),F(2)))
QUOTIENTS=frozenset(x/y for x in LOW for y in LOW
                    if F(1,2)<=x/y<=2)
EXPECTED_QUOTIENTS=frozenset((F(1,2),F(9,16),F(2,3),F(3,4),F(8,9),
                              F(1),F(9,8),F(4,3),F(3,2),F(16,9),F(2)))
EXCEPTIONAL=(F(3,5),F(5,3))
EDGE_COUNT=9

def require(condition,detail):
 if not condition:raise RuntimeError(detail)

def degree(t,x):
 return int(x in LOW)+int(x/t in LOW)

def inside_row(t,x):
 # Only positive-degree optional leaves can improve the low-edge count.
 candidates=tuple(sorted((set(LOW)|{t*s for s in LOW})-{F(1),t,x}))
 best=None
 for size in range(4):
  for chosen in combinations(candidates,size):
   levels=(F(1),t,x,*chosen)
   if len(set(levels))!=len(levels) or max(levels)>2*min(levels):continue
   low=int(t in LOW)+degree(t,x)+sum(degree(t,y) for y in chosen)
   row=(low,chosen)
   if best is None or row>best:best=row
 if best is None:return None
 return EDGE_COUNT-best[0],best,candidates

def main():
 require(QUOTIENTS==EXPECTED_QUOTIENTS,(QUOTIENTS,EXPECTED_QUOTIENTS))
 require(LOW<=QUOTIENTS,(LOW,QUOTIENTS))
 require(all(x not in LOW and x.numerator+x.denominator==8
             for x in EXCEPTIONAL),EXCEPTIONAL)
 rows=[]
 for x in EXCEPTIONAL:
  for t in sorted(QUOTIENTS-{F(1),x}):
   row=inside_row(t,x)
   if row is not None:rows.append((x,t,*row))
 # If t is outside LOW/LOW, the pivot edge is high and no leaf is low to
 # both pivots.  Thus at most four of the four leaf incidences are low,
 # giving the stronger outside bound of five high edges.
 inside_min=min(row[2] for row in rows)
 outside_min=5
 require(inside_min==4,("inside minimum",inside_min))
 require(outside_min>=4,outside_min)
 sharp=tuple(row for row in rows if row[2]==inside_min)
 semantic=repr((tuple(sorted(LOW)),tuple(sorted(QUOTIENTS)),rows,
                inside_min,outside_min,sharp)).encode()
 print("LRC EXCEPTIONAL 3:5 EDGE -- EXACT HIGH-COUNT AUDIT")
 print(f"low_ratios={tuple(sorted(LOW))}")
 print(f"pivot_quotient_rows={len(QUOTIENTS)-1};exceptional_orientations={EXCEPTIONAL}")
 print(f"inside_rows={len(rows)};inside_high_min={inside_min};outside_high_min={outside_min}")
 print(f"sharp_rows={sharp}")
 print("consequence=exceptional edge implies at least three other regular high edges")
 print(f"semantic_sha256={sha256(semantic).hexdigest()}")

if __name__=="__main__":main()
