"""Every arc from the published strong sources at bounded offsets from E_61."""
from seam_quarter import arc_changes, inflate, transitive_blocks, two_slot_kernel
from seam_hosts import canonical, decode
from math import comb
from itertools import permutations, product

# offset, upper triangular word, marked sites, width shifts relative to p-1
RAILS = [
    (-5, "111110111111111111011101011111111101", (1,8), (-2,-3)),
    (-9, "110111111101111110111111111110101111", (0,8), (-4,-2)),
    (-15, "110111111101101111010111111111111111", (0,8), (-3,-3)),
    (-41, "101011001111111111111111110111111111111111110", (2,8), (-3,-3)),
]
OLD_RAILS = [
    (-5,(990,1020,1016,976,928,137,688,768,608,32),(1,9),(57,57)),
    (-9,(998,1020,1016,977,929,456,656,768,576,32),(1,9),(55,58)),
    (-15,(974,1020,1016,976,929,201,528,832,608,32),(1,9),(56,57)),
]


def arc_census(e, rows, ports, widths):
    kernel=two_slot_kernel(rows,ports)
    changes=list(arc_changes(rows,ports,widths,kernel))
    n=len(rows)-len(ports)+sum(widths)
    assert len(changes)==comb(n,2)
    near=[(edge,d,e+d) for edge,d in changes if d and abs(e+d)<=200]
    positive=min(d for _,d in changes if d>0)
    negative=max(d for _,d in changes if d<0)
    _,_,maximal_widths,_=transitive_blocks(inflate(rows,ports,widths))
    assert not near and sum(d==0 for _,d in changes)==sum(w-1 for w in maximal_widths)
    print(e, n, len(changes), "closest changes",negative,positive,
          "neutral",sum(d==0 for _,d in changes),
          "near",near,flush=True)


def color_code(graph, ports):
 n=len(graph)
 cells=[tuple(i for i in range(n) if i not in ports),tuple(ports)]
 while True:
  groups=[]
  for cc in cells:
   d={}
   for i in cc:
    k=tuple(sum(graph[i][j] for j in c) for c in cells)
    d.setdefault(k,[]).append(i)
   groups.extend(d[k] for k in sorted(d))
  if len(groups)==len(cells):break
  cells=groups
 return min(tuple(graph[i][j] for i in p for j in p) for ordering in product(*(permutations(c) for c in cells)) if (p:=sum(ordering,())))


def main():
    for e, word, ports, shifts in RAILS:
        rows=tuple(sum(b<<j for j,b in enumerate(row)) for row in decode(word))
        board = decode(word)
        n = len(board)
        parent_codes, split_codes = set(), set()
        for g in (board, list(zip(*board))):
            parent_codes.add(color_code(g, ports))
            for a in ports:
                labels = list(range(n)) + [a]
                h = [[int(i < j) if u == v else g[u][v] for j, v in enumerate(labels)]
                     for i, u in enumerate(labels)]
                for t in (a, n):
                    marked = [t if v == a else v for v in ports]
                    split_codes.add(color_code(h, marked))
        assert len(parent_codes) == 2 and len(split_codes) == 8
        sizes = [1]*n
        for site, shift in zip(ports, shifts):
            sizes[site] = 60 + shift
        assert canonical(board, sizes) != canonical(list(zip(*board)), sizes)
        arc_census(e, rows, ports, tuple(60+s for s in shifts))
    for entry in OLD_RAILS:
        arc_census(*entry)


if __name__=="__main__":
    main()
