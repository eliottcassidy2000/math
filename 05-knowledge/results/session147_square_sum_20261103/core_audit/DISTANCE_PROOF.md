# A sharp distance threshold for 2 and 4

In the simple positive-integer square-sum graph, 2 and 4 are not adjacent since 2+4=6 is not a square. Nor can they have a common neighbor w: (w+4)−(w+2)=2, but a difference of integer squares is never 2 modulo 4. Thus any path between 2 and 4 has length at least 3, for every n.

They first become connected in G_14. Direct enumeration of G_13 gives the two components {1,3,4,5,6,8,10,11,12,13} and {2,7,9}. In G_14, the path 2–14–11–5–4 has square sums 16,25,16,9. As the graphs are nested, the pair is disconnected for n≤13 and their distance is at most 4 for n≥14.

Suppose a 3-edge path 2–u–v–4 exists with u,v≤76. There are positive integers s,t,r with

    u=s²−2,  v=t²−4,  r²=u+v=s²+t²−6.

Here s≥2 and t≥3 because u,v are positive. Since u,v≤76, s≤8 and t≤8. Modulo 4, the equation forces s and t both odd and r even: odd/even and even/odd give r²≡3, and even/even gives r²≡2. Therefore s,t∈{3,5,7}. The six possible values of s²+t²−6, accounting for symmetry, are 12, 28, 44, 52, 68, and 92, none of which is a square. Hence no 3-edge path exists in G_76 or any smaller G_n.

In G_77, the distinct positive vertices 2–23–77–4 form a 3-edge path with sums 25,100,81. Equivalently s=5, t=9, r=10. No shorter path is possible by the first paragraph. Therefore the distance from 2 to 4 is infinite for n≤13, exactly 4 for 14≤n≤76, and exactly 3 for n≥77.

The same common-neighbor argument shows that 1 and 3 have none. The pair 4,5 has none either: their common neighbor would give two positive squares whose difference is 1, which is impossible. Thus neither 1–3 nor 4–5 lies in a triangle, even though both are edges.
