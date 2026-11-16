# The Collatz pair 27,63 has a literal coalescing cylinder

Use the positive shortcut T(n)=n/2 for even n and (3n+1)/2 for odd n.
The two paths are

    27,41,62,31,47,71,107,161,242,121,182
    63,95,143,215,323,485,728,364,182

and meet there for the first time. Under the unshortened rule their
first common point is 364, at steps 17 and 13, respectively.
The shortcut words are 1101111101 and 11111100. Shifting a start by
a multiple of 2^length keeps its entire parity word, with exact
output shift 3^(number of odd steps) times the input multiplier.
Consequently the two positive cylinders can be *paired*:

    u=27+1024j, v=63+2304j, j≥0
    T^10(u)=T^8(v)=182+6561j.

Indeed the second shift is 2^8·9j and there are 6 odd steps
rather than 8. In these coordinates v=(9u+9)/4=T²(u)+1.
Thus a genuine size-1 offset after two shortcut steps has a
controlled exact return for this family, with input address density
1/1024 on the u side. This is a coalescence theorem, not a claim that
all endpoints on that progression reach 1. The 1024 in the input
period is required by the ten parity steps, not inferred from the
decimal identity. The graph-count rewiring and the lack of
factored 27→63 containment are distinct constructions.
