# Cyclic composition of transitive tournaments: exact path enumerator

For nonnegative integers a,b,c, take three ordered transitive tournaments A,B,C of those sizes and orient all interblock edges A→B→C→A. Write H(a,b,c) for the number of directed Hamiltonian paths, with H(0,0,0)=1 for the empty path. If any block is empty the tournament is transitive, so H=1.

Put u=e^x−1, v=e^y−1, w=e^z−1. The trivariate exponential generating function is

    Σ H(a,b,c) x^a y^b z^c / (a! b! c!)
      = (e^(x+y+z)+uvw)/(1−uvw).                    (1)

In particular, for all nonnegative a,b,c,

    H(a,b,c) = 1 + Σ_{k=1}^{min(a,b,c)} (k!)^3
       [ S(a+1,k+1) S(b+1,k+1) S(c+1,k+1)
         + S(a,k) S(b,k) S(c,k) ],                    (2)

where S(n,k) is a Stirling number of the second kind. The convention makes an empty sum zero.

Proof. Uniquely decompose a nonempty directed path into maximal runs within a block. Its block labels cycle A,B,C,A,..., because each cross-block arc must follow the cyclic orientation. A run with j chosen vertices has exactly one permitted internal order and EGF x^j/j!, hence u counts an A-run, and similarly v,w. Summing all cyclic words of runs gives

    1 + (u+v+w+uv+vw+wu+3uvw)/(1−uvw)
      = (e^(x+y+z)+uvw)/(1−uvw).

The 3uvw term counts the three possible starts of a full turn. The empty path gives the first 1. This run decomposition also covers empty blocks; the corresponding variable is zero. Expand the right-hand side of (1) as e^(x+y+z)+Σ_{k≥1}(e^(x+y+z)+1)(uvw)^k. The coefficients satisfy [x^a/a!] (e^x−1)^k=k!S(a,k) and [x^a/a!] e^x(e^x−1)^k=k!S(a+1,k+1). Extracting coefficients yields (2). No assumption that vertex labels preserve order between separate runs is made: choices of which block vertices lie in each run are essential.

Three exact one-parameter consequences for a≥1 are

    H(a,1,1) = 2^a+1,
    H(a,2,1) = 3·2^a−1,
    H(a,2,2) = 4·3^a+5·2^a−11.

Each formula requires a≥1; the zero-block case equals 1. The first is the existing almost-transitive family, not a novelty claim. The derived formulas furnish explicit witnesses but do not settle the global spectrum. These constructions with positive blocks are strongly connected, and all can be padded with a source or sink without changing H.

Verification: `python check_composition.py` independently constructs each block tournament and uses a subset/last-vertex path DP, comparing against Stirling formula (2) for all 220 nonnegative triples of total order at most nine, including the empty case and all degenerate two-block cases. It checks the three displayed families whenever their blocks are positive. See check.out. The test is finite corroboration of the proof, not an exhaustive spectrum claim.

Prior-art scope: Repository search of canon/results for both cyclic/transitive-block and Stirling/Hamiltonian terms did not locate this exact EGF or Stirling expression. The repository already records H(a,1,1)=2^a+1 and broad blow-up/path-cover machinery; no claim of mathematical novelty or priority is made. No external literature priority search was performed.
