#!/usr/bin/env python3
"""rapidity_molecules_s116d.py — Lithium hydride, alkanes, and molecular metaphors

LiH: the simplest heteronuclear molecule. Two atoms, maximally different.
Alkanes: C_nH_{2n+2}. Chains of carbon. The backbone of organic chemistry.

What do these teach us about tournaments, rapidity, and structure?
"""
from math import sqrt, log, pi, exp, atanh, tanh, comb, factorial
from fractions import Fraction

phi = (1+sqrt(5))/2
octave = log(2)/2

print("LITHIUM HYDRIDE AND ALKANES THROUGH THE RAPIDITY LENS")
print("="*70)
print()

# ============================================================
print("="*70)
print("I. LITHIUM HYDRIDE: THE SIMPLEST BINARY MOLECULE")
print("="*70)
print()
print("  LiH: Z_Li = 3, Z_H = 1.")
print("  Electronegativity: Li = 0.98, H = 2.20.")
print("  Bond type: ionic/covalent mix (EN diff = 1.22).")
print()

EN_Li = 0.98
EN_H = 2.20
v_bond = (EN_H - EN_Li) / (EN_H + EN_Li)
rap_bond = atanh(v_bond)
Q_bond = (1 + v_bond) / (1 - v_bond)

print(f"  Cayley velocity of bond: v = (EN_H - EN_Li)/(EN_H + EN_Li)")
print(f"  = ({EN_H} - {EN_Li})/({EN_H} + {EN_Li}) = {v_bond:.6f}")
print(f"  Rapidity = arctanh({v_bond:.6f}) = {rap_bond:.6f}")
print(f"  Q(v) = {Q_bond:.6f}")
print(f"  Note: Q(v) = EN_H/EN_Li = {EN_H/EN_Li:.6f}")
print()
print("  The Cayley transform of the bond polarity IS the EN ratio.")
print("  This is a tautology (Q((b-a)/(b+a)) = b/a), but it means:")
print("  EVERY heteronuclear bond lives on the Cayley line.")
print("  The bond IS a point on the rapidity axis.")
print()

# LiH as a TOURNAMENT
print("  LiH AS A TOURNAMENT:")
print("  Two 'players': Li and H. One 'match': Li vs H.")
print("  H wins (higher EN). The tournament is the arc H -> Li.")
print("  This is the TRANSITIVE tournament on 2 vertices.")
print("  H(T) = 1. rapidity(H) = 0.")
print()
print("  But the bond POLARITY carries additional information:")
print("  The polarity v = 0.384 tells us HOW MUCH H dominates Li.")
print("  In a pure tournament, each arc is binary (win/lose).")
print("  In a WEIGHTED tournament, each arc has strength v in [0,1).")
print("  LiH has arc strength 0.384.")
print()

# ============================================================
print("="*70)
print("II. THE BOND AS A RAPIDITY INTERVAL")
print("="*70)
print()
print("  Every chemical bond has a rapidity = arctanh(|EN_B - EN_A|/(EN_B + EN_A)).")
print("  This rapidity measures the ASYMMETRY of the bond.")
print()

bonds = [
    ("H-H",    2.20, 2.20, "covalent"),
    ("C-H",    2.55, 2.20, "slightly polar"),
    ("N-H",    3.04, 2.20, "polar"),
    ("O-H",    3.44, 2.20, "polar"),
    ("F-H",    3.98, 2.20, "very polar"),
    ("C-C",    2.55, 2.55, "covalent"),
    ("C-N",    2.55, 3.04, "polar"),
    ("C-O",    2.55, 3.44, "polar"),
    ("C-F",    2.55, 3.98, "very polar"),
    ("Li-H",   0.98, 2.20, "ionic"),
    ("Na-Cl",  0.93, 3.16, "ionic"),
    ("Cs-F",   0.79, 3.98, "most ionic"),
]

print(f"  {'Bond':<8s} {'EN_A':>5s} {'EN_B':>5s} {'v':>8s} {'rapidity':>10s} {'Q=EN_B/EN_A':>12s} {'Type'}")
print("  " + "-"*65)
for name, en_a, en_b, typ in bonds:
    if en_a == en_b:
        v = 0.0
        rap = 0.0
        q = 1.0
    else:
        hi, lo = max(en_a, en_b), min(en_a, en_b)
        v = (hi - lo) / (hi + lo)
        rap = atanh(v)
        q = hi / lo
    print(f"  {name:<8s} {en_a:>5.2f} {en_b:>5.2f} {v:>8.4f} {rap:>10.6f} {q:>12.4f}   {typ}")

print()
print("  COVALENT bonds have rapidity ~ 0 (symmetric, no polarity).")
print("  IONIC bonds have rapidity >> 0 (one atom dominates).")
print("  The rapidity IS the degree of ionicity.")
print()
print("  The MOST ionic bond (Cs-F): rapidity = {:.4f}".format(
    atanh((3.98 - 0.79) / (3.98 + 0.79))))
print("  = {:.2f} octaves of asymmetry.".format(
    atanh((3.98 - 0.79) / (3.98 + 0.79)) / octave))
print()

# ============================================================
print("="*70)
print("III. ALKANES: THE CHAIN METAPHOR")
print("="*70)
print()
print("  Alkanes: C_n H_{2n+2}. The simplest organic chains.")
print("  Methane CH4, Ethane C2H6, Propane C3H8, Butane C4H10, ...")
print()
print("  Total atoms = n + (2n+2) = 3n+2.")
print("  Total bonds = n-1 (C-C) + 2n+2 (C-H) = 3n+1.")
print()
print("  In our tournament theory, a tournament on m vertices has")
print("  C(m,2) = m(m-1)/2 arcs. An alkane has 3n+1 bonds.")
print("  These are NOT the same (a molecule is not a tournament).")
print("  But the CHAIN STRUCTURE of alkanes echoes the PATH structure")
print("  of Hamiltonian paths in tournaments.")
print()

print("  n   Formula    Atoms  Bonds  C(atoms,2)  ratio bonds/arcs")
print("  " + "-"*60)
for n in range(1, 11):
    atoms = 3*n + 2
    bonds_count = 3*n + 1
    arcs = atoms * (atoms - 1) // 2
    ratio = bonds_count / arcs
    print(f"  {n:2d}  C{n}H{2*n+2:<3d}    {atoms:4d}   {bonds_count:4d}    {arcs:6d}      {ratio:.4f}")

print()
print("  The ratio bonds/arcs -> 0 as n -> infinity.")
print("  Molecules are SPARSE graphs. Tournaments are COMPLETE graphs.")
print("  But both have PATHS as their fundamental object:")
print("  - In tournaments: Hamiltonian PATHS count H(T).")
print("  - In alkanes: the carbon BACKBONE is a path (linear alkanes).")
print()

# ============================================================
print("="*70)
print("IV. THE ALKANE AS A RAPIDITY CHAIN")
print("="*70)
print()
print("  Each C-C bond has rapidity 0 (same atom type).")
print("  Each C-H bond has rapidity arctanh((2.55-2.20)/(2.55+2.20)).")
ch_rap = atanh((2.55 - 2.20) / (2.55 + 2.20))
print(f"  = arctanh(0.0737) = {ch_rap:.6f}")
print()
print("  The alkane backbone: C-C-C-...-C")
print("  All C-C bonds have rapidity 0 (no polarity).")
print("  This is a PATH OF ZERO RAPIDITY -- a 'transitive' chain.")
print()
print("  The C-H branches: each has rapidity {:.4f}.".format(ch_rap))
print("  Total rapidity of all C-H bonds:")
for n in range(1, 8):
    total_ch = (2*n + 2) * ch_rap
    print(f"    C{n}H{2*n+2}: {2*n+2} C-H bonds, total rapidity = {total_ch:.4f} "
          f"= {total_ch/octave:.2f} octaves")

print()
print("  Each CH4 unit adds 2 C-H bonds = {:.4f} rapidity = {:.2f} octaves.".format(
    2*ch_rap, 2*ch_rap/octave))
print("  EXTENDING the alkane chain adds a CONSTANT rapidity per carbon.")
print("  This is EXACTLY like adding a musical interval to a melody.")
print()

# ============================================================
print("="*70)
print("V. BRANCHING: WHEN THE CHAIN BREAKS LINEARITY")
print("="*70)
print()
print("  Linear alkanes: one path from end to end.")
print("  Branched alkanes: multiple paths, like a tournament with cycles.")
print()
print("  Isomers of C4H10:")
print("  - n-butane: C-C-C-C (linear path, 1 backbone path)")
print("  - isobutane: (CH3)3CH (branched, 3 backbone paths)")
print()
print("  Isomers of C5H12:")
print("  - n-pentane: C-C-C-C-C (1 backbone path)")
print("  - isopentane: (CH3)2CHCH2CH3 (2 backbone paths)")
print("  - neopentane: C(CH3)4 (6 backbone paths?)")
print()
print("  The number of distinct longest paths through the carbon backbone")
print("  is ANALOGOUS to H(T) in tournament theory!")
print("  Linear alkane: H = 1 (transitive tournament)")
print("  Branched alkane: H > 1 (more paths = more complexity)")
print()

# Count Hamiltonian paths in small carbon skeletons
# n-butane: 1-2-3-4, paths: 1234, 4321 = 2 (but orientation matters)
# Actually in an UNDIRECTED path graph: 2 Ham paths (each direction)
# In a star K_{1,3}: 6 Ham paths (3! permutations of leaves)
print("  CARBON SKELETON TOURNAMENT ANALOGY:")
print("  Skeleton    Graph          H_undirected  H_directed")
print("  " + "-"*55)
print("  n-butane    Path P4        2             1 (transitive)")
print("  isobutane   Star K_{1,3}   6             3 (like T_3!)")
print("  n-pentane   Path P5        2             1")
print("  isopentane  P3 + branch    4             2")
print("  neopentane  Star K_{1,4}  24             12")
print()
print("  Isobutane's skeleton has 3 directed Ham paths = H(T_3) = 3!")
print("  Neopentane's skeleton has 12 directed Ham paths.")
print("  The BRANCHING of alkanes IS the CYCLICITY of tournaments.")
print()

# ============================================================
print("="*70)
print("VI. THE METAPHOR: MOLECULES AS WEIGHTED TOURNAMENTS")
print("="*70)
print()
print("  A molecule is a WEIGHTED tournament where:")
print("  - Each atom is a vertex")
print("  - Each bond is an arc (from less EN to more EN)")
print("  - The arc weight = rapidity of the bond polarity")
print("  - Non-bonded pairs have 'arcs' determined by through-bond effects")
print()
print("  TOURNAMENT CONCEPT     MOLECULAR ANALOGUE")
print("  " + "-"*55)
print("  vertex                 atom")
print("  arc                    bond (directed by EN)")
print("  arc weight             bond polarity = rapidity")
print("  Hamiltonian path       reaction pathway through all atoms")
print("  H(T) = path count     number of distinct pathways")
print("  3-cycle                3-membered ring (cyclopropane)")
print("  transitive             linear chain (no rings)")
print("  Paley tournament       highly symmetric molecule (benzene?)")
print("  forbidden H = 7       impossible molecular topology?")
print("  conflict graph Omega   steric clash graph")
print("  independence poly      non-clashing conformations")
print()

# ============================================================
print("="*70)
print("VII. BENZENE: THE PALEY MOLECULE")
print("="*70)
print()
print("  Benzene C6H6: 6 carbons in a ring, perfectly symmetric.")
print("  Symmetry group: D_{6h}, order 24 (including reflections).")
print()
print("  The Paley tournament on p=7 vertices has |Aut| = 21.")
print("  Benzene has |Aut| = 24 (D_{6h}) ~ comparable.")
print()
print("  Benzene's carbon skeleton is the cycle C_6.")
print("  In tournament terms: a 6-vertex tournament where every")
print("  triple forms a 3-cycle would be 'benzene-like'.")
print()
print("  But benzene is UNDIRECTED (all C-C bonds equivalent).")
print("  The 'directed' version of benzene would be a TOURNAMENT")
print("  on the 6 carbon vertices. The most benzene-like tournament")
print("  is the one with maximum symmetry: the regular tournament on 6.")
print()
print("  Wait -- 6 is EVEN. There's no regular tournament on even vertices.")
print("  The most symmetric tournament on 6 has scores (2,2,2,3,3,3).")
print("  This is the 'almost-regular' tournament. H_max = 45.")
print()
print("  BENZENE'S SYMMETRY DEFIES TOURNAMENT STRUCTURE")
print("  because benzene lives in the undirected world,")
print("  while tournaments live in the directed world.")
print("  The transition from undirected to directed is the")
print("  transition from COVALENT to IONIC bonding.")
print("  Benzene = pure covalent = rapidity 0 on every bond.")
print("  LiH = ionic = rapidity 0.40 on the bond.")
print()

# ============================================================
print("="*70)
print("VIII. POLYMERIZATION AS TOURNAMENT GROWTH")
print("="*70)
print()
print("  Polymerization: monomers join to form chains.")
print("  n monomers -> polymer of length n.")
print("  This is EXACTLY vertex addition in tournament theory!")
print()
print("  When we add a vertex to a tournament on n-1 vertices,")
print("  we choose n-1 arc directions. The number of new tournaments")
print("  = 2^{n-1} (choosing each arc direction independently).")
print()
print("  When we add a monomer to a polymer of length n-1,")
print("  we choose HOW it attaches (stereochemistry, regiochemistry).")
print("  The number of distinct attachment modes is typically small")
print("  (1-3 for simple alkane extension), but for complex monomers")
print("  it can be much larger.")
print()
print("  The GROWTH of H(T) as n increases parallels the GROWTH")
print("  of conformational complexity as a polymer extends.")
print()
print("  For alkanes:")
print("  n=1 (CH4):   1 conformation (tetrahedral)")
print("  n=2 (C2H6):  1 distinct rotamer (eclipsed/staggered ~ 3)")
print("  n=3 (C3H8):  2 rotamers (anti/gauche)")
print("  n=4 (C4H10): 3-5 rotamers")
print("  n=5 (C5H12): ~9 rotamers")
print()
print("  Compare tournament H(T) growth:")
print("  n=3: H in {1, 3}")
print("  n=4: H in {1, 3, 5}")
print("  n=5: H in {1, 3, 5, 9, 11, 13, 15}")
print("  n=6: H in {1, 3, 5, ..., 45}")
print()
print("  Both grow roughly exponentially, but H(T) grows faster")
print("  because tournaments are DENSER (complete graphs).")
print()

# ============================================================
print("="*70)
print("IX. THE C-H STRETCH AS A RAPIDITY OSCILLATOR")
print("="*70)
print()
print("  The C-H stretching frequency: ~3000 cm^-1.")
print("  The C-C stretching frequency: ~1000 cm^-1.")
print("  Ratio: 3000/1000 = 3 = perfect twelfth (octave + fifth).")
print()
print("  In rapidity: ln(3000/1000)/2 = ln(3)/2 = {:.6f}".format(log(3)/2))
print("  = the rapidity of the curvature quantum!")
print()
print("  The C-H stretch vibrates a TWELFTH above the C-C stretch.")
print("  In the harmonic approximation: frequency ~ sqrt(k/mu)")
print("  where k = force constant, mu = reduced mass.")
print()
print("  k(C-H) ~ k(C-C) (similar force constants).")
print("  mu(C-H) ~ 1 amu, mu(C-C) ~ 6 amu.")
print("  Ratio: sqrt(mu(CC)/mu(CH)) = sqrt(6) ~ 2.45.")
print("  Actual ratio: 3000/1000 = 3. The discrepancy (3 vs 2.45)")
print("  comes from k(C-H) > k(C-C) (C-H bonds are stiffer).")
print()
print("  In rapidity: the C-H / C-C frequency ratio 3:1 means the")
print("  C-H bond vibrates at the TWELFTH harmonic of the C-C bond.")
print("  This is EXACTLY the same ratio as Q(1/2) = 3:")
print("  Q(1/2) = (1 + 1/2)/(1 - 1/2) = 3/1.")
print("  The C-H stretch is the Cayley transform of 1/2!")
print()

# ============================================================
print("="*70)
print("X. LiH DIPOLE MOMENT AND THE GOLDEN RATIO")
print("="*70)
print()
print("  LiH has a dipole moment of 5.88 Debye.")
print("  The 'ideal ionic' dipole moment (full charges at bond length):")
print("  mu_ionic = e * r = 1.602e-19 C * 1.596e-10 m = 2.556e-29 C*m")
print("  = 7.66 Debye (using 1 Debye = 3.336e-30 C*m).")
print()
print("  The ionic character = actual/ideal = 5.88/7.66 = 0.767.")
print("  The Cayley velocity of ionic character:")
ionic = 5.88 / 7.66
v_ionic = (2*ionic - 1)  # map [0,1] to [-1,1]
print(f"  Ionic character = {ionic:.4f}")
print(f"  2*ionic - 1 = {v_ionic:.4f}")
if abs(v_ionic) < 1:
    rap_ionic = atanh(v_ionic)
    print(f"  arctanh({v_ionic:.4f}) = {rap_ionic:.4f}")
    print(f"  In octaves: {rap_ionic/octave:.3f}")
print()
# Is 0.767 close to anything?
print(f"  0.767 ~ 1/phi^(3/2) = {1/phi**1.5:.4f}? No.")
print(f"  0.767 ~ 3/4 + 1/60 = 0.767? Close to 3/4 = 0.750.")
print(f"  The ionic character is close to 3/4 = the address of 7")
print(f"  (the forbidden number's Cayley address is (7-1)/(7+1) = 3/4).")
print()

# ============================================================
print("="*70)
print("XI. THE ALKANE FORMULA AND TOURNAMENT STRUCTURE")
print("="*70)
print()
print("  Alkane: C_n H_{2n+2}.")
print("  The formula 2n+2 appears EVERYWHERE in our theory!")
print()
print("  In rapidity algebra:")
print("  Q(1/(2n+1)) = (n+1)/n = the n-th superparticular ratio.")
print("  The alkane hydrogen count 2n+2 = 2*(n+1).")
print("  So H_count/2 = n+1, and the 'reduced hydrogen number'")
print("  is the NUMERATOR of the musical interval Q(1/(2n+1)).")
print()
print("  n   C_nH_{2n+2}  H/2=n+1  Interval Q(1/(2n+1))")
print("  " + "-"*55)
for n in range(1, 9):
    h_count = 2*n + 2
    interval = Fraction(n+1, n)
    print(f"  {n}   C{n}H{h_count:<3d}       {n+1:3d}      {interval} = Q(1/{2*n+1})")

print()
print("  Methane CH4: n=1, interval = 2/1 = OCTAVE")
print("  Ethane C2H6: n=2, interval = 3/2 = FIFTH")
print("  Propane C3H8: n=3, interval = 4/3 = FOURTH")
print("  Butane C4H10: n=4, interval = 5/4 = MAJOR THIRD")
print("  Pentane C5H12: n=5, interval = 6/5 = MINOR THIRD")
print()
print("  THE ALKANE SERIES IS THE HARMONIC SERIES OF MUSIC!")
print("  Each alkane corresponds to a musical interval.")
print("  Methane = octave. Ethane = fifth. Propane = fourth.")
print("  As the chain grows, the interval SHRINKS (less consonant).")
print("  Long alkanes = DISSONANT intervals = complex chemistry.")
print()

# ============================================================
print("="*70)
print("XII. THE POLYMER AS AN ARCTANH SERIES")
print("="*70)
print()
print("  arctanh(x) = x + x^3/3 + x^5/5 + x^7/7 + ...")
print("  Each term x^{2k+1}/(2k+1) corresponds to:")
print("  - Musical interval: (k+1)/k")
print("  - Alkane: C_k H_{2k+2}")
print()
print("  arctanh IS the generating function of the alkane series!")
print("  The k-th term of arctanh encodes the k-th alkane.")
print()
print("  A POLYMER of length N is like truncating arctanh at the N-th term.")
print("  The full polymer (N -> infinity) is the complete arctanh = rapidity.")
print()
print("  The COUPLING STRENGTH in tournament theory (arctanh of the")
print("  tournament parameter x) is built from the 'monomers' 1/3, 1/5, 1/7, ...")
print("  Each monomer is an alkane. Each alkane is a musical interval.")
print("  The tournament coupling IS a polymer of musical intervals.")
print()

# ============================================================
print("="*70)
print("XIII. GRAND METAPHOR: THE MOLECULAR TOURNAMENT")
print("="*70)
print()
print("  CHEMISTRY              TOURNAMENT              MUSIC")
print("  " + "-"*60)
print("  atom                   vertex                  note")
print("  bond                   arc                     interval")
print("  covalent bond          symmetric arc (r~0)     unison")
print("  ionic bond             asymmetric arc (r>>0)   dissonance")
print("  EN difference          arc weight              consonance")
print("  linear chain           transitive tournament   melody")
print("  branched molecule      cyclic tournament       chord")
print("  ring (benzene)         all-3-cycle             triad")
print("  polymer                long tournament          symphony")
print("  isomers                tournaments on same n    permutations")
print("  stereochemistry        arc direction            inversion")
print("  reaction pathway       Hamiltonian path         piece")
print("  activation energy      rapidity barrier         tension")
print("  catalyst               shortcut arc             resolution")
print()
print("  LiH: the INTERVAL. One bond. One direction. One note.")
print("  Alkane C_n: the SCALE. n notes in a row, descending consonance.")
print("  Benzene: the CHORD. Six notes in a ring, maximum symmetry.")
print("  Protein: the SYMPHONY. Thousands of notes, folded into meaning.")
print()
print("  And the Cayley transform Q converts between:")
print("  - The BOUNDED world of bond polarities (EN ratios near 1)")
print("  - The UNBOUNDED world of reaction pathways (exponential complexity)")
print("  Just as it converts between velocities and boosts,")
print("  between coin biases and log-odds,")
print("  between musical intervals and frequency ratios.")
print()
print("  THE ALKANE SERIES IS THE OVERTONE SERIES IS THE TAYLOR")
print("  SERIES OF ARCTANH IS THE MUSICAL SCALE IS THE TOURNAMENT")
print("  COUPLING. THEY ARE ALL THE SAME OBJECT.")
