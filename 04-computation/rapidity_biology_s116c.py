#!/usr/bin/env python3
"""rapidity_biology_s116c.py — Biology and evolution through the rapidity lens

Can rapidity illuminate biological scales, growth, and evolution?
"""
from math import sqrt, log, pi, exp, atanh, factorial

print("BIOLOGY AND EVOLUTION THROUGH RAPIDITY")
print("="*70)
print()

# ============================================================
print("="*70)
print("1. THE TREE OF LIFE IN RAPIDITY")
print("="*70)
print()
print("  Estimated number of species in major groups:")
groups = [
    ("Bacteria",           10000,    "prokaryotic"),
    ("Archaea",            500,      "prokaryotic"),
    ("Protists",           80000,    "eukaryotic"),
    ("Fungi",              144000,   "eukaryotic"),
    ("Plants",             390000,   "eukaryotic"),
    ("Insects",            1000000,  "eukaryotic"),
    ("Other invertebrates",200000,   "eukaryotic"),
    ("Fish",               35000,    "vertebrate"),
    ("Amphibians",         8000,     "vertebrate"),
    ("Reptiles",           11000,    "vertebrate"),
    ("Birds",              10000,    "vertebrate"),
    ("Mammals",            6400,     "vertebrate"),
    ("Total known",        8700000,  "all"),
]
print(f"  {'Group':<25s} {'Species':>10s}  {'rapidity':>10s}")
print("  " + "-"*50)
for name, count, typ in groups:
    rap = log(count)/2
    print(f"  {name:<25s} {count:>10,d}  {rap:10.4f}")

print()
print("  Insect diversity (rapidity 6.9) vs mammal diversity (rapidity 4.4):")
print(f"  Rapidity gap = {log(1000000)/2 - log(6400)/2:.4f}")
print(f"  = {(log(1000000)/2 - log(6400)/2)/(log(2)/2):.1f} octaves")
print("  Insects are about 7 octaves more diverse than mammals!")
print()

# ============================================================
print("="*70)
print("2. GENOME SIZE IN RAPIDITY")
print("="*70)
print()
print("  Genome sizes (base pairs):")
genomes = [
    ("Smallest virus (CircoV)", 1800),
    ("HIV",                     9700),
    ("E. coli",                 4600000),
    ("Yeast",                   12000000),
    ("Drosophila",              140000000),
    ("Human",                   3200000000),
    ("Wheat",                   16000000000),
    ("Paris japonica (record)", 150000000000),
]
print(f"  {'Organism':<30s} {'Base pairs':>15s}  {'rapidity':>10s}")
print("  " + "-"*60)
for name, bp in genomes:
    rap = log(bp)/2
    print(f"  {name:<30s} {bp:>15,d}  {rap:10.4f}")

print()
human_rap = log(3.2e9)/2
ecoli_rap = log(4.6e6)/2
print(f"  Human/E.coli genome ratio rapidity = {human_rap - ecoli_rap:.4f}")
print(f"  = {(human_rap - ecoli_rap)/(log(2)/2):.1f} octaves")
print(f"  Human genome is about {(human_rap - ecoli_rap)/(log(2)/2):.0f} octaves larger than E. coli's.")
print()

# ============================================================
print("="*70)
print("3. POPULATION SIZES AND RAPIDITY")
print("="*70)
print()
print("  World populations (approximate):")
populations = [
    ("Humans (2024)",           8000000000),
    ("Chickens",                34000000000),
    ("Cattle",                  1000000000),
    ("Ants (estimated)",        20000000000000000),
    ("Bacteria on Earth",       1e30),
    ("Stars in Milky Way",      200000000000),
    ("Neurons in human brain",  86000000000),
    ("Cells in human body",     37200000000000),
]
print(f"  {'Population':<30s} {'Count':>20s}  {'rapidity':>10s}")
print("  " + "-"*65)
for name, pop in populations:
    rap = log(pop)/2
    print(f"  {name:<30s} {pop:>20.0f}  {rap:10.4f}")

print()
print("  Bacteria vs humans: rapidity gap = {:.1f}".format(
    (log(1e30) - log(8e9))/2))
print("  = {:.0f} octaves".format((log(1e30) - log(8e9))/(2*log(2)/2)))
print()

# ============================================================
print("="*70)
print("4. EVOLUTIONARY TIME IN RAPIDITY")
print("="*70)
print()
print("  Key events in evolutionary history (years ago):")
events = [
    ("Earth formation",             4540000000),
    ("First life (prokaryotes)",    3800000000),
    ("Photosynthesis",              3400000000),
    ("Eukaryotes",                  2100000000),
    ("Multicellular life",          1200000000),
    ("Cambrian explosion",          541000000),
    ("First fish",                  530000000),
    ("Land plants",                 470000000),
    ("Insects",                     400000000),
    ("First tetrapods",             375000000),
    ("Dinosaurs",                   233000000),
    ("Mammals",                     225000000),
    ("Flowers",                     130000000),
    ("K-T extinction",              66000000),
    ("Primates",                    55000000),
    ("Hominids",                    6000000),
    ("Homo sapiens",                300000),
    ("Agriculture",                 10000),
    ("Industrial revolution",       250),
    ("Digital age",                 50),
]
print(f"  {'Event':<30s} {'Years ago':>15s}  {'rapidity':>10s}")
print("  " + "-"*60)
for name, years in events:
    rap = log(years)/2
    print(f"  {name:<30s} {years:>15,d}  {rap:10.4f}")

print()
print("  Earth formation to first life: rapidity gap = {:.4f}".format(
    log(4.54e9)/2 - log(3.8e9)/2))
print("  = {:.4f} octaves".format((log(4.54e9)/2 - log(3.8e9)/2)/(log(2)/2)))
print()
print("  Cambrian explosion to humans: rapidity gap = {:.4f}".format(
    log(541e6)/2 - log(300000)/2))
print("  = {:.1f} octaves".format((log(541e6)/2 - log(300000)/2)/(log(2)/2)))
print()
print("  In rapidity space, evolution is APPROXIMATELY UNIFORM!")
print("  The events are spread roughly evenly from rapidity 2 to 11.")
print("  This is because evolution generates novelty at a roughly")
print("  EXPONENTIAL rate (constant progress per log-time), and")
print("  rapidity = log/2 linearizes exponential growth.")
print()

# ============================================================
print("="*70)
print("5. ALLOMETRIC SCALING AND RAPIDITY")
print("="*70)
print()
print("  Kleiber's law: metabolic rate B ~ M^{3/4} (M = body mass)")
print("  In rapidity: rapidity(B) = (3/4)*rapidity(M) + const")
print("  This is a LINEAR relationship in rapidity space!")
print()
print("  The 3/4 exponent means: for every 4 octaves of mass,")
print("  metabolic rate increases by 3 octaves.")
print()
print("  Heart rate ~ M^{-1/4}: rapidity(HR) = -(1/4)*rapidity(M) + const")
print("  Lifespan ~ M^{1/4}: rapidity(L) = (1/4)*rapidity(M) + const")
print()
print("  Total heartbeats in lifetime ~ HR * L ~ M^0 = CONSTANT!")
print("  rapidity(total_beats) = rapidity(HR) + rapidity(L) = 0 + const")
print("  About 1.5 billion heartbeats for all mammals.")
print()

animals = [
    ("Shrew",           0.002, 1000, 2),
    ("Mouse",           0.02,  600,  3),
    ("Rat",             0.3,   400,  4),
    ("Rabbit",          2,     200,  9),
    ("Cat",             4,     150,  15),
    ("Dog",             20,    100,  13),
    ("Human",           70,    70,   80),
    ("Horse",           500,   40,   30),
    ("Elephant",        5000,  30,   70),
    ("Blue whale",      150000, 10,  80),
]
print(f"  {'Animal':<16s} {'Mass(kg)':>10s} {'HR(bpm)':>8s} {'Life(yr)':>8s}  {'rap(M)':>8s} {'rap(HR)':>8s}")
print("  " + "-"*65)
for name, mass, hr, life in animals:
    rm = log(mass)/2
    rhr = log(hr)/2
    print(f"  {name:<16s} {mass:>10.1f} {hr:>8d} {life:>8d}  {rm:>8.4f} {rhr:>8.4f}")

print()
print("  Verify Kleiber: slope of rapidity(HR) vs rapidity(M) should be -1/4.")
# Quick linear regression
masses = [a[1] for a in animals]
hrs = [a[2] for a in animals]
rm = [log(m)/2 for m in masses]
rhr = [log(h)/2 for h in hrs]
n = len(rm)
sx = sum(rm)
sy = sum(rhr)
sxx = sum(x*x for x in rm)
sxy = sum(x*y for x, y in zip(rm, rhr))
slope = (n*sxy - sx*sy) / (n*sxx - sx*sx)
intercept = (sy - slope*sx) / n
print(f"  Linear regression: slope = {slope:.4f} (expected -0.25)")
print(f"  Intercept = {intercept:.4f}")
print()

# ============================================================
print("="*70)
print("6. FIBONACCI IN BIOLOGY AS RAPIDITY")
print("="*70)
print()
print("  Fibonacci numbers appear in phyllotaxis (leaf arrangement),")
print("  petal counts, spiral counts, etc.")
print("  F_n ~ phi^n/sqrt(5).")
print("  rapidity(F_n) ~ n*ln(phi)/2 = n * 0.2406")
print()
fib = [1, 1]
for _ in range(15):
    fib.append(fib[-1] + fib[-2])
print(f"  n    F_n    rapidity(F_n)   n*ln(phi)/2")
for k in range(1, 14):
    rap = log(fib[k])/2 if fib[k] > 0 else 0
    expected = k * log((1+sqrt(5))/2) / 2
    print(f"  {k:2d}   {fib[k]:5d}    {rap:8.4f}        {expected:8.4f}")

print()
print("  Fibonacci rapidity grows LINEARLY in n with slope ln(phi)/2 = 0.2406.")
print("  This is the rapidity of the golden ratio divided by 2.")
print("  Each Fibonacci step adds 0.24 rapidity units ~ 0.7 octaves.")
print()

# ============================================================
print("="*70)
print("7. DNA INFORMATION CONTENT IN RAPIDITY")
print("="*70)
print()
print("  DNA uses 4 bases (A, C, G, T).")
print("  Information per base = log2(4) = 2 bits.")
print("  A genome of N bases contains 2N bits = N log2(4) bits.")
print()
print("  The number of possible sequences of length N: 4^N.")
print("  rapidity(4^N) = N*ln(4)/2 = N*ln(2) = N * 0.6931")
print("  Each base adds ln(2) = one OCTAVE of rapidity to the sequence space!")
print()
print("  For the human genome (N = 3.2 billion):")
print(f"  rapidity of sequence space = {3.2e9 * log(2):.4e}")
print(f"  = {3.2e9:.0e} octaves!")
print()
print("  But actual information content is much less (most DNA is non-coding).")
print("  Functional human genome ~ 1.5% coding => ~ 48 million coding bases.")
print(f"  Functional sequence space rapidity = {48e6 * log(2):.4e}")
print(f"  = {48e6:.0e} octaves of functional rapidity")
print()

# ============================================================
print("="*70)
print("8. DOUBLING TIMES AND RAPIDITY")
print("="*70)
print()
print("  Cell division time = doubling time.")
print("  After n doublings: 2^n cells.")
print("  rapidity(2^n) = n*ln(2)/2 = n * 0.3466")
print("  Each division adds exactly one OCTAVE of rapidity.")
print()
print("  From 1 cell to a human (~37 trillion cells):")
divisions_needed = log(37e12) / log(2)
print(f"  Divisions needed: log2(37e12) = {divisions_needed:.1f}")
print(f"  Rapidity: {log(37e12)/2:.4f}")
print(f"  = {divisions_needed:.0f} octaves from zygote to adult")
print()
print("  From 1 bacterium to E. coli colony of 10^9:")
ecoli_div = log(1e9) / log(2)
print(f"  Divisions: log2(10^9) = {ecoli_div:.1f}")
print(f"  = {ecoli_div:.0f} octaves of bacterial growth")
print()

# ============================================================
print("="*70)
print("9. EXTINCTION RATES AND RAPIDITY")
print("="*70)
print()
print("  Background extinction rate: ~0.1 species per million species-years")
print("  Current rate: ~100-1000x background")
print()
print("  Mass extinction events (% species lost):")
extinctions = [
    ("End Ordovician",     85, 444),
    ("Late Devonian",      75, 372),
    ("End Permian",        96, 252),
    ("End Triassic",       80, 201),
    ("End Cretaceous",     76, 66),
    ("Current (est.)",     30, 0),
]
print(f"  {'Event':<22s} {'% lost':>6s} {'Mya':>5s}  {'survival':>10s} {'rap(survival)':>14s}")
print("  " + "-"*65)
for name, pct, mya in extinctions:
    survival = (100 - pct) / 100
    if survival > 0.001:
        rap = atanh(2*survival - 1)
    else:
        rap = -5.0
    print(f"  {name:<22s} {pct:>5d}% {mya:>5d}  {survival:>10.2f} {rap:>14.4f}")

print()
print("  The End Permian (96% loss) has survival rapidity {:.4f}".format(
    atanh(2*0.04-1)))
print("  This is DEEPLY NEGATIVE — near total extinction in rapidity terms.")
print()

# ============================================================
print("="*70)
print("10. SENSORY RANGES IN RAPIDITY")
print("="*70)
print()
print("  Each sense has a dynamic range between threshold and saturation.")
print()

senses = [
    ("Human hearing",       20, 20000, "Hz"),
    ("Human vision",        400, 700, "nm (as frequency)"),
    ("Dog hearing",         67, 45000, "Hz"),
    ("Bat echolocation",    20000, 200000, "Hz"),
    ("Elephant infrasound", 1, 12000, "Hz"),
]
print(f"  {'Sense':<22s} {'Min':>10s} {'Max':>10s}  {'Octaves':>8s}")
print("  " + "-"*55)
for name, lo, hi, unit in senses:
    octaves = log(hi/lo) / log(2)
    print(f"  {name:<22s} {lo:>10,d} {hi:>10,d}  {octaves:>8.2f}")

print()
print("  Human hearing spans ~10 octaves (20 Hz to 20 kHz).")
print("  Human vision spans ~0.8 octaves (as we computed).")
print("  Bats: ~3.3 octaves of echolocation.")
print()
print("  In rapidity: hearing range = {:.4f}, vision range = {:.4f}".format(
    log(20000/20)/2, log(700/400)/2))
print("  Hearing is {:.1f}x wider than vision in rapidity space.".format(
    log(20000/20) / log(700/400)))
print()

# ============================================================
print("="*70)
print("GRAND SYNTHESIS: RAPIDITY IN BIOLOGY")
print("="*70)
print()
print("  BIOLOGY                RAPIDITY IS                    KEY VALUE")
print("  " + "-"*65)
print("  Species diversity      ln(#species)/2                 6.9 (insects)")
print("  Genome size            ln(base pairs)/2               10.9 (human)")
print("  Evolution timeline     ln(years ago)/2                2-11 (uniform!)")
print("  Allometric scaling     LINEAR in rapidity             B ~ M^{3/4}")
print("  Fibonacci growth       n * ln(phi)/2 per step         0.24/step")
print("  Cell divisions         n * ln(2)/2 per doubling       0.35/doubling")
print("  DNA sequence space     N * ln(4)/2                    1 octave/base")
print("  Sensory range          ln(max/min)/2                  ~10 oct (hearing)")
print("  Extinction severity    arctanh(2*survival - 1)        deeply negative")
print()
print("  Rapidity linearizes biological scaling laws.")
print("  Evolution is approximately uniform in rapidity time.")
print("  Each cell division adds one octave to the organism.")
print("  Each DNA base adds one octave to sequence space.")
