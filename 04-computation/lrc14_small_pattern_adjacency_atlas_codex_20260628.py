#!/usr/bin/env python3
"""HYP-3226 scout: small-pattern adjacency atlas for the LRC14 frontier.

The goal is not to prove LRC14.  It is to keep a large set of tempting small
patterns honest by attaching each one to typed proof payloads:

  preserved coordinate, destroyed coordinate, repair sidecar, and risk.

The scout is deliberately repo-local.  It scans recent hypotheses, results,
reflections, coordination notes, and forum drafts for keyword support, then
builds a ranked motif ledger and a tournament over motif families.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
import itertools
import math
import re
from typing import Iterable


ROOT = Path(__file__).resolve().parents[1]
SELF_GENERATED_OUTPUT = (
    ROOT / "05-knowledge/results/lrc14_small_pattern_adjacency_atlas_codex_20260628.out"
)


PAYLOAD_WEIGHTS = {
    "AP_NORMAL": 7,
    "TOEPLITZ": 7,
    "COV_LAYER": 6,
    "ORDERED_TAIL": 6,
    "TRAP_BOUNDARY": 7,
    "HB_PERRON": 6,
    "GREEN_LORENTZIAN": 5,
    "ANALYTIC_EQ": 7,
    "QUOTIENT_LEGALITY": 5,
    "PGF_ROOT": 4,
    "EDGE_PACKET": 4,
    "CHEBYSHEV": 5,
    "P_ADIC": 3,
    "SELBERG": 3,
    "CIRCUIT": 3,
    "GEOMETRY": 3,
    "SIDE_CARRIER": 2,
}

RISK_PENALTY = {
    "direct": 0,
    "sidecar": 3,
    "analogy": 7,
    "raw": 12,
}


@dataclass(frozen=True)
class Motif:
    ident: str
    label: str
    family: str
    pattern: str
    payloads: tuple[str, ...]
    destroyed: str
    repair: str
    risk: str
    keywords: tuple[str, ...]


MOTIFS: tuple[Motif, ...] = (
    Motif("M001", "consecutive plus doubled AP", "equality-face",
          "two-row skyline / primitive uniqueness",
          ("AP_NORMAL", "TOEPLITZ", "COV_LAYER", "ORDERED_TAIL"),
          "primitive-vs-dilation distinction",
          "dilation equality sidecar",
          "direct",
          ("doubled AP", "Pareto skyline", "primitive normal form")),
    Motif("M002", "11 non-AP exchange traps plus AP", "trap-boundary",
          "11/12 local maxima ledger",
          ("TRAP_BOUNDARY", "TOEPLITZ", "COV_LAYER"),
          "bulk exchange monotonicity",
          "moment-cone trap discharge",
          "direct",
          ("11 non-AP", "exchange traps", "lambda_min")),
    Motif("M003", "HYP-3204 exchange-rate ratio", "coefficient-pricing",
          "q3 gain priced by q0+q6 loss",
          ("ORDERED_TAIL", "TRAP_BOUNDARY", "AP_NORMAL"),
          "raw central mass q3",
          "ordered-state bimodality sidecar",
          "direct",
          ("12882/17161", "q3", "q0+q6")),
    Motif("M004", "AP support projection", "root-locus",
          "polarized cyclotomic support, not raw norm",
          ("AP_NORMAL", "PGF_ROOT", "CHEBYSHEV"),
          "raw closest-to-uniform cyclotomic norm",
          "AP residual polarization",
          "direct",
          ("AP_support", "39766/540225", "cyclotomic")),
    Motif("M005", "Toeplitz lambda-min margin", "moment-cone",
          "Caratheodory-Toeplitz interior slack",
          ("TOEPLITZ", "TRAP_BOUNDARY", "AP_NORMAL"),
          "PSD yes/no scalar",
          "lambda-min margin / Schur sidecar",
          "direct",
          ("Toeplitz", "lambda_min", "Caratheodory")),
    Motif("M006", "D1/D2/D3 layer split", "covariance",
          "three cyclic-distance covariance faces",
          ("COV_LAYER", "GREEN_LORENTZIAN", "TRAP_BOUNDARY"),
          "total covariance only",
          "distance-layer profile",
          "direct",
          ("D1", "D2", "D3", "distance layer")),
    Motif("M007", "Perron all-ones mode", "spectral-mode",
          "uniform coherent mode vs AFM mode shift",
          ("HB_PERRON", "COV_LAYER", "GREEN_LORENTZIAN"),
          "raw Perron alignment as terminal maximum",
          "boundary-aware covariance matrix sidecar",
          "sidecar",
          ("Perron", "all-ones", "AFM")),
    Motif("M008", "Hermite-Biehler E/O legs", "stability-glue",
          "E=x^2+5x+4, O=x^2+4x+1",
          ("HB_PERRON", "CHEBYSHEV", "PGF_ROOT"),
          "full miss-PGF self-inversive defect",
          "Joukowski off-circle sidecar",
          "direct",
          ("Hermite-Biehler", "Wronskian", "E=x^2+5x+4")),
    Motif("M009", "Chebyshev V7 double root", "magic-function",
          "V7(u)-2 has squared de Moivre cubic",
          ("CHEBYSHEV", "AP_NORMAL", "ANALYTIC_EQ"),
          "raw cyclotomic numerology",
          "Delsarte/Cohn-Elkies dual certificate",
          "sidecar",
          ("V_7", "Chebyshev", "Delsarte")),
    Motif("M010", "Lee-Yang apex zero at 1/p=7", "analytic-obstruction",
          "config-blind certificates fail at apex 7",
          ("ANALYTIC_EQ", "PGF_ROOT", "SIDE_CARRIER"),
          "algebraic closure without equidistribution",
          "effective equidistribution / Weyl sidecar",
          "direct",
          ("Lee-Yang", "apex", "1/p")),
    Motif("M011", "K3 flip kernel", "edge-packet",
          "two-state edge flip kernel with minority edge debt",
          ("EDGE_PACKET", "QUOTIENT_LEGALITY", "PGF_ROOT"),
          "state-level edge roles",
          "minority-edge / Worpitzky sidecar",
          "direct",
          ("K3", "edge flip", "minority edge")),
    Motif("M012", "Worpitzky rows 1,4,1 and 1,11,11,1", "worpitzky",
          "Eulerian ascent payload",
          ("EDGE_PACKET", "HB_PERRON", "QUOTIENT_LEGALITY"),
          "ordered descent word",
          "Eulerian / odd sidecar",
          "direct",
          ("Worpitzky", "1,4,1", "1,11,11,1")),
    Motif("M013", "pair functions a+b,a*b,a^b,b^a", "function-quotient",
          "two unordered-safe, two ordered-sensitive",
          ("QUOTIENT_LEGALITY", "EDGE_PACKET"),
          "ordered channel",
          "ordered-pair sidecar",
          "direct",
          ("a+b", "a*b", "a^b", "b^a")),
    Motif("M014", "n=4 OR compression", "finite-chart",
          "x=a OR c, y=b OR c",
          ("QUOTIENT_LEGALITY", "EDGE_PACKET", "SIDE_CARRIER"),
          "S fiber PGF and canary coordinate",
          "filler/canary/deletion sidecar",
          "direct",
          ("x=a OR c", "canary", "S fiber")),
    Motif("M015", "degree-4 resolvent ceiling", "resolvent",
          "u^4-5u^2+4 below quintic wall",
          ("HB_PERRON", "QUOTIENT_LEGALITY", "CHEBYSHEV"),
          "odd coordinate debt",
          "biquadratic plus Worpitzky sidecar",
          "sidecar",
          ("u^4-5u^2+4", "degree 4", "resolvent")),
    Motif("M016", "duodecimal trap/equality count", "duodecimal",
          "12 local endpoints / 12 source witnesses",
          ("TRAP_BOUNDARY", "EDGE_PACKET", "SIDE_CARRIER"),
          "why 12 appears in different charts",
          "chart id and fiber sidecar",
          "sidecar",
          ("duodecimal", "12 local", "local maxima")),
    Motif("M017", "doubled triangular numbers", "sequence-shape",
          "2,6,12,20,30,42 = n(n+1)",
          ("GEOMETRY", "SIDE_CARRIER"),
          "raw sequence coincidence",
          "partial-cube/simplex bridge-rank sidecar",
          "analogy",
          ("2,6,12,20,30,42", "triangular", "partial cube")),
    Motif("M018", "A000568 perspective break", "perspective-functor",
          "node/edge perspective count defect",
          ("QUOTIENT_LEGALITY", "EDGE_PACKET", "SIDE_CARRIER"),
          "unrooted orientation and cross-sector data",
          "edge-envelope / sector-deck sidecar",
          "sidecar",
          ("A000568", "perspective", "edge-envelope")),
    Motif("M019", "Green-current bottleneck", "network",
          "covariance as response matrix",
          ("GREEN_LORENTZIAN", "COV_LAYER", "TRAP_BOUNDARY"),
          "current location",
          "effective-resistance profile",
          "sidecar",
          ("Green-current", "effective resistance", "Schur")),
    Motif("M020", "Lorentzian exchange chamber", "valuated-exchange",
          "co-emptiness finite differences as exchange slacks",
          ("GREEN_LORENTZIAN", "TRAP_BOUNDARY", "AP_NORMAL"),
          "orbit/mirror/two-block sidecars",
          "valuated matroid / tropical Plucker sidecar",
          "sidecar",
          ("Lorentzian", "valuated", "Plucker")),
    Motif("M021", "Selberg trace regularization", "spectral-regularizer",
          "regularize spectral error terms",
          ("SELBERG", "ANALYTIC_EQ", "SIDE_CARRIER"),
          "noncompact spectral leakage",
          "principal-part / trace sidecar",
          "analogy",
          ("Selberg", "trace regularization", "Lindelof")),
    Motif("M022", "p-adic tau valuation stability", "p-adic",
          "valuation sidecar for modular arithmetic",
          ("P_ADIC", "ANALYTIC_EQ", "SIDE_CARRIER"),
          "motivic speculation without LRC payload",
          "Hensel/local root stability sidecar",
          "analogy",
          ("tau", "p-adic", "Hensel")),
    Motif("M023", "Skewes sign-change warning", "prime-oscillation",
          "first crossing as extreme oscillation marker",
          ("ANALYTIC_EQ", "SIDE_CARRIER"),
          "prime-counting scalar imported raw",
          "discrepancy/zero-free sidecar",
          "analogy",
          ("Skewes", "prime-counting", "li(x)")),
    Motif("M024", "Helfgott-Ruzsa additive constant", "additive-comb",
          "small sumset as compression warning",
          ("QUOTIENT_LEGALITY", "ANALYTIC_EQ", "SIDE_CARRIER"),
          "sumset analogy without packet map",
          "Freiman/Ruzsa model sidecar",
          "analogy",
          ("Helfgott", "Ruzsa", "sumset")),
    Motif("M025", "Collatz two-block determinant", "two-block-adic",
          "2^E - 3^k mirrors LRC two-block residual",
          ("P_ADIC", "QUOTIENT_LEGALITY", "ANALYTIC_EQ"),
          "density-blind residual set",
          "Baker/linear-forms sidecar",
          "sidecar",
          ("Collatz", "two-block", "2-adic")),
    Motif("M026", "Roth-Vaughan discrepancy sidecar", "discrepancy",
          "higher-dimensional discrepancy as analytic repair",
          ("ANALYTIC_EQ", "SIDE_CARRIER"),
          "using discrepancy without endpoint packet",
          "height/lattice sidecar",
          "sidecar",
          ("Roth-Vaughan", "discrepancy", "height")),
    Motif("M027", "BDH variance packet", "prime-residue",
          "mean-square residue regularity",
          ("ANALYTIC_EQ", "SELBERG", "SIDE_CARRIER"),
          "prime residue theorem as raw import",
          "conductor/residue-family sidecar",
          "analogy",
          ("Barban-Davenport-Halberstam", "residues", "variance")),
    Motif("M028", "Hensel-Newton 2-adic edge case", "p-adic-lift",
          "local singular lifting guardrail",
          ("P_ADIC", "QUOTIENT_LEGALITY", "SIDE_CARRIER"),
          "assuming lift without singular residue audit",
          "Hensel-Newton/Krasner sidecar",
          "sidecar",
          ("Hensel-Newton", "2-adic", "lifting")),
    Motif("M029", "Fejer-Riesz square", "positive-trig",
          "PSD Toeplitz slack factorization",
          ("TOEPLITZ", "AP_NORMAL", "CHEBYSHEV"),
          "factor existence without endpoint slack",
          "explicit square factor sidecar",
          "sidecar",
          ("Fejer", "Riesz", "square")),
    Motif("M030", "Verblunsky/Schur parameters", "opuc",
          "moment-cone boundary coordinates",
          ("TOEPLITZ", "TRAP_BOUNDARY", "SIDE_CARRIER"),
          "PSD margin compressed to one eigenvalue",
          "Schur parameter word",
          "sidecar",
          ("Verblunsky", "Schur", "OPUC")),
    Motif("M031", "random-current coupling order", "ising-current",
          "covariance comparison via currents",
          ("COV_LAYER", "GREEN_LORENTZIAN", "ANALYTIC_EQ"),
          "plain positive association",
          "current switching sidecar",
          "sidecar",
          ("random-current", "Aizenman", "GKS")),
    Motif("M032", "Asano contraction warning", "lee-yang",
          "zero-free polydisk contraction",
          ("PGF_ROOT", "ANALYTIC_EQ", "SIDE_CARRIER"),
          "config-blind zero-free proof",
          "apex-7 obstruction sidecar",
          "sidecar",
          ("Asano", "zero-free", "polydisk")),
    Motif("M033", "Beurling-Selberg minorant floor", "minorant",
          "analytic lower envelope for wide arcs",
          ("ANALYTIC_EQ", "TOEPLITZ", "SIDE_CARRIER"),
          "minorant scalar without leakage audit",
          "Gaussian leakage / endpoint sidecar",
          "sidecar",
          ("Beurling-Selberg", "minorant", "Gaussian")),
    Motif("M034", "phi4 kappa4 stabilizer", "field-shape",
          "quartic stabilizer as sidecar, not minimum",
          ("SIDE_CARRIER", "COV_LAYER"),
          "raw kappa4 minimization",
          "stabilizer-only label",
          "sidecar",
          ("phi4", "kappa4", "quartic")),
    Motif("M035", "2002 Pascal/binomial currency", "pascal",
          "C(14,5)=2*7*11*13",
          ("SIDE_CARRIER", "ANALYTIC_EQ"),
          "entropy proof handle",
          "Pascal pair-mass sidecar",
          "analogy",
          ("2002", "Pascal", "binomial")),
    Motif("M036", "1/7 scalar refutation", "apex-prime",
          "near 1/7 but no exact law",
          ("ANALYTIC_EQ", "SIDE_CARRIER"),
          "universal apex-prime scalar",
          "named associator residual",
          "direct",
          ("1/7", "associator", "ratio")),
    Motif("M037", "odd/even ratio 3.149364", "odd-dominance",
          "odd Worpitzky term dominates even fold",
          ("HB_PERRON", "EDGE_PACKET", "SIDE_CARRIER"),
          "even-only proof",
          "odd associator sidecar",
          "direct",
          ("3.149364", "odd", "even")),
    Motif("M038", "n=5 score-compression failure", "compression-wall",
          "score->iso first fails at n=5",
          ("QUOTIENT_LEGALITY", "EDGE_PACKET"),
          "score sequence scalar",
          "iso/fiber sidecar",
          "direct",
          ("n=5", "score", "compression")),
    Motif("M039", "K4 Einheit exact section", "finite-basis",
          "two-free-arc matching section",
          ("QUOTIENT_LEGALITY", "EDGE_PACKET", "SIDE_CARRIER"),
          "fixed-path cover fiber collision",
          "exact subbasis sidecar",
          "direct",
          ("Einheit", "Klein", "two-free")),
    Motif("M040", "S fiber multiplicity five", "fiber-warning",
          "T,+,- singletons but S has five fibers",
          ("QUOTIENT_LEGALITY", "EDGE_PACKET"),
          "class count compression",
          "fiber PGF sidecar",
          "direct",
          ("S=5", "fiber", "fixed path")),
    Motif("M041", "edge-perspective tip-tail packet", "edge-perspective",
          "directed edge as observer with tip/tail children",
          ("EDGE_PACKET", "QUOTIENT_LEGALITY", "SIDE_CARRIER"),
          "node-only perspective",
          "tip-tail sector word",
          "sidecar",
          ("tip-tail", "edge perspective", "tail")),
    Motif("M042", "observer-cut payload", "observer",
          "cut payload with named forgotten coordinate",
          ("QUOTIENT_LEGALITY", "SIDE_CARRIER", "CIRCUIT"),
          "unnamed observer loss",
          "observer gluing sidecar",
          "sidecar",
          ("observer-cut", "payload", "gluing")),
    Motif("M043", "normal-fan Cech barcode", "component-packet",
          "component/barcode packet for chart overlap",
          ("AP_NORMAL", "SIDE_CARRIER", "GEOMETRY"),
          "component scalar",
          "Cech/barcode certificate",
          "sidecar",
          ("normal-fan", "Cech", "barcode")),
    Motif("M044", "Desargues median route center", "median",
          "route compatibility via median center",
          ("GEOMETRY", "QUOTIENT_LEGALITY", "SIDE_CARRIER"),
          "non-median route compatibility",
          "route-state closure sidecar",
          "analogy",
          ("Desargues", "median", "route center")),
    Motif("M045", "Minkowski q-body threshold", "geometry-number",
          "lattice/height fence for Diophantine estimates",
          ("GEOMETRY", "ANALYTIC_EQ", "SIDE_CARRIER"),
          "lattice analogy without height packet",
          "Roth-Minkowski sidecar",
          "analogy",
          ("Minkowski", "Roth", "height")),
    Motif("M046", "circuit complexity ledger", "proof-circuit",
          "proof-DAG missing input vector",
          ("CIRCUIT", "QUOTIENT_LEGALITY", "SIDE_CARRIER"),
          "black-box proof step",
          "proof-circuit gate sidecar",
          "sidecar",
          ("circuit", "proof-DAG", "missing-input")),
    Motif("M047", "Fibonacci/fibbinary partial cube", "partial-cube",
          "no-adjacent support as cube subgraph",
          ("GEOMETRY", "QUOTIENT_LEGALITY", "SIDE_CARRIER"),
          "using fibbinary as numerology",
          "partial-cube bridge-rank sidecar",
          "analogy",
          ("fibbinary", "partial cube", "Zeckendorf")),
    Motif("M048", "Moser-de Bruijn sparse basis", "sparse-basis",
          "unique base-4 digit support",
          ("QUOTIENT_LEGALITY", "SIDE_CARRIER"),
          "basis analogy without target function",
          "basis uniqueness sidecar",
          "analogy",
          ("Moser", "de Bruijn", "sparse")),
    Motif("M049", "Markov numbers / Hurwitz orbit", "markov-hurwitz",
          "finite-address seeds for modular orbit",
          ("GEOMETRY", "ANALYTIC_EQ", "SIDE_CARRIER"),
          "Markov numerology",
          "Hurwitz finite-address seed",
          "analogy",
          ("Markov", "Hurwitz", "finite-address")),
    Motif("M050", "cannonball/Pell recurrence", "pell",
          "quadratic recurrence warning",
          ("GEOMETRY", "P_ADIC", "SIDE_CARRIER"),
          "Pell analogy without residual equation",
          "quadratic-form sidecar",
          "analogy",
          ("Pell", "cannonball", "square pyramidal")),
    Motif("M051", "square-peg four-witness gate", "four-witness",
          "rectangle/hourglass residue packet",
          ("TOEPLITZ", "GEOMETRY", "SIDE_CARRIER"),
          "four-point analogy without scale",
          "Toeplitz square-peg scale gate",
          "sidecar",
          ("square-peg", "rectangle", "hourglass")),
    Motif("M052", "S217 fixed-path diagonal flow", "diagonal-flow",
          "rectangle/hourglass cycle residues",
          ("EDGE_PACKET", "GEOMETRY", "SIDE_CARRIER"),
          "cycle residue without path id",
          "fixed-path diagonal-flow sidecar",
          "sidecar",
          ("S217", "diagonal flow", "hourglass")),
    Motif("M053", "A000568 n<=7 tameness", "tameness-window",
          "bounded-core below Abel-Ruffini wall",
          ("QUOTIENT_LEGALITY", "EDGE_PACKET", "HB_PERRON"),
          "class-count scalar",
          "perspective-depth sidecar",
          "sidecar",
          ("n<=7", "A000568", "Abel-Ruffini")),
    Motif("M054", "series-parallel nested ear", "ear-decomposition",
          "nested/odd/directed ear grammar",
          ("EDGE_PACKET", "TRAP_BOUNDARY", "SIDE_CARRIER"),
          "ear fact as graph folklore",
          "ear payload id",
          "sidecar",
          ("ear decomposition", "nested ear", "odd ear")),
    Motif("M055", "factor-critical odd ear", "matching",
          "odd ear as parity repair",
          ("EDGE_PACKET", "SIDE_CARRIER"),
          "matching analogy without odd payload",
          "parity/odd-ear sidecar",
          "analogy",
          ("factor-critical", "odd ear", "matching")),
    Motif("M056", "Tutte/Potts deletion-contraction", "dc-recursion",
          "tie induction as deletion-contraction",
          ("CIRCUIT", "EDGE_PACKET", "PGF_ROOT"),
          "recursion without target partition function",
          "Potts/Tutte packet sidecar",
          "sidecar",
          ("Tutte", "Potts", "deletion-contraction")),
    Motif("M057", "relation lattice Poisson sum", "lattice-fourier",
          "lonely time as lattice-resonance sum",
          ("ANALYTIC_EQ", "PGF_ROOT", "GEOMETRY"),
          "measure-only proof",
          "relation-lattice short-vector sidecar",
          "sidecar",
          ("relation lattice", "Poisson", "short-vector")),
    Motif("M058", "permanent/determinant sign symmetry", "sign-symmetry",
          "need sign-reversing symmetry to solve cancellation",
          ("CIRCUIT", "ANALYTIC_EQ", "SIDE_CARRIER"),
          "unsigned permanent scalar",
          "sign-involution sidecar",
          "analogy",
          ("permanent", "determinant", "sign-reversing")),
    Motif("M059", "Paley/QR eigenvalue packet", "paley",
          "quadratic-residue spectral shadow",
          ("EDGE_PACKET", "PGF_ROOT", "SIDE_CARRIER"),
          "Paley analogy without quotient map",
          "Gauss-sum/eigenvalue sidecar",
          "analogy",
          ("Paley", "quadratic residue", "Gauss")),
    Motif("M060", "Beraha/Mahler/subshift hint", "symbolic-dynamics",
          "root loci as subshift/Mahler boundary",
          ("PGF_ROOT", "CHEBYSHEV", "SIDE_CARRIER"),
          "symbolic-dynamics name without root packet",
          "root-locus packet sidecar",
          "analogy",
          ("Beraha", "Mahler", "subshift")),
    Motif("M061", "Stark L-value / Gamma0(7)", "modular-magic",
          "possible level-7 magic-function arithmetic",
          ("CHEBYSHEV", "SELBERG", "ANALYTIC_EQ"),
          "modular form analogy without principal part",
          "modular cusp principal-part gate",
          "analogy",
          ("Stark", "Gamma_0(7)", "L-value")),
    Motif("M062", "Lehmer tau nonvanishing", "modular-nonvanishing",
          "tau valuation as nonvanishing-density sidecar",
          ("P_ADIC", "SELBERG", "SIDE_CARRIER"),
          "Lehmer conjecture as raw import",
          "valuation stability sidecar",
          "raw",
          ("Lehmer", "tau", "non-vanishing")),
    Motif("M063", "Lindelof horizontal zeros", "lindelof",
          "zero-distribution regularization guardrail",
          ("SELBERG", "ANALYTIC_EQ", "SIDE_CARRIER"),
          "assuming LH as proof input",
          "explicit hypothesis boundary sidecar",
          "raw",
          ("Lindelof", "zeros", "zeta")),
    Motif("M064", "HYP-2990 no-free-slider", "meta-rule",
          "forget only with zero defect or repair sidecar",
          ("QUOTIENT_LEGALITY", "SIDE_CARRIER", "CIRCUIT"),
          "all compressed scalar routes",
          "residual defect meter",
          "direct",
          ("No-Free-Slider", "HYP-2990", "sidecar")),
    Motif("M065", "controlled-forgetting ladder", "meta-rule",
          "k-depth observers require retained sidecars",
          ("QUOTIENT_LEGALITY", "EDGE_PACKET", "SIDE_CARRIER"),
          "observer depth collapse",
          "observer-extension payload",
          "direct",
          ("controlled-forgetting", "observer", "payload")),
    Motif("M066", "certificate-Helly separation", "convex-certificate",
          "small dictionary intersection forces AP",
          ("AP_NORMAL", "TOEPLITZ", "TRAP_BOUNDARY"),
          "coordinate-by-coordinate maxima only",
          "Helly/Farkas sidecar",
          "direct",
          ("Helly", "Farkas", "dictionary")),
    Motif("M067", "normal-cone dual slack", "dual-certificate",
          "one dual controls many visible faces",
          ("AP_NORMAL", "TOEPLITZ", "ORDERED_TAIL", "COV_LAYER"),
          "another isolated scalar",
          "shared slack decomposition",
          "direct",
          ("normal cone", "dual", "slack")),
    Motif("M068", "multi-chart proof split", "chart-proof",
          "bulk exchange plus boundary moment cone plus odd gluing",
          ("TRAP_BOUNDARY", "TOEPLITZ", "HB_PERRON", "ORDERED_TAIL"),
          "single global move theorem",
          "chart-overlap certificate",
          "direct",
          ("multi-chart", "bulk", "boundary")),
    Motif("M069", "raw famous-problem magnet", "risk",
          "Skewes/tau/LH/Collatz as name-only imports",
          ("SIDE_CARRIER",),
          "all LRC payload",
          "reject or type as sidecar",
          "raw",
          ("Skewes", "Lindelof", "Collatz")),
    Motif("M070", "forum base64 draft warning", "coordination",
          "forum posts may carry motifs but require decoding and typing",
          ("CIRCUIT", "QUOTIENT_LEGALITY"),
          "opaque payload transport",
          "decode/source-id sidecar",
          "sidecar",
          ("base64", "forum", "post_178")),
    Motif("M071", "PFR additive-model integrity", "pfr",
          "additive resonance must emit a Freiman model sidecar",
          ("QUOTIENT_LEGALITY", "ANALYTIC_EQ", "SIDE_CARRIER"),
          "additive-combinatorics name without LRC packet",
          "Freiman-model endpoint / normal-fan sidecar",
          "analogy",
          ("PFR", "Polynomial Freiman-Ruzsa", "Freiman-Ruzsa")),
    Motif("M072", "conductance/Fiedler trap graph", "conductance-graph",
          "trap-discharge graph with algebraic-connectivity defect island",
          ("GREEN_LORENTZIAN", "TRAP_BOUNDARY", "CIRCUIT", "SIDE_CARRIER"),
          "scalar conductance ranking",
          "Fiedler-cut / Schur-complement sidecar",
          "direct",
          ("conductance graph", "Fiedler", "M-matrix", "Schur-complement")),
    Motif("M073", "comb-overlap Gram kernel", "magic-gram",
          "K(p,q)=meas(D_p cap D_q) as PSD/Bochner kernel",
          ("TOEPLITZ", "CHEBYSHEV", "AP_NORMAL", "ANALYTIC_EQ", "COV_LAYER"),
          "spatial minorant LP scalar",
          "Gram/Bochner autocorrelation sidecar",
          "direct",
          ("comb-overlap", "Gram kernel", "Bochner", "K(p,q)")),
    Motif("M074", "single-arc peeling recursion", "single-arc",
          "speed-1 tooth gives meas(intersection)=1/(7 max S)",
          ("AP_NORMAL", "ANALYTIC_EQ", "QUOTIENT_LEGALITY", "SIDE_CARRIER"),
          "assuming all teeth behave like the speed-1 tooth",
          "peeling-recursion endpoint sidecar",
          "direct",
          ("single-arc", "peeling recursion", "1/(7 max", "D_1")),
    Motif("M075", "order-3 overlap correction constants", "triple-overlap",
          "binding rows leave clean third-order rational residues",
          ("ANALYTIC_EQ", "PGF_ROOT", "SIDE_CARRIER"),
          "order-2 Gram kernel as complete proof",
          "triple-overlap / associator residue sidecar",
          "sidecar",
          ("-37/1092", "-61/588", "triple-overlap", "order-3")),
    Motif("M076", "induction-base critical flag", "literature-base",
          "published frontier may be n<=10 rather than n<=13",
          ("ANALYTIC_EQ", "QUOTIENT_LEGALITY", "SIDE_CARRIER"),
          "conditional base hidden inside a lift",
          "verified-base or n=11..13 sidecar",
          "direct",
          ("induction base", "n<=10", "n=11,12,13", "Rosenfeld")),
    Motif("M077", "Chen-Cusick 23-to-14 lift", "modulus-bridge",
          "baseline 1/23 is a floor to lift toward the 1/14 apex target",
          ("ANALYTIC_EQ", "GEOMETRY", "SIDE_CARRIER"),
          "starting the proof from scratch",
          "mod-23 floor sidecar; 23/M=2/23 is bounded-bank coincidence only",
          "sidecar",
          ("Chen-Cusick", "1/23", "1/14", "bounded-D coincidence")),
    Motif("M078", "polyhedron-zonotope flatness route", "geometry-of-numbers",
          "LRC as lattice point / covering-radius flatness certificate",
          ("GEOMETRY", "ANALYTIC_EQ", "QUOTIENT_LEGALITY", "SIDE_CARRIER"),
          "volume or polyhedron name without certificate map",
          "flatness / covering-radius sidecar",
          "sidecar",
          ("polyhedron", "flatness", "zonotope", "covering radius")),
    Motif("M079", "Rosenfeld exponential-sum Node-3 route", "far-speed",
          "large-speed branch controlled by modern exponential-sum estimates",
          ("ANALYTIC_EQ", "SELBERG", "SIDE_CARRIER"),
          "bounded-core certificate reused on far speeds",
          "Node-3 exponential-sum sidecar",
          "sidecar",
          ("Rosenfeld", "exponential-sum", "Node-3", "Tao")),
    Motif("M080", "shell L_y magic quartic", "finite-shell-magic",
          "f(n)=((n-1)(n-2)(n-4)(n-5))/4 gives 10q0+q3+10q6",
          ("TOEPLITZ", "CHEBYSHEV", "AP_NORMAL", "ORDERED_TAIL", "ANALYTIC_EQ"),
          "cyclic PSD positivity as false terminal certificate",
          "finite Delsarte/Newton shell-dual sidecar",
          "direct",
          ("10q0 + q3 + 10q6", "shell values", "magic_deficit", "rho >= 18.019")),
    Motif("M081", "Gamma0(7) Eisenstein coefficient engine", "modular-coefficients",
          "E_7=(7E_2(7tau)-E_2(tau))/6 supplies level-7 divisor fibers",
          ("CHEBYSHEV", "SELBERG", "ANALYTIC_EQ", "SIDE_CARRIER"),
          "raw modular-form name without finite LP rows",
          "match q-coefficients to comb-overlap Gram and Toeplitz rows",
          "sidecar",
          ("Gamma0(7)", "E_7(tau)", "sigma_1", "7 E_2(7 tau)")),
    Motif("M082", "Beraha/Mahler height gauge", "beraha-height",
          "B7 cubic and Mahler(m)=B7-1 monitor de Moivre perturbation height",
          ("CHEBYSHEV", "PGF_ROOT", "SIDE_CARRIER"),
          "height constant mistaken for inequality",
          "height-to-slack comparison sidecar",
          "sidecar",
          ("Beraha", "B^3 - 5B^2", "Mahler(m)", "B_7 - 1")),
    Motif("M083", "subshift transfer Perron defect", "transfer-operator",
          "AP autocorrelation 7-|n| is rank-one transfer state",
          ("PGF_ROOT", "HB_PERRON", "ANALYTIC_EQ", "SIDE_CARRIER"),
          "symbolic dynamics analogy without retained predicate",
          "transfer-matrix slack / Perron-defect sidecar",
          "sidecar",
          ("subshift", "transfer-operator", "7-|n|", "rank-one transfer")),
    Motif("M084", "Dirichlet-L/Stark denominator guardrail", "l-value-guardrail",
          "mod-7 L(-1,chi) denominators contract to 7 despite discriminant 7^2",
          ("SELBERG", "ANALYTIC_EQ", "P_ADIC", "SIDE_CARRIER"),
          "conductor-49 intuition promoted to cap formula",
          "explicit L-value denominator audit sidecar",
          "sidecar",
          ("Dirichlet-L", "Stark", "denominator 7", "conductor-49")),
    Motif("M085", "three-gap Stern-Brocot cap-kernel recursion", "cap-kernel-recursion",
          "K(a,b)=g(a,b)/(7ab) with g piecewise-linear by continued fractions",
          ("TOEPLITZ", "CHEBYSHEV", "ANALYTIC_EQ", "GEOMETRY", "SIDE_CARRIER"),
          "pairwise kernel mistaken for full cap proof",
          "order-3 overlap / inclusion-exclusion sidecar",
          "direct",
          ("K(a,b)=g(a,b)/(7ab)", "Stern-Brocot", "three-gap", "continued-fraction")),
    Motif("M086", "scale-normal recursion ledger", "scale-normal-recursion",
          "normalize scale, expose first surviving coordinate, attach sidecar, recurse",
          ("QUOTIENT_LEGALITY", "ANALYTIC_EQ", "SIDE_CARRIER", "CIRCUIT", "GEOMETRY"),
          "scale quotient without emitted cocycle",
          "scale-fiber exactness / omega_Q sidecar",
          "sidecar",
          ("scale-normal", "renormalization_depth", "omega_Q", "primitive projective shape")),
    Motif("M087", "LRC(2p) moment-order ladder", "moment-order-ladder",
          "apex depth (p+1)/2 and ladder depth (p-1)/2 match cyclotomic degree",
          ("CHEBYSHEV", "ANALYTIC_EQ", "SIDE_CARRIER"),
          "family-depth law used as local certificate",
          "finite p=7 certificate plus base-flag sidecar",
          "sidecar",
          ("moment-order DEPTH", "apex depth", "(p+1)/2", "cap_k = cap_{k-1}")),
    Motif("M088", "2-adic reflection fold", "reflection-fold",
          "s -> 6-s complement halves the k=8 degree-4 packet to a quadratic fold",
          ("HB_PERRON", "CHEBYSHEV", "QUOTIENT_LEGALITY", "SIDE_CARRIER"),
          "degree halving without odd-coordinate debt",
          "even/odd Worpitzky-HB gluing sidecar",
          "sidecar",
          ("2-adic reflection fold", "degree halving", "s -> 6-s", "u^4-5u^2+4")),
    Motif("M089", "modulus-covariance apex break", "modulus-covariance",
          "K^(2n)/K^(n)=1/2 until the apex n/2 break in the antipode half",
          ("TOEPLITZ", "ANALYTIC_EQ", "GEOMETRY", "TRAP_BOUNDARY", "SIDE_CARRIER"),
          "clean scale law used past its apex range",
          "antipode-half deviation / order-3 correction sidecar",
          "direct",
          ("modulus-covariance", "K^(2n)/K^(n)", "apex break", "antipode half")),
    Motif("M090", "cyclotomic subfield mode lattice", "subfield-mode-lattice",
          "Mobius/Eisenstein/Legendre/cubic/sextic modes form the Q(zeta7) tower",
          ("CHEBYSHEV", "SELBERG", "ANALYTIC_EQ", "P_ADIC", "SIDE_CARRIER"),
          "character-mode analogy without signed recursion packet",
          "chi_3 cubic-mode / L-value ladder sidecar",
          "sidecar",
          ("subfield lattice", "chi_3", "cubic de Moivre", "Gaussian periods")),
    Motif("M091", "cyclotomic factor grading", "cyclotomic-factor-grading",
          "mode=(x-1)^depth*Phi_d separates moment depth from character factor",
          ("CHEBYSHEV", "ANALYTIC_EQ", "QUOTIENT_LEGALITY", "SIDE_CARRIER"),
          "recurrence formula used without its character factor",
          "Phi_d character / moment-depth grading sidecar",
          "direct",
          ("cyclotomic factors", "(x-1)^depth", "Phi_d", "character factor")),
    Motif("M092", "signed address chart-change sheaf", "signed-chart-sheaf",
          "A..G recurrences are local chart addresses, not global letters",
          ("QUOTIENT_LEGALITY", "EDGE_PACKET", "CIRCUIT", "SIDE_CARRIER"),
          "same sign word reused after forgetting local slots",
          "signed_address_chart / chart_change_map sidecar",
          "direct",
          ("signed address", "chart-change", "A+B+D-C-E-F+G", "local address slots")),
    Motif("M093", "AP self-dual Fejer equidistribution certificate", "fejer-equidistribution",
          "AP autocorrelation is both maximal coherence and Fejer/Vaaler rescue",
          ("ANALYTIC_EQ", "CHEBYSHEV", "AP_NORMAL", "PGF_ROOT", "SIDE_CARRIER"),
          "coherence-vs-arithmetic split treated as two unrelated objects",
          "signed Fejer/Vaaler tail with phi(p) reserve sidecar",
          "direct",
          ("coherence-vs-arithmetic", "Gauss sum", "phi(p)", "Fejer/Vaaler")),
    Motif("M094", "totally-real cap field conductor packet", "totally-real-cap-field",
          "cap/dip live in Q(cos2pi/7) and binding rows carry conductor 7^1,2",
          ("CHEBYSHEV", "P_ADIC", "ANALYTIC_EQ", "TOEPLITZ", "SIDE_CARRIER"),
          "field-name import without binding-row denominator control",
          "cyclotomic-unit / totally-positive trace sidecar",
          "direct",
          ("Q(cos2pi/7)", "disc=49", "7^2-conductor", "totally-positive square")),
)


def candidate_files() -> list[Path]:
    patterns = [
        "00-navigation/*.md",
        "05-knowledge/hypotheses/HYP-31*.md",
        "05-knowledge/hypotheses/HYP-32*.md",
        "05-knowledge/results/*.out",
        "07-reflections/*.md",
        "comms/*.md",
        "poke-forum/post_178*.md",
    ]
    files: list[Path] = []
    for pattern in patterns:
        for path in ROOT.glob(pattern):
            if path == SELF_GENERATED_OUTPUT:
                continue
            if path.is_file() and path.stat().st_size <= 600_000:
                files.append(path)
    return sorted(set(files))


def read_text(path: Path) -> str:
    try:
        return path.read_text(encoding="utf-8", errors="ignore")
    except OSError:
        return ""


def build_corpus(paths: Iterable[Path]) -> dict[Path, str]:
    return {path: read_text(path) for path in paths}


def count_hits(motif: Motif, corpus: dict[Path, str]) -> tuple[int, list[str]]:
    total = 0
    sources: list[tuple[int, str]] = []
    for path, text in corpus.items():
        lower = text.lower()
        score = 0
        for kw in motif.keywords:
            score += lower.count(kw.lower())
        if score:
            total += score
            sources.append((score, str(path.relative_to(ROOT))))
    sources.sort(reverse=True)
    return total, [f"{p}:{s}" for s, p in sources[:4]]


def payload_score(motif: Motif, hit_count: int) -> int:
    score = sum(PAYLOAD_WEIGHTS[p] for p in motif.payloads)
    score += min(10, int(math.log2(hit_count + 1)) if hit_count else 0)
    score -= RISK_PENALTY[motif.risk]
    return score


def tournament_order(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    return sorted(
        rows,
        key=lambda r: (
            int(r["score"]),
            int(r["payload_count"]),
            int(r["hit_count"]),
            str(r["id"]),
        ),
        reverse=True,
    )


def count_directed_triangles(order: list[dict[str, object]]) -> int:
    # The orientation is induced by the total order, so this is a sanity check.
    rank = {row["id"]: i for i, row in enumerate(order)}
    cycles = 0
    ids = [row["id"] for row in order]
    for a, b, c in itertools.combinations(ids, 3):
        # Edges point from lower rank index to higher payload priority.
        edges = {
            (a, b): rank[a] < rank[b],
            (b, c): rank[b] < rank[c],
            (c, a): rank[c] < rank[a],
        }
        if all(edges.values()) or not any(edges.values()):
            cycles += 1
    return cycles


def main() -> None:
    corpus = build_corpus(candidate_files())
    rows: list[dict[str, object]] = []
    for motif in MOTIFS:
        hit_count, sources = count_hits(motif, corpus)
        rows.append(
            {
                "id": motif.ident,
                "label": motif.label,
                "family": motif.family,
                "pattern": motif.pattern,
                "payloads": motif.payloads,
                "payload_count": len(motif.payloads),
                "destroyed": motif.destroyed,
                "repair": motif.repair,
                "risk": motif.risk,
                "hit_count": hit_count,
                "sources": sources,
                "score": payload_score(motif, hit_count),
            }
        )

    order = tournament_order(rows)
    payload_counter: Counter[str] = Counter()
    family_counter: Counter[str] = Counter()
    risk_counter: Counter[str] = Counter()
    for row in rows:
        family_counter[str(row["family"])] += 1
        risk_counter[str(row["risk"])] += 1
        for payload in row["payloads"]:  # type: ignore[union-attr]
            payload_counter[payload] += 1

    score_hist = Counter(int(row["score"]) for row in rows)
    directed_3cycles = count_directed_triangles(order)
    hamiltonian_path = " -> ".join(str(row["id"]) for row in order)

    print("HYP-3226 small-pattern adjacency atlas")
    print("=" * 72)
    print(f"repo_files_scanned={len(corpus)}")
    print(f"motifs={len(rows)}")
    print(f"families={len(family_counter)}")
    print(f"risk_hist={dict(sorted(risk_counter.items()))}")
    print(f"score_hist={dict(sorted(score_hist.items()))}")
    print(f"directed_3cycles={directed_3cycles}")
    print("hamiltonian_path_count=1")
    print(f"priority_path={hamiltonian_path}")
    print()

    print("Payload Coverage")
    print("-" * 72)
    for payload, count in payload_counter.most_common():
        print(f"{payload:18s} {count:2d}")
    print()

    print("Top Motifs By Payload-Retention Score")
    print("-" * 72)
    print("keyword_hit_count=repo-local triage support; generic counts are capped and not proof strength.")
    for row in order[:24]:
        payloads = ",".join(row["payloads"])  # type: ignore[arg-type]
        print(
            f"{row['id']:>4s} score={row['score']:>2} risk={row['risk']:<7s} "
            f"hits={row['hit_count']:<4} family={row['family']} :: {row['label']}"
        )
        print(f"     payloads={payloads}")
        print(f"     repair={row['repair']}")
    print()

    print("Pattern Ledger")
    print("-" * 72)
    for row in order:
        sources = "; ".join(row["sources"]) if row["sources"] else "no local keyword hit"
        payloads = ",".join(row["payloads"])  # type: ignore[arg-type]
        print(f"{row['id']} | {row['label']} | family={row['family']} | risk={row['risk']}")
        print(f"  pattern: {row['pattern']}")
        print(f"  payloads: {payloads}")
        print(f"  destroyed: {row['destroyed']}")
        print(f"  repair_sidecar: {row['repair']}")
        print(f"  hits={row['hit_count']} sources={sources}")
    print()

    print("Small Numeric / Structural Signals")
    print("-" * 72)
    signals = [
        ("2-row skyline", "AP and doubled AP are the all-bank equality face; primitive form leaves AP."),
        ("11/12 traps", "11 non-AP arbitrary-exchange traps plus AP; Toeplitz slack discharges all 11."),
        ("3432/3431", "anchored bounded k=8 rows / primitive rows in the current exact bank."),
        ("39766/540225", "AP support projection value for consecutive speeds."),
        ("12882/17161", "worst central exchange-rate ratio for q3 gain vs q0+q6 loss."),
        ("6237419/8643600", "total covariance Sigma_kappa2 for consecutive k=8."),
        ("6237419/25930800", "ideal C6 Perron quotient lambda0=(1^T C 1)/6."),
        ("0.042304730706", "Toeplitz lambda-min margin at the AP row."),
        ("1/7", "seductive scalar refuted as exact associator law; survives as apex obstruction."),
        ("3.149364", "odd Worpitzky contribution dominates even fold in the k=8 packet."),
        ("1,4,1", "K3/Worpitzky Eulerian row; edge quotient warning."),
        ("1,11,11,1", "next Worpitzky row; odd-sidecar payload."),
        ("2,6,12,20,30,42", "doubled triangular numbers; usable only with partial-cube/simplex sidecar."),
        ("E/O legs", "E=x^2+5x+4 and O=x^2+4x+1 strictly interlace."),
        ("V7 double root", "Chebyshev/de Moivre cubic squared suggests magic-function dual."),
        ("F7 Fejer kernel", "positive-definite Delsarte sidecar with weights (7-|n|)_+."),
        ("K(1,q)=1/(7q)", "comb-overlap Gram kernel has an exact speed-1 row; q=13 is least resonant."),
        ("peeling recursion", "cap(P)=cap(P\\{1})-(1/7)(1-1/min(P\\{1})) for speed-1 rows <=13."),
        ("-37/1092,-61/588", "order-3 overlap residues for the binding j=4 and j=5 cap rows."),
        ("10q0+q3+10q6", "finite shell L_y magic dual; quartic contact at shells 1,2,4,5."),
        ("rho>=18.019", "cyclic PSD positivity overprices the central repair; guardrail only."),
        ("E7 Eisenstein", "E_7=(7E_2(7tau)-E_2(tau))/6 is a level-7 coefficient sidecar."),
        ("B7 cubic", "B^3-5B^2+6B-1=0 is a Beraha/Mahler height gauge, not an inequality."),
        ("7-|n|", "AP autocorrelation is the Fejer/subshift rank-one transfer state."),
        ("K(a,b)=g/(7ab)", "three-gap/Stern-Brocot recursion for the order-2 cap kernel."),
        ("K(a,13)", "antipode column K(a,13)=(2a-1)/(91a) is verified for a=1..12."),
        ("cap_k ladder", "cap_k=cap_{k-1}+k/C(2p,2) is the Faulhaber moment-order ladder."),
        ("depth=(p+1)/2", "LRC(2p) apex moment depth matches the cyclotomic-degree law for p=3,5,7."),
        ("scale-normal", "primitive projective shape plus first surviving coordinate is the route recursion."),
        ("K fold x1/2", "modulus-covariance gives K^(2n)/K^(n)=1/2 until the apex break."),
        ("chi3 cosets", "cubic mode cosets {1,6}/{2,5}/{3,4} are the de Moivre angles."),
        ("Phi_d grading", "recursion mode=(x-1)^depth*Phi_d separates moment depth from character."),
        ("A..G local", "signed recurrences require chart addresses before cancellation is legal."),
        ("|g(chi7)|^2=7", "Gauss-sum modulus equals the Lee-Yang apex zero and Fejer reserve."),
        ("disc=49", "Q(cos2pi/7) puts the binding cap rows on the 7^2 conductor."),
        ("2 heads", "n=14 has a 7-cap head and a 3^3 witness head, both depth 3."),
        ("1/23 -> 1/14", "Chen-Cusick supplies a floor-to-target lift; the 23/M=2/23 link is only bounded-bank coincidence."),
    ]
    for key, meaning in signals:
        print(f"{key:18s} {meaning}")
    print()

    print("Synthesis")
    print("-" * 72)
    print(
        "The useful small patterns cluster around seven proof payloads: "
        "normal-fan exposure, Toeplitz/moment curvature, covariance layers, "
        "ordered-tail pricing, finite trap discharge, HB/Perron gluing, and "
        "analytic equidistribution.  Patterns outside those payloads are not "
        "discarded, but they remain sidecars until they name the coordinate "
        "they preserve and the coordinate they destroy."
    )
    print(
        "With HYP-3225 now supplying the first trap-fingerprint table, the "
        "next high-value computation is symbolic discharge: prove the trap "
        "classes by first failed dictionary coordinate, Toeplitz slack, "
        "Green-current bottleneck type, conductance/Fiedler cut, "
        "Lorentzian/Plucker defect, M-matrix/Schur-complement debt, "
        "Worpitzky/HB sidecar debt, and Fejer/Delsarte F7 slack.  HYP-3227 "
        "makes that conductance graph payload concrete: it is a sidecar "
        "discharge ledger, not a terminal scalar ranking.  Incoming S75/S31ap "
        "adds a second layer of concrete payloads: comb-overlap Gram kernels, "
        "speed-1 peeling, order-3 overlap residues, the induction-base audit, "
        "and the Chen-Cusick/polyhedron-zonotope/Rosenfeld route split.  Incoming "
        "HYP-3228/HYP-3229 adds the finite shell magic coordinate "
        "10q0+q3+10q6, the Gamma0(7) Eisenstein coefficient sidecar, the "
        "Beraha/Mahler height gauge, the subshift transfer Perron-defect "
        "sidecar, and the Dirichlet-L/Stark denominator guardrail.  Incoming "
        "HYP-3230/HYP-3231/HYP-3216 then names the recursion layer: the "
        "three-gap/Stern-Brocot order-2 cap kernel, the scale-normal packet "
        "tower, the LRC(2p) moment-order depth law, and the 2-adic reflection "
        "fold.  Incoming HYP-3232/HYP-3217 adds the modulus-covariance apex "
        "break and the cyclotomic subfield/character-mode lattice, including "
        "the cubic de Moivre chi_3 mode.  Incoming HYP-3233/HYP-3234/HYP-3218/"
        "HYP-3235 then tightens that layer: recursion modes factor as "
        "(x-1)^depth*Phi_d, signed A..G recurrences need chart-change sidecars, "
        "AP autocorrelation is the self-dual Fejer equidistribution certificate, "
        "and the cap lives in the totally-real conductor-7 field with the binding "
        "rows carrying the 7^2 ramification debt."
    )


if __name__ == "__main__":
    main()
