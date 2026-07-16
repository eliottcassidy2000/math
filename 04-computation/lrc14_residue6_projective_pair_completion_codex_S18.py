#!/usr/bin/env python3
"""Solver-free exact verifier for THM-909 and HYP-7065.

The seven-shift quotient of THM-908 leaves a finite function S_r on affine
lines in F_7^5.  Zero-coordinate uniformity, three unordered pair-ray rows,
and twenty affine-invariant labeled pair-ray rows cover every projective
residue direction.  All certificate weights below are literal integers;
SciPy is not imported or used.

Tournament Analysis is diagnostic.  Vertices are the twenty nonzero
collision orbits.  The pairwise observable is the difference of rigorous
upper bounds, the switch replaces the raw affine ceiling by the completed
pair-ray bound, and the gauge points toward the larger obstruction with
lexicographic tie completion.  This preserves worst-case residue ordering
but destroys within-cell wall chronology and the higher relation lattice.
Alternative vertex sets considered were runners, gaps, fixed sections,
boundaries, wall crossings, residues, affine cosets, Fourier modes, matroid
circuits, primitive relations, and proof obligations.  Projective directions
plus determinant-class pair sidecars are used because this is the smallest
quotient found that retains both the signed kernel and exact pair-ray laws.
"""

from collections import Counter
from fractions import Fraction
from functools import cache
from heapq import heapify, heappop, heappush
from itertools import combinations, combinations_with_replacement, permutations
from math import gcd, lcm

import numpy as np

from lrc14_residue6_seven_shift_sieve_codex_S18 import (
    PROPAGATION_SLACK,
    ceiling_census,
    projective_directions,
    scalar_tournament,
    state_arrays,
)


FIELD_SIZE = 7
MOVER_COUNT = 5
PAIR_TYPES = tuple(combinations_with_replacement(range(FIELD_SIZE), 2))
PAIR_TYPE_INDEX = {pair_type: index for index, pair_type in enumerate(PAIR_TYPES)}
PAIR_RAY_MINIMIZERS = (
    (1, 8), (1, 2), (1, 3), (1, 4), (1, 5), (1, 6), (2, 9),
    (2, 3), (2, 11), (2, 5), (2, 13), (3, 10), (3, 4), (3, 5),
    (3, 13), (4, 11), (4, 5), (4, 13), (5, 12), (5, 6), (6, 13),
)
EXPECTED_ZERO_MARGINAL_HISTOGRAM = {
    380: 5,
    252: 20,
    248: 10,
    204: 60,
    202: 60,
    200: 60,
    198: 20,
    196: 60,
    194: 60,
    188: 60,
    187: 120,
    184: 15,
    176: 30,
    170: 60,
    166: 120,
    164: 90,
    160: 120,
    159: 180,
    156: 110,
    148: 70,
    145: 20,
    138: 50,
    136: 30,
    126: 50,
    115: 20,
    81: 5,
}
SPARSE_PAIR_CERTIFICATES = {
    (1, 0, 0, 0, 0): (
        361,
        (
            408, -408, 408, -408, 408, -408, 204, 1224, -408, 2996,
            -408, 1224, 2276, 408, 2480, -204, -408, 5240, 1224, 204,
            1224, 14662, 408, -408, 0, 1224, 0, 0,
        ),
        Fraction(9626, 2527),
    ),
    (1, 6, 0, 0, 0): (
        2,
        (
            0, 11, 7, 24, 19, 11, 35, 22, -11, 13, -11, 26, 13, 22,
            17, 22, -11, 13, 22, 13, 13, 37, 22, -11, 13, 22, 13, 22,
        ),
        Fraction(97, 14),
    ),
    (1, 2, 4, 0, 0): (
        1,
        (
            0, 36, 36, 36, 12, 12, 24, 36, 3, 36, 3, 30, 30, 36,
            36, 30, 3, 30, 36, 36, 36, 36, 12, -3, 30, 12, 30, 24,
        ),
        Fraction(912, 35),
    ),
}
EXPECTED_NONZERO_ORBIT_COUNTS = {
    (1, 1, 1, 1, 4): 5,
    (1, 1, 1, 2, 2): 10,
    (1, 1, 1, 2, 3): 20,
    (1, 1, 1, 2, 4): 20,
    (1, 1, 1, 3, 4): 20,
    (1, 1, 1, 3, 6): 20,
    (1, 1, 1, 4, 4): 10,
    (1, 1, 1, 4, 5): 20,
    (1, 1, 1, 4, 6): 20,
    (1, 1, 1, 5, 6): 20,
    (1, 1, 2, 2, 3): 30,
    (1, 1, 2, 3, 3): 30,
    (1, 1, 2, 3, 4): 60,
    (1, 1, 2, 3, 5): 60,
    (1, 1, 2, 3, 6): 60,
    (1, 1, 2, 4, 5): 60,
    (1, 1, 2, 4, 6): 60,
    (1, 1, 2, 5, 5): 30,
    (1, 1, 3, 4, 5): 60,
    (1, 1, 3, 5, 6): 60,
}

INVARIANT_WEIGHT_DENOMINATOR = 100000
INVARIANT_WEIGHT_NUMERATORS = {
    (1, 1, 1, 1, 4): (200001, -399999, 28572, 114286, 114286, 28572, -399999, 200001, -399999, 28572, 114286, 114286, 28572, -399999, 200001, -399999, 28572, 114286, 114286, 28572, -399999, 171429, 85715, -171428, 85715, -599999, 171429, 257143, 200001, -399999, 28572, 114286, 114286, 28572, -399999, 200001, -399999, 28572, 114286, 114286, 28572, -399999, 171429, 85715, -171428, 85715, -599999, 171429, 257143, 200001, -399999, 28572, 114286, 114286, 28572, -399999, 171429, 85715, -171428, 85715, -599999, 171429, 257143, 2485715, 2400001, 2142858, 2400001, 1714286, 2485715, 2571429),
    (1, 1, 1, 2, 2): (-367914, -59893, -124064, 367915, 367915, -124064, -59893, -367914, -59893, -124064, 367915, 367915, -124064, -59893, -81890, 31602, 54185, 6061, 119553, -37662, -91847, -45382, -67965, -181457, -133333, 262266, 201299, -35425, -367914, -59893, -124064, 367915, 367915, -124064, -59893, -81890, 31602, 54185, 6061, 119553, -37662, -91847, 2060959, 2038376, 1924884, 1973008, 2368607, 2307640, 2070916, -81890, 31602, 54185, 6061, 119553, -37662, -91847, -45382, -67965, -181457, -133333, 262266, 201299, -35425, -620779, 185689, -258912, 20194, 281, -248955, 922486),
    (1, 1, 1, 2, 3): (144648, -41105, 19249, -122790, -122790, 19249, -41105, 144648, -41105, 19249, -122790, -122790, 19249, -41105, 302980, 164943, 312310, 24074, -264162, -189015, -351126, 200662, -34336, -148633, 191333, 34337, -222011, -21349, 144648, -41105, 19249, -122790, -122790, 19249, -41105, 2469899, 2331861, 2479228, 2190992, 1902756, 1977903, 1815793, 200662, -34336, -148633, 191333, 34337, -222011, -21349, 302980, 164943, 312310, 24074, -264162, -189015, -351126, 200662, -34336, -148633, 191333, 34337, -222011, -21349, -140280, -211573, 386751, -435027, -386750, 1081629, -294746),
    (1, 1, 1, 2, 4): (-202342, -186467, 152882, 100507, 100507, 152882, -186467, -202342, -186467, 152882, 100507, 100507, 152882, -186467, 26730, 420489, 75170, -104884, -75169, -315603, -26729, 2261865, 1969192, 2047347, 2069735, 1758473, 1976519, 2118175, -202342, -186467, 152882, 100507, 100507, 152882, -186467, 26730, 420489, 75170, -104884, -75169, -315603, -26729, 233108, -59565, 18589, 40977, -270284, -52238, 89417, 26730, 420489, 75170, -104884, -75169, -315603, -26729, 233108, -59565, 18589, 40977, -270284, -52238, 89417, -500426, 664974, 951134, -296220, 96589, -530684, -385365),
    (1, 1, 1, 3, 4): (-187489, -244005, 246804, 15117, 15117, 246804, -244005, -187489, -244005, 246804, 15117, 15117, 246804, -244005, 2125171, 2371549, 2267813, 2431109, 2150647, 2277546, 2203851, 262937, 264094, -139689, -124403, -343000, -45685, 125750, -187489, -244005, 246804, 15117, 15117, 246804, -244005, -135926, 110452, 6716, 170012, -110451, 16448, -57247, 262937, 264094, -139689, -124403, -343000, -45685, 125750, -135926, 110452, 6716, 170012, -110451, 16448, -57247, 262937, 264094, -139689, -124403, -343000, -45685, 125750, -994062, -636096, -497883, 1, -610849, -261845, 220903),
    (1, 1, 1, 3, 6): (-236008, -54338, 181670, -9326, -9326, 181670, -54338, -236008, -54338, 181670, -9326, -9326, 181670, -54338, -23425, 181292, 96559, 138617, -181291, -122465, -89284, -78943, 1, -227900, 28860, -42622, -35093, 219697, -236008, -54338, 181670, -9326, -9326, 181670, -54338, -23425, 181292, 96559, 138617, -181291, -122465, -89284, 1695723, 1774667, 1546767, 1803526, 1732044, 1739573, 1994363, -23425, 181292, 96559, 138617, -181291, -122465, -89284, -78943, 1, -227900, 28860, -42622, -35093, 219697, -359635, 1197, -162128, 160932, -15545, 606300, -231118),
    (1, 1, 1, 4, 4): (-503584, -397654, 655445, -310323, -94959, 569461, -96307, -490566, -83288, 539084, -125336, -297304, 668464, -384636, -7547, -13207, 32076, -18867, -141509, 43397, 105661, 2141510, 2135850, 2181133, 2130189, 2007548, 2192453, 2254717, -477547, -371617, 638087, -327681, -112318, 552103, -70269, -7547, -13207, 32076, -18867, -141509, 43397, 105661, -7547, -13207, 32076, -18867, -141509, 43397, 105661, -7547, -13207, 32076, -18867, -141509, 43397, 105661, -7547, -13207, 32076, -18867, -141509, 43397, 105661, 171699, -64150, -213207, 105661, 105661, -213207, -64150),
    (1, 1, 1, 4, 5): (-95963, -315343, 178526, 85325, 85325, 178526, -315343, -95963, -315343, 178526, 85325, 85325, 178526, -315343, 1989236, 2154716, 1850596, 2039752, 1683707, 2330472, 2056670, -89346, -248886, 248887, -299402, 181250, -67154, 274655, -95963, -315343, 178526, 85325, 85325, 178526, -315343, -25784, 139695, -164425, 24731, -331314, 315451, 41649, -89346, -248886, 248887, -299402, 181250, -67154, 274655, -25784, 139695, -164425, 24731, -331314, 315451, 41649, -89346, -248886, 248887, -299402, 181250, -67154, 274655, -531724, -465148, -335645, -176821, 335646, 818791, 354904),
    (1, 1, 1, 4, 6): (-153744, -337202, 386172, -114504, -114504, 386172, -337202, -153744, -337202, 386172, -114504, -114504, 386172, -337202, 40768, 35831, -37842, 2012, -322947, -149558, 431739, 2355971, 2857622, 2323951, 2389790, 2400141, 2609638, 2826769, -153744, -337202, 386172, -114504, -114504, 386172, -337202, 40768, 35831, -37842, 2012, -322947, -149558, 431739, -501651, 1, -533671, -467831, -457481, -247983, -30853, 40768, 35831, -37842, 2012, -322947, -149558, 431739, -501651, 1, -533671, -467831, -457481, -247983, -30853, -309717, -77407, 49927, 60278, 77408, 254616, -55101),
    (1, 1, 1, 5, 6): (-61557, -46167, 15390, 61558, 61558, 15390, -46167, -61557, -46167, 15390, 61558, 61558, 15390, -46167, -63797, 211487, 106739, -313198, 230884, -551224, 379112, -658226, 1, -501777, -263751, -642223, -634581, -536672, -61557, -46167, 15390, 61558, 61558, 15390, -46167, -63797, 211487, 106739, -313198, 230884, -551224, 379112, 2774371, 3432597, 2930819, 3168846, 2790374, 2798015, 2895925, -63797, 211487, 106739, -313198, 230884, -551224, 379112, -658226, 1, -501777, -263751, -642223, -634581, -536672, 15766, -88918, 309495, -291508, 477739, 380428, -802998),
    (1, 1, 2, 2, 3): (1, 66950, -66949, 1, 1, -66949, 66950, -269317, -200318, 348206, -45408, 14020, 245727, -92907, -269317, -200318, 348206, -45408, 14020, 245727, -92907, 101638, 556365, -514278, 859542, -258980, -274268, -470014, -269317, -200318, 348206, -45408, 14020, 245727, -92907, 1849534, 1918532, 2467056, 2073442, 2132870, 2364577, 2025943, 101638, 556365, -514278, 859542, -258980, -274268, -470014, 402934, 207986, -383121, -331789, -331789, -383121, 207986, -53844, -2457, -101119, -12442, 101120, 367954, -299209, -53844, -2457, -101119, -12442, 101120, 367954, -299209),
    (1, 1, 2, 3, 3): (1, 140563, 19058, -159620, -159620, 19058, 140563, 103669, -379917, 504403, 73046, -358312, 306873, -249758, 91375, 43843, 48483, 362921, 15367, -297076, -264910, 37806, 25083, 94680, 417493, 26692, -293892, -307860, 103669, -379917, 504403, 73046, -358312, 306873, -249758, 91375, 43843, 48483, 362921, 15367, -297076, -264910, 2239488, 2226764, 2296362, 2619175, 2228374, 1907790, 1893822, -161561, -122646, 28181, -634303, 344891, 645114, -99671, -221504, -39399, 84464, -593832, 152736, 752668, -135129, 123157, -80406, -17552, -154875, -14459, 17553, 3428),
    (1, 1, 2, 3, 4): (-101578, -31498, 121935, -52610, -52610, 121935, -31498, 76387, 161053, 307556, -284316, -5278, 123264, -378663, -119312, 88060, -46186, 583886, -88059, -149537, -268849, -205165, 106707, 27358, -134063, -387937, 336221, 256883, 76387, 161053, 307556, -284316, -5278, 123264, -378663, -119312, 88060, -46186, 583886, -88059, -149537, -268849, -205165, 106707, 27358, -134063, -387937, 336221, 256883, -22111, -146295, 82479, -359539, -82478, 865375, -337428, -476348, 188133, 31629, 244011, 43729, 223839, -254990, 1623579, 1997507, 2099316, 2130540, 2058439, 1832725, 1647446),
    (1, 1, 2, 3, 5): (-72742, 7110, 79852, -50589, -50589, 79852, 7110, 102744, 120324, 272934, -107488, -203830, -12835, -171846, 231114, 185589, -87725, 419622, -185588, -397061, -165948, -163613, -264418, 351698, -62815, 80414, -295535, 354272, 102744, 120324, 272934, -107488, -203830, -12835, -171846, 231114, 185589, -87725, 419622, -185588, -397061, -165948, 1876614, 1775810, 2391926, 1977413, 2120642, 1744693, 2394500, 35964, -107223, -239127, -321194, 289417, 649036, -306869, -276542, -129261, 1, 182741, -106445, -268520, -399466, -152596, 53717, -83013, 434, 111444, -12562, 82580),
    (1, 1, 2, 3, 6): (-173646, -102914, 70732, 119007, 119007, 70732, -102914, -234387, -95027, 353751, 106609, -140533, -11581, 21171, 2552943, 2904150, 2876773, 3283688, 2586696, 2603096, 2410616, -520169, 1, -358613, -594098, -662681, 166415, -348783, -234387, -95027, 353751, 106609, -140533, -11581, 21171, -192480, 158727, 131350, 538265, -158726, -142326, -334807, -520169, 1, -358613, -594098, -662681, 166415, -348783, 525173, -276301, 142179, -120485, -142178, 517273, -645657, -410290, 26976, 102460, 204556, -48088, -335169, 459559, 61699, 141649, -126530, -15118, -132547, -103023, 173874),
    (1, 1, 2, 4, 5): (-204240, -14396, -87723, 204241, 204241, -87723, -14396, 167655, -10749, 94207, -63588, 9982, 74339, -271842, -102995, 324133, -172060, 6474, -196721, 5555, 135617, 2025818, 1944740, 2235573, 1814287, 2034964, 1897765, 2677949, 167655, -10749, 94207, -63588, 9982, 74339, -271842, -102995, 324133, -172060, 6474, -196721, 5555, 135617, -64338, -145416, 145417, -275869, -55192, -192391, 587793, -183186, 294559, 63513, 264035, 176224, -287595, -327547, -139663, -696204, 1, -85874, -277976, 370206, -539124, 19795, -300034, -289165, -211682, 289166, 723399, -231477),
    (1, 1, 2, 4, 6): (-303186, -83431, 211402, 23623, 23623, 211402, -83431, -158078, 20322, 240252, -318037, -15669, 297716, -66503, 111804, 293152, -183832, -109318, -208479, -12361, 109037, 2394931, 2996363, 2261428, 2592764, 2127520, 2708260, 2697129, -158078, 20322, 240252, -318037, -15669, 297716, -66503, 111804, 293152, -183832, -109318, -208479, -12361, 109037, -601432, 1, -734934, -403599, -868842, -288102, -299234, -647540, 667908, -15662, 1714, 300279, -320646, 13949, -69984, -175443, 255105, -199615, 306190, -349261, 233012, -480236, 335975, -33550, 100519, -75541, 324714, -171878),
    (1, 1, 2, 5, 5): (-383223, -150835, 407672, -109044, -109044, 407672, -150835, -168123, -123539, 320598, -362878, -64562, 569723, -171215, -53326, 6354, -6353, -100554, 44052, -153880, 263711, -53326, 6354, -6353, -100554, 44052, -153880, 263711, -168123, -123539, 320598, -362878, -64562, 569723, -171215, -53326, 6354, -6353, -100554, 44052, -153880, 263711, -53326, 6354, -6353, -100554, 44052, -153880, 263711, -411684, -207484, 1, 396290, -314136, 178124, 358892, 1706551, 1910751, 2118235, 2514525, 1804099, 2296359, 2477126, -362089, 331792, 85224, -235970, -235970, 85224, 331792),
    (1, 1, 3, 4, 5): (-416578, -93068, -106923, 408282, 408282, -106923, -93068, -204830, 582947, -3702, 366737, -253003, 23314, -511460, -99288, 144485, -103049, -41435, -141360, 97854, 142796, -136385, 215035, -117881, 131572, -33742, -101966, 43370, -204830, 582947, -3702, 366737, -253003, 23314, -511460, -99288, 144485, -103049, -41435, -141360, 97854, 142796, -136385, 215035, -117881, 131572, -33742, -101966, 43370, 1534998, 2614494, 2188407, 2639957, 2060139, 1671454, 2015245, 152399, -237758, 108711, 163734, 245774, -160413, -272444, -189808, 557880, -151284, -394478, 185626, 162396, -170329),
    (1, 1, 3, 5, 6): (-73121, -37015, 456, 73122, 73122, 456, -37015, -137253, 167153, -162533, 330612, 12274, 53215, -263465, -240548, -198452, 198453, 17207, 145467, -223342, 301218, -209917, 1, -388956, -51038, -280401, 74975, 361442, 1863039, 2167445, 1837759, 2330904, 2012566, 2053507, 1736827, -240548, -198452, 198453, 17207, 145467, -223342, 301218, -209917, 1, -388956, -51038, -280401, 74975, 361442, -13592, -173123, -283694, 148032, 485069, -334216, 171528, 6771, 176616, -496891, 325977, -400634, 300558, 87607, -8044, 75833, 195179, -107570, -163644, 31739, -23488),
}


@cache
def labelled_pair_distribution(first_speed, second_speed):
    common_divisor = gcd(first_speed, second_speed)
    speeds = (
        first_speed // common_divisor,
        second_speed // common_divisor,
    )
    period_scale = lcm(*speeds)
    sectors = [0, 0]
    events = []
    for runner_index, speed in enumerate(speeds):
        event_step = period_scale // speed
        events.append((event_step, runner_index, 1, speed, event_step))
    heapify(events)
    numerators = [[0] * FIELD_SIZE for _ in range(FIELD_SIZE)]
    previous = 0
    while events:
        event_position = events[0][0]
        numerators[sectors[0]][sectors[1]] += event_position - previous
        while events and events[0][0] == event_position:
            _, runner_index, event_index, speed, event_step = heappop(events)
            sectors[runner_index] = (sectors[runner_index] + 1) % FIELD_SIZE
            if event_index < FIELD_SIZE * speed:
                next_index = event_index + 1
                heappush(
                    events,
                    (
                        next_index * event_step,
                        runner_index,
                        next_index,
                        speed,
                        event_step,
                    ),
                )
        previous = event_position
    return tuple(
        tuple(
            Fraction(numerators[first][second], FIELD_SIZE * period_scale)
            for second in range(FIELD_SIZE)
        )
        for first in range(FIELD_SIZE)
    )


def unordered_pair_distribution(first_speed, second_speed):
    labelled = labelled_pair_distribution(first_speed, second_speed)
    return {
        (first, second): labelled[first][second]
        if first == second
        else labelled[first][second] + labelled[second][first]
        for first, second in PAIR_TYPES
    }


def ordered_endpoint_representatives():
    endpoints = {}
    for first_residue in range(1, FIELD_SIZE):
        for second_residue in range(1, FIELD_SIZE):
            candidates = [
                (first * second, first, second)
                for first in range(1, 30)
                for second in range(1, 30)
                if first != second
                and gcd(first, second) == 1
                and first % FIELD_SIZE == first_residue
                and second % FIELD_SIZE == second_residue
            ]
            minimum_product = min(product_value for product_value, _, _ in candidates)
            first, second = min(
                (first, second)
                for product_value, first, second in candidates
                if product_value == minimum_product
            )
            endpoints[first_residue, second_residue] = (first, second)
    return endpoints


def labelled_ray_vertices(endpoints):
    uniform = tuple(
        tuple(Fraction(1, 49) for _ in range(FIELD_SIZE))
        for _ in range(FIELD_SIZE)
    )
    return {
        residue_pair: (
            first * second,
            labelled_pair_distribution(first, second),
        )
        for residue_pair, (first, second) in endpoints.items()
    }, uniform


def verify_labelled_ray_law(vertices, uniform, diameter=60):
    checked = 0
    for first_speed in range(1, diameter + 1):
        for second_speed in range(1, diameter + 1):
            if first_speed == second_speed:
                continue
            common_divisor = gcd(first_speed, second_speed)
            first = first_speed // common_divisor
            second = second_speed // common_divisor
            actual = labelled_pair_distribution(first, second)
            residues = (first % FIELD_SIZE, second % FIELD_SIZE)
            if 0 in residues:
                predicted = uniform
            else:
                endpoint_product, endpoint = vertices[residues]
                interpolation = Fraction(endpoint_product, first * second)
                predicted = tuple(
                    tuple(
                        (1 - interpolation) * uniform[source][target]
                        + interpolation * endpoint[source][target]
                        for target in range(FIELD_SIZE)
                    )
                    for source in range(FIELD_SIZE)
                )
            if actual != predicted:
                raise AssertionError((first_speed, second_speed, residues))
            checked += 1
    return checked


def line_data(direction, digits, powers, signed_kernel):
    pivot = next(index for index, residue in enumerate(direction) if residue)
    representatives = digits[digits[:, pivot] == 0]
    shifts = np.arange(FIELD_SIZE, dtype=np.int16)[:, None, None]
    line_indices = (
        representatives[None, :, :]
        + shifts * np.array(direction, dtype=np.int16)[None, None, :]
    ) % FIELD_SIZE @ powers
    line_sums = signed_kernel[line_indices].sum(axis=0)
    return representatives, line_sums


def canonical_direction(direction):
    return min(
        tuple((scalar * direction[index]) % FIELD_SIZE for index in permutation)
        for scalar in range(1, FIELD_SIZE)
        for permutation in permutations(range(MOVER_COUNT))
    )


def zero_marginal_numerator(direction, representatives, line_sums):
    return min(
        sum(
            int(line_sums[representatives[:, coordinate] == sector].max())
            for sector in range(FIELD_SIZE)
        )
        for coordinate, residue in enumerate(direction)
        if residue == 0
    )


def unordered_pair_vertices():
    uniform = {
        pair_type: Fraction(1 if pair_type[0] == pair_type[1] else 2, 49)
        for pair_type in PAIR_TYPES
    }
    return (uniform,) + tuple(
        unordered_pair_distribution(*speeds)
        for speeds in PAIR_RAY_MINIMIZERS
    )


def verify_sparse_certificate(direction, certificate, digits, powers, signed_kernel, pair_vertices):
    denominator, numerators, expected_pair_cap = certificate
    representatives, line_sums = line_data(direction, digits, powers, signed_kernel)
    zero_pairs = tuple(
        pair
        for pair in combinations(range(MOVER_COUNT), 2)
        if direction[pair[0]] == direction[pair[1]] == 0
    )
    minimum_slack = min(
        sum(
            numerators[
                PAIR_TYPE_INDEX[
                    tuple(sorted((int(state[first]), int(state[second]))))
                ]
            ]
            for first, second in zero_pairs
        )
        - denominator * int(line_sum)
        for state, line_sum in zip(representatives, line_sums)
    )
    pair_cap = max(
        sum(
            Fraction(numerator, denominator) * distribution[pair_type]
            for numerator, pair_type in zip(numerators, PAIR_TYPES)
        )
        for distribution in pair_vertices
    )
    if pair_cap != expected_pair_cap:
        raise AssertionError((direction, pair_cap, expected_pair_cap))
    return minimum_slack, len(zero_pairs) * pair_cap / 343


def invariant_distribution(direction, first, second, distribution):
    masses = [Fraction(0) for _ in range(FIELD_SIZE)]
    for source in range(FIELD_SIZE):
        for target in range(FIELD_SIZE):
            invariant = (
                direction[second] * source - direction[first] * target
            ) % FIELD_SIZE
            masses[invariant] += distribution[source][target]
    return tuple(masses)


def invariant_pair_vertices(direction, first, second, ordered_vertices, uniform):
    vertices = [invariant_distribution(direction, first, second, uniform)]
    for scalar in range(1, FIELD_SIZE):
        residues = (
            scalar * direction[first] % FIELD_SIZE,
            scalar * direction[second] % FIELD_SIZE,
        )
        vertices.append(
            invariant_distribution(
                direction,
                first,
                second,
                ordered_vertices[residues][1],
            )
        )
    return tuple(vertices)


def verify_invariant_certificate(direction, numerators, digits, powers, signed_kernel, ordered_vertices, uniform):
    representatives, line_sums = line_data(direction, digits, powers, signed_kernel)
    coordinate_pairs = tuple(combinations(range(MOVER_COUNT), 2))
    minimum_slack = min(
        sum(
            numerators[FIELD_SIZE * pair_index + (
                direction[second] * int(state[first])
                - direction[first] * int(state[second])
            ) % FIELD_SIZE]
            for pair_index, (first, second) in enumerate(coordinate_pairs)
        )
        - INVARIANT_WEIGHT_DENOMINATOR * int(line_sum)
        for state, line_sum in zip(representatives, line_sums)
    )
    pair_caps = []
    for pair_index, (first, second) in enumerate(coordinate_pairs):
        weights = tuple(
            Fraction(
                numerators[FIELD_SIZE * pair_index + invariant],
                INVARIANT_WEIGHT_DENOMINATOR,
            )
            for invariant in range(FIELD_SIZE)
        )
        vertices = invariant_pair_vertices(
            direction, first, second, ordered_vertices, uniform
        )
        pair_caps.append(
            max(
                sum(weight * mass for weight, mass in zip(weights, distribution))
                for distribution in vertices
            )
        )
    return minimum_slack, sum(pair_caps), sum(pair_caps) / 343


def check(label, condition):
    if not condition:
        raise AssertionError(label)
    print(f"PASS  {label}")


def main():
    print("THM-909 / HYP-7065: PROJECTIVE PAIR-RAY COMPLETION")
    print("=" * 78)
    digits, powers, signed_kernel, positive_majorant = state_arrays()
    census = ceiling_census(digits, powers, signed_kernel, positive_majorant)
    records = census["records"]

    print("\n[1] Exact labeled pair-ray law")
    endpoints = ordered_endpoint_representatives()
    ordered_vertices, uniform_labelled = labelled_ray_vertices(endpoints)
    checked_pairs = verify_labelled_ray_law(ordered_vertices, uniform_labelled)
    check("3,540 ordered speed pairs through 60", checked_pairs == 3540)
    check("36 ordered nonzero residue endpoints", len(endpoints) == 36)
    print("  product-minimal endpoint representatives:", endpoints)

    print("\n[2] Zero-coordinate uniform-marginal census")
    zero_histogram = Counter()
    zero_directions = []
    sparse_directions = []
    for direction, _ in projective_directions():
        if 0 not in direction:
            continue
        representatives, line_sums = line_data(direction, digits, powers, signed_kernel)
        numerator = zero_marginal_numerator(direction, representatives, line_sums)
        zero_histogram[numerator] += 1
        zero_directions.append(direction)
        if numerator > 204:
            sparse_directions.append(direction)
    check("1,505 zero-containing directions", len(zero_directions) == 1505)
    check("exact zero-marginal histogram", dict(zero_histogram) == EXPECTED_ZERO_MARGINAL_HISTOGRAM)
    check("1,470 directions satisfy U<=204", len(zero_directions) - len(sparse_directions) == 1470)
    print("  U(r) numerator histogram:", dict(sorted(zero_histogram.items(), reverse=True)))
    print("  uniform-marginal bound: 204/2401 =", float(Fraction(204, 2401)))

    print("\n[3] Three sparse zero-direction pair certificates")
    sparse_orbits = Counter(canonical_direction(direction) for direction in sparse_directions)
    sparse_certificate_orbits = {
        canonical_direction(direction): direction
        for direction in SPARSE_PAIR_CERTIFICATES
    }
    check("35 sparse directions", len(sparse_directions) == 35)
    check("three sparse permutation/projective orbits", set(sparse_orbits) == set(sparse_certificate_orbits))
    pair_vertices = unordered_pair_vertices()
    sparse_bounds = {}
    for direction, certificate in SPARSE_PAIR_CERTIFICATES.items():
        minimum_slack, bound = verify_sparse_certificate(
            direction, certificate, digits, powers, signed_kernel, pair_vertices
        )
        check(f"sparse pointwise row {direction}", minimum_slack >= 0)
        check(f"sparse ray bound {direction}", bound < PROPAGATION_SLACK)
        sparse_bounds[direction] = bound
        print("   ", direction, "orbit", sparse_orbits[canonical_direction(direction)], "minimum slack", minimum_slack, "-F6 bound", bound)

    print("\n[4] Twenty nonzero affine-invariant pair certificates")
    bad_nonzero = [
        direction
        for direction, (signed_ceiling, _) in records.items()
        if 0 not in direction and signed_ceiling >= 34
    ]
    orbit_counts = Counter(canonical_direction(direction) for direction in bad_nonzero)
    check("675 bad nonzero directions", len(bad_nonzero) == 675)
    check("twenty exact symmetry-orbit counts", dict(orbit_counts) == EXPECTED_NONZERO_ORBIT_COUNTS)
    check("twenty literal invariant certificates", set(INVARIANT_WEIGHT_NUMERATORS) == set(orbit_counts))
    invariant_bounds = {}
    invariant_ceiling_bounds = {}
    for direction in sorted(INVARIANT_WEIGHT_NUMERATORS):
        minimum_slack, ceiling_bound, residue_bound = verify_invariant_certificate(
            direction,
            INVARIANT_WEIGHT_NUMERATORS[direction],
            digits,
            powers,
            signed_kernel,
            ordered_vertices,
            uniform_labelled,
        )
        check(f"invariant pointwise row {direction}", minimum_slack >= 0)
        check(f"invariant ray bound {direction}", residue_bound < PROPAGATION_SLACK)
        invariant_ceiling_bounds[direction] = ceiling_bound
        invariant_bounds[direction] = residue_bound
        print("   ", direction, "orbit", orbit_counts[direction], "minimum numerator slack", minimum_slack, "E[S] cap", ceiling_bound, "-F6 cap", residue_bound)
    worst_direction = max(invariant_ceiling_bounds, key=invariant_ceiling_bounds.get)
    check("worst invariant orbit is (1,1,1,1,4)", worst_direction == (1, 1, 1, 1, 4))
    check("worst invariant ceiling bound", invariant_ceiling_bounds[worst_direction] == Fraction(3240009, 140000))

    print("\n[5] Universal limiting closure")
    raw_nonzero_closures = sum(
        0 not in direction and Fraction(signed_ceiling, 343) < PROPAGATION_SLACK
        for direction, (signed_ceiling, _) in records.items()
    )
    check("621 raw no-zero closures", raw_nonzero_closures == 621)
    check("all 2,801 projective directions covered", 1505 + raw_nonzero_closures + len(bad_nonzero) == 2801)
    branch_bounds = {
        "raw affine ceiling": Fraction(32, 343),
        "zero marginal": Fraction(204, 2401),
        "sparse pair rows": max(sparse_bounds.values()),
        "invariant pair rows": max(invariant_bounds.values()),
    }
    universal_bound = max(branch_bounds.values())
    print("  branch bounds:", branch_bounds)
    print("  universal limiting bound:", universal_bound, float(universal_bound))
    check("universal bound is 32/343", universal_bound == Fraction(32, 343))
    check("universal bound clears 0.097", universal_bound < PROPAGATION_SLACK)

    print("\n[6] Tournament Analysis: raw ceiling -> completed pair bound")
    vertices = tuple(sorted(INVARIANT_WEIGHT_NUMERATORS))
    raw_values = {direction: Fraction(records[direction][0]) for direction in vertices}
    completed_values = {direction: invariant_ceiling_bounds[direction] for direction in vertices}
    raw_tournament = scalar_tournament(vertices, raw_values)
    completed_tournament = scalar_tournament(vertices, completed_values)
    edge_flips = len(raw_tournament["edges"] ^ completed_tournament["edges"]) // 2
    print("  observable: rigorous obstruction-bound difference")
    print("  switch/gauge: raw affine ceiling -> completed pair-ray cap; larger wins")
    print("  raw fingerprint:", {key: value for key, value in raw_tournament.items() if key != "edges"})
    print("  completed fingerprint:", {key: value for key, value in completed_tournament.items() if key != "edges"})
    print("  switch edge flips:", edge_flips)
    check("both bound tournaments are transitive", raw_tournament["directed_triangles"] == completed_tournament["directed_triangles"] == 0)
    check("both tie Hamiltonian paths are unique", raw_tournament["hamiltonian_path_count"] == completed_tournament["hamiltonian_path_count"] == 1)
    print("  Assumption challenge: raw pair marginals fail before signed averaging;")
    print("  determinant-class pair sidecars suffice only after the seven-shift quotient.")

    print("\nVERDICT")
    print("  PROVED: every zero-containing projective direction closes.")
    print("  PROVED: all 675 nonzero collision directions close in twenty orbits.")
    print("  PROVED: universal limiting -F6<=32/343<0.097.")
    print("  OPEN: compose this limiting sign closure with the finite-t wall remainder.")


if __name__ == "__main__":
    main()
