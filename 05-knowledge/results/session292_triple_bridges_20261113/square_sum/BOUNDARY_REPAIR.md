# Prepublication simple-cycle boundary repair

The independent manuscript reviewer found that the initial all-word compiler's cycle flag imposed closure and distinct vertices before closure but omitted the undirected simple-cycle minimum length. Thus word `(4)` at N2 accepted the loop at2, and `(4,4)` at N3 accepted traversing one edge twice. The hexagon theorem was unaffected.

Repair: require word length at least3 in the theorem and both the direct and compiled functions, with explicit length-one and length-two hostile controls. This distinguishes simple cycles from closed walks; all other range and collision guards remain unchanged. Normal and optimized replays use the corrected source.
