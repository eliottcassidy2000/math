import ProcgenSelfieEdim

-- ListLemmas
#print axioms ProcgenSelfieEdim.nodup_range
#print axioms ProcgenSelfieEdim.length_filter_lt
#print axioms ProcgenSelfieEdim.pigeonhole
#print axioms ProcgenSelfieEdim.length_eq_of_same_members
#print axioms ProcgenSelfieEdim.eq_of_map_nodup
#print axioms ProcgenSelfieEdim.nodup_of_map_nodup
#print axioms ProcgenSelfieEdim.nodup_flatMap
#print axioms ProcgenSelfieEdim.nodup_of_allDistinct
#print axioms ProcgenSelfieEdim.all_of_all_eq_true
#print axioms ProcgenSelfieEdim.forall_lt_of_all_range
#print axioms ProcgenSelfieEdim.filter_cons_true
#print axioms ProcgenSelfieEdim.filter_cons_false
#print axioms ProcgenSelfieEdim.exists_of_length_filter_ne_zero
#print axioms ProcgenSelfieEdim.length_filter_or_le
#print axioms ProcgenSelfieEdim.length_filter_mono
#print axioms ProcgenSelfieEdim.length_filter_add_not
#print axioms ProcgenSelfieEdim.length_filter_false
#print axioms ProcgenSelfieEdim.length_filter_any_le
#print axioms ProcgenSelfieEdim.sum_le_mul
#print axioms ProcgenSelfieEdim.nodup_map_of_inj_on
#print axioms ProcgenSelfieEdim.memR_iff
#print axioms ProcgenSelfieEdim.nodup_of_distinctR
#print axioms ProcgenSelfieEdim.appR_eq
#print axioms ProcgenSelfieEdim.flatR_eq
#print axioms ProcgenSelfieEdim.downR_succ
#print axioms ProcgenSelfieEdim.mem_downR
#print axioms ProcgenSelfieEdim.length_downR
#print axioms ProcgenSelfieEdim.nodup_downR
#print axioms ProcgenSelfieEdim.countR_eq
#print axioms ProcgenSelfieEdim.allR_eq_true
#print axioms ProcgenSelfieEdim.anyR_eq_true
#print axioms ProcgenSelfieEdim.rsum_congr
#print axioms ProcgenSelfieEdim.rsum_add
#print axioms ProcgenSelfieEdim.rsum_zero
#print axioms ProcgenSelfieEdim.rsum_single

-- Hypercube
#print axioms ProcgenSelfieEdim.hamming_succ
#print axioms ProcgenSelfieEdim.hamming_le
#print axioms ProcgenSelfieEdim.hamming_comm
#print axioms ProcgenSelfieEdim.hamming_self
#print axioms ProcgenSelfieEdim.two_pow_succ
#print axioms ProcgenSelfieEdim.eq_of_hamming_eq_zero
#print axioms ProcgenSelfieEdim.hamming_congr_of_zero
#print axioms ProcgenSelfieEdim.hamming_adjacent
#print axioms ProcgenSelfieEdim.hamming_triangle
#print axioms ProcgenSelfieEdim.antipode_lt
#print axioms ProcgenSelfieEdim.antipode_antipode
#print axioms ProcgenSelfieEdim.hamming_antipode
#print axioms ProcgenSelfieEdim.hamming_antipode_antipode
#print axioms ProcgenSelfieEdim.isEdge_antipode
#print axioms ProcgenSelfieEdim.edgeDist_comm
#print axioms ProcgenSelfieEdim.hist_comm
#print axioms ProcgenSelfieEdim.edgeDist_le
#print axioms ProcgenSelfieEdim.edgeDist_antipode
#print axioms ProcgenSelfieEdim.hist_antipode
#print axioms ProcgenSelfieEdim.bit_zero
#print axioms ProcgenSelfieEdim.bit_succ
#print axioms ProcgenSelfieEdim.hamming_flip
#print axioms ProcgenSelfieEdim.flip_lt
#print axioms ProcgenSelfieEdim.edge_cases
#print axioms ProcgenSelfieEdim.mem_edgeList
#print axioms ProcgenSelfieEdim.isEdge_of_mem_edgeList
#print axioms ProcgenSelfieEdim.canonical_mem
#print axioms ProcgenSelfieEdim.hist_nil
#print axioms ProcgenSelfieEdim.hist_cons
#print axioms ProcgenSelfieEdim.hist_eq_zero_of_le
#print axioms ProcgenSelfieEdim.weighted_sum_eq
#print axioms ProcgenSelfieEdim.sum_map_one
#print axioms ProcgenSelfieEdim.rsum_hist
#print axioms ProcgenSelfieEdim.resolving_of_keys
#print axioms ProcgenSelfieEdim.hamR_eq
#print axioms ProcgenSelfieEdim.fmin_eq
#print axioms ProcgenSelfieEdim.keyR_eq
#print axioms ProcgenSelfieEdim.resolving_of_check
#print axioms ProcgenSelfieEdim.rsum_shift
#print axioms ProcgenSelfieEdim.rsum_palindrome_even
#print axioms ProcgenSelfieEdim.hist_not_palindrome
#print axioms ProcgenSelfieEdim.no_antipodal_collision

-- EdimQ6
#print axioms ProcgenSelfieEdim.paperSet_is_vertex_set
#print axioms ProcgenSelfieEdim.check_paperSet
#print axioms ProcgenSelfieEdim.edgeList_six
#print axioms ProcgenSelfieEdim.paperSet_resolving

-- CountingBound
#print axioms ProcgenSelfieEdim.eq_antipode_of_hamming_eq
#print axioms ProcgenSelfieEdim.touch_count_check
#print axioms ProcgenSelfieEdim.touch_count
#print axioms ProcgenSelfieEdim.touched_of_bad
#print axioms ProcgenSelfieEdim.mem_comps
#print axioms ProcgenSelfieEdim.comps_small
#print axioms ProcgenSelfieEdim.rsum_six
#print axioms ProcgenSelfieEdim.histVec6_mem_comps
#print axioms ProcgenSelfieEdim.not_resolving_of_length_le_six

-- Stabilizer
#print axioms ProcgenSelfieEdim.geodesic_step
#print axioms ProcgenSelfieEdim.hamming_map_le
#print axioms ProcgenSelfieEdim.automorphism_isometry
#print axioms ProcgenSelfieEdim.surj_of_inj
#print axioms ProcgenSelfieEdim.length_filter_comp
#print axioms ProcgenSelfieEdim.hist_invariant
#print axioms ProcgenSelfieEdim.isEdge_flipBit
#print axioms ProcgenSelfieEdim.flipBit_dist
#print axioms ProcgenSelfieEdim.trivial_stabilizer
#print axioms ProcgenSelfieEdim.paperSet_trivial_stabilizer

-- KeyBlocks
#print axioms ProcgenSelfieEdim.keyBR_eq
#print axioms ProcgenSelfieEdim.resolving_of_edgeKeys
#print axioms ProcgenSelfieEdim.blk_succ
#print axioms ProcgenSelfieEdim.eq_of_blocks
#print axioms ProcgenSelfieEdim.eq_of_eqListR
#print axioms ProcgenSelfieEdim.mergeR_perm
#print axioms ProcgenSelfieEdim.mergePairsR_perm
#print axioms ProcgenSelfieEdim.flatten_singletons
#print axioms ProcgenSelfieEdim.msortR_perm
#print axioms ProcgenSelfieEdim.strictFrom_lt
#print axioms ProcgenSelfieEdim.nodup_of_strictFrom
#print axioms ProcgenSelfieEdim.nodup_of_strictIncR
#print axioms ProcgenSelfieEdim.nodup_of_msortR

-- EdimQ7Data
#print axioms ProcgenSelfieEdim.q7set_is_vertex_set
#print axioms ProcgenSelfieEdim.q7_lengths

-- EdimQ7P0
#print axioms ProcgenSelfieEdim.q7_block_0
#print axioms ProcgenSelfieEdim.q7_block_1
#print axioms ProcgenSelfieEdim.q7_block_2
#print axioms ProcgenSelfieEdim.q7_block_3

-- EdimQ7P1
#print axioms ProcgenSelfieEdim.q7_block_4
#print axioms ProcgenSelfieEdim.q7_block_5
#print axioms ProcgenSelfieEdim.q7_block_6

-- EdimQ7
#print axioms ProcgenSelfieEdim.q7_distinct
#print axioms ProcgenSelfieEdim.q7_keys_eq
#print axioms ProcgenSelfieEdim.q7set_resolving

-- EdimQ8Data
#print axioms ProcgenSelfieEdim.q8set_is_vertex_set
#print axioms ProcgenSelfieEdim.q8_lengths

-- EdimQ8P0
#print axioms ProcgenSelfieEdim.q8_block_0
#print axioms ProcgenSelfieEdim.q8_block_1
#print axioms ProcgenSelfieEdim.q8_block_2
#print axioms ProcgenSelfieEdim.q8_block_3

-- EdimQ8P1
#print axioms ProcgenSelfieEdim.q8_block_4
#print axioms ProcgenSelfieEdim.q8_block_5
#print axioms ProcgenSelfieEdim.q8_block_6
#print axioms ProcgenSelfieEdim.q8_block_7

-- EdimQ8P2
#print axioms ProcgenSelfieEdim.q8_block_8
#print axioms ProcgenSelfieEdim.q8_block_9
#print axioms ProcgenSelfieEdim.q8_block_10
#print axioms ProcgenSelfieEdim.q8_block_11

-- EdimQ8P3
#print axioms ProcgenSelfieEdim.q8_block_12
#print axioms ProcgenSelfieEdim.q8_block_13
#print axioms ProcgenSelfieEdim.q8_block_14
#print axioms ProcgenSelfieEdim.q8_block_15

-- EdimQ8
#print axioms ProcgenSelfieEdim.q8_distinct
#print axioms ProcgenSelfieEdim.q8_keys_eq
#print axioms ProcgenSelfieEdim.q8set_resolving

-- EdimQ9Data
#print axioms ProcgenSelfieEdim.q9set_is_vertex_set
#print axioms ProcgenSelfieEdim.q9_lengths

-- EdimQ9P0
#print axioms ProcgenSelfieEdim.q9_block_0
#print axioms ProcgenSelfieEdim.q9_block_1
#print axioms ProcgenSelfieEdim.q9_block_2

-- EdimQ9P1
#print axioms ProcgenSelfieEdim.q9_block_3
#print axioms ProcgenSelfieEdim.q9_block_4
#print axioms ProcgenSelfieEdim.q9_block_5

-- EdimQ9P2
#print axioms ProcgenSelfieEdim.q9_block_6
#print axioms ProcgenSelfieEdim.q9_block_7
#print axioms ProcgenSelfieEdim.q9_block_8

-- EdimQ9P3
#print axioms ProcgenSelfieEdim.q9_block_9
#print axioms ProcgenSelfieEdim.q9_block_10
#print axioms ProcgenSelfieEdim.q9_block_11

-- EdimQ9P4
#print axioms ProcgenSelfieEdim.q9_block_12
#print axioms ProcgenSelfieEdim.q9_block_13
#print axioms ProcgenSelfieEdim.q9_block_14

-- EdimQ9P5
#print axioms ProcgenSelfieEdim.q9_block_15
#print axioms ProcgenSelfieEdim.q9_block_16
#print axioms ProcgenSelfieEdim.q9_block_17

-- EdimQ9P6
#print axioms ProcgenSelfieEdim.q9_block_18
#print axioms ProcgenSelfieEdim.q9_block_19
#print axioms ProcgenSelfieEdim.q9_block_20

-- EdimQ9P7
#print axioms ProcgenSelfieEdim.q9_block_21
#print axioms ProcgenSelfieEdim.q9_block_22
#print axioms ProcgenSelfieEdim.q9_block_23

-- EdimQ9P8
#print axioms ProcgenSelfieEdim.q9_block_24
#print axioms ProcgenSelfieEdim.q9_block_25
#print axioms ProcgenSelfieEdim.q9_block_26

-- EdimQ9P9
#print axioms ProcgenSelfieEdim.q9_block_27
#print axioms ProcgenSelfieEdim.q9_block_28
#print axioms ProcgenSelfieEdim.q9_block_29

-- EdimQ9P10
#print axioms ProcgenSelfieEdim.q9_block_30
#print axioms ProcgenSelfieEdim.q9_block_31
#print axioms ProcgenSelfieEdim.q9_block_32

-- EdimQ9P11
#print axioms ProcgenSelfieEdim.q9_block_33
#print axioms ProcgenSelfieEdim.q9_block_34
#print axioms ProcgenSelfieEdim.q9_block_35

-- EdimQ9
#print axioms ProcgenSelfieEdim.q9_distinct
#print axioms ProcgenSelfieEdim.q9_keys_eq
#print axioms ProcgenSelfieEdim.q9set_resolving

-- Tournament
#print axioms ProcgenSelfieEdim.tiling_down
#print axioms ProcgenSelfieEdim.tiling_up
#print axioms ProcgenSelfieEdim.tiling_lt
#print axioms ProcgenSelfieEdim.tiling_gt
#print axioms ProcgenSelfieEdim.four_cases
#print axioms ProcgenSelfieEdim.tiling_isTournament
#print axioms ProcgenSelfieEdim.switch_isTournament
#print axioms ProcgenSelfieEdim.switch_switch
#print axioms ProcgenSelfieEdim.not_eq_not
#print axioms ProcgenSelfieEdim.switch_compl
#print axioms ProcgenSelfieEdim.switch_congr
#print axioms ProcgenSelfieEdim.switch_loops_congr
#print axioms ProcgenSelfieEdim.tiling_congr
#print axioms ProcgenSelfieEdim.path_arc
#print axioms ProcgenSelfieEdim.tiles_of_switch_eq
#print axioms ProcgenSelfieEdim.loops_relation
#print axioms ProcgenSelfieEdim.selfie_fibre
#print axioms ProcgenSelfieEdim.loops_ne_compl
#print axioms ProcgenSelfieEdim.gauge_path
#print axioms ProcgenSelfieEdim.selfie_onto
#print axioms ProcgenSelfieEdim.reversedCount_parity
#print axioms ProcgenSelfieEdim.selfie_parameter_count
#print axioms ProcgenSelfieEdim.pairs_eq

-- HamPath
#print axioms ProcgenSelfieEdim.isPath_cons_cons
#print axioms ProcgenSelfieEdim.isPath_append_right
#print axioms ProcgenSelfieEdim.length_append_two
#print axioms ProcgenSelfieEdim.usesArc_iff
#print axioms ProcgenSelfieEdim.usesArc_mem
#print axioms ProcgenSelfieEdim.usesArc_arc
#print axioms ProcgenSelfieEdim.usesArc_self
#print axioms ProcgenSelfieEdim.extR_zero
#print axioms ProcgenSelfieEdim.extR_succ_cons
#print axioms ProcgenSelfieEdim.hamPathsR_zero
#print axioms ProcgenSelfieEdim.hamPathsR_succ
#print axioms ProcgenSelfieEdim.extR_shape
#print axioms ProcgenSelfieEdim.valid_cons
#print axioms ProcgenSelfieEdim.mem_extR
#print axioms ProcgenSelfieEdim.nodup_extR
#print axioms ProcgenSelfieEdim.mem_hamPathsR
#print axioms ProcgenSelfieEdim.nodup_hamPathsR
#print axioms ProcgenSelfieEdim.hpCount_unique
#print axioms ProcgenSelfieEdim.arcCount_unique

-- ArcParity
#print axioms ProcgenSelfieEdim.beq_false_of_ne
#print axioms ProcgenSelfieEdim.length_filter_eq_sum
#print axioms ProcgenSelfieEdim.rsum_listSum
#print axioms ProcgenSelfieEdim.sum_map_const
#print axioms ProcgenSelfieEdim.rsum_ite_out
#print axioms ProcgenSelfieEdim.rsum_rsum_point
#print axioms ProcgenSelfieEdim.ind_usesArc_cons_cons
#print axioms ProcgenSelfieEdim.arcs_of_path
#print axioms ProcgenSelfieEdim.sum_arcCount
#print axioms ProcgenSelfieEdim.rsum_one
#print axioms ProcgenSelfieEdim.arcTotal
#print axioms ProcgenSelfieEdim.rsum_mod_two_congr
#print axioms ProcgenSelfieEdim.arcCount_eq_zero_of_not
#print axioms ProcgenSelfieEdim.arcCount_mod_two_of_allOdd
#print axioms ProcgenSelfieEdim.pairs_parity
#print axioms ProcgenSelfieEdim.sum_arcCount_mod_two
#print axioms ProcgenSelfieEdim.oddArcs_parity
#print axioms ProcgenSelfieEdim.allOdd_mod_four
#print axioms ProcgenSelfieEdim.allEven_odd
#print axioms ProcgenSelfieEdim.isPath_deleteArc
#print axioms ProcgenSelfieEdim.shave_arc

-- AltSum
#print axioms ProcgenSelfieEdim.sgn_add
#print axioms ProcgenSelfieEdim.sgn_congr
#print axioms ProcgenSelfieEdim.hamming_eq_rsum
#print axioms ProcgenSelfieEdim.bit_lt_two
#print axioms ProcgenSelfieEdim.bit_flip
#print axioms ProcgenSelfieEdim.edgeDist_projection
#print axioms ProcgenSelfieEdim.alt_sum
#print axioms ProcgenSelfieEdim.alt_sum_hist
#print axioms ProcgenSelfieEdim.sgn_cases
#print axioms ProcgenSelfieEdim.alt_sum_two_values

-- TournamentCode
#print axioms ProcgenSelfieEdim.tourC_lt
#print axioms ProcgenSelfieEdim.tourC_gt
#print axioms ProcgenSelfieEdim.tourC_self
#print axioms ProcgenSelfieEdim.tourC_isTournament
#print axioms ProcgenSelfieEdim.nthB_append_left
#print axioms ProcgenSelfieEdim.nthB_append_right
#print axioms ProcgenSelfieEdim.bitB_bitsToNat
#print axioms ProcgenSelfieEdim.bitsToNat_lt
#print axioms ProcgenSelfieEdim.ascList_length
#print axioms ProcgenSelfieEdim.nthB_ascList
#print axioms ProcgenSelfieEdim.codeList_length
#print axioms ProcgenSelfieEdim.pairs_mono
#print axioms ProcgenSelfieEdim.codeList_nthB
#print axioms ProcgenSelfieEdim.pairs_formula
#print axioms ProcgenSelfieEdim.tourC_code

-- SelfieFinite
#print axioms ProcgenSelfieEdim.isPath_congr
#print axioms ProcgenSelfieEdim.isHamPath_congr
#print axioms ProcgenSelfieEdim.hpCount_congr
#print axioms ProcgenSelfieEdim.arcCount_congr
#print axioms ProcgenSelfieEdim.usesArcFrom_eq
#print axioms ProcgenSelfieEdim.usesArcR_eq
#print axioms ProcgenSelfieEdim.arcCountR_eq
#print axioms ProcgenSelfieEdim.lenR_eq
#print axioms ProcgenSelfieEdim.checkArcs_sound
#print axioms ProcgenSelfieEdim.qr7_isTournament_aux
#print axioms ProcgenSelfieEdim.qr7_isTournament
#print axioms ProcgenSelfieEdim.qr7_check
#print axioms ProcgenSelfieEdim.qr7_counts
#print axioms ProcgenSelfieEdim.lift_lt
#print axioms ProcgenSelfieEdim.lift_ne
#print axioms ProcgenSelfieEdim.lift_inj
#print axioms ProcgenSelfieEdim.qr7del_isTournament
#print axioms ProcgenSelfieEdim.qr7del_check
#print axioms ProcgenSelfieEdim.qr7del_counts
#print axioms ProcgenSelfieEdim.qr7del_antipodal_arcs
#print axioms ProcgenSelfieEdim.qr7del_allOdd
#print axioms ProcgenSelfieEdim.hasEvenArc_sound
#print axioms ProcgenSelfieEdim.checkFrom_sound
#print axioms ProcgenSelfieEdim.exists_even_arc_of_codes
#print axioms ProcgenSelfieEdim.three_codes
#print axioms ProcgenSelfieEdim.four_codes
#print axioms ProcgenSelfieEdim.five_chunk_0
#print axioms ProcgenSelfieEdim.five_chunk_1
#print axioms ProcgenSelfieEdim.five_chunk_2
#print axioms ProcgenSelfieEdim.five_chunk_3
#print axioms ProcgenSelfieEdim.five_chunk_4
#print axioms ProcgenSelfieEdim.five_chunk_5
#print axioms ProcgenSelfieEdim.five_chunk_6
#print axioms ProcgenSelfieEdim.five_chunk_7
#print axioms ProcgenSelfieEdim.five_codes
#print axioms ProcgenSelfieEdim.no_allOdd_three_to_five

-- Shaved
#print axioms ProcgenSelfieEdim.c3c3_isTournament_aux
#print axioms ProcgenSelfieEdim.c3c3_isTournament
#print axioms ProcgenSelfieEdim.lastOr_nil
#print axioms ProcgenSelfieEdim.lastOr_cons
#print axioms ProcgenSelfieEdim.eq_append_lastOr
#print axioms ProcgenSelfieEdim.lastOr_append
#print axioms ProcgenSelfieEdim.containsH_of_mem
#print axioms ProcgenSelfieEdim.hasH_sound
#print axioms ProcgenSelfieEdim.anyR_of_mem
#print axioms ProcgenSelfieEdim.anyR_headBeatsLast_of_containsH
#print axioms ProcgenSelfieEdim.containsH_congr
#print axioms ProcgenSelfieEdim.four_H_codes
#print axioms ProcgenSelfieEdim.every_four_contains_H4
#print axioms ProcgenSelfieEdim.nthR_mem
#print axioms ProcgenSelfieEdim.nthR_inj
#print axioms ProcgenSelfieEdim.isPath_complete
#print axioms ProcgenSelfieEdim.findIso_sound
#print axioms ProcgenSelfieEdim.five_H_chunk_0
#print axioms ProcgenSelfieEdim.five_H_chunk_1
#print axioms ProcgenSelfieEdim.five_H_chunk_2
#print axioms ProcgenSelfieEdim.five_H_chunk_3
#print axioms ProcgenSelfieEdim.five_H_chunk_4
#print axioms ProcgenSelfieEdim.five_H_chunk_5
#print axioms ProcgenSelfieEdim.five_H_chunk_6
#print axioms ProcgenSelfieEdim.five_H_chunk_7
#print axioms ProcgenSelfieEdim.five_H_codes
#print axioms ProcgenSelfieEdim.c3c3_avoids_check
#print axioms ProcgenSelfieEdim.c3c3_avoids_H5
#print axioms ProcgenSelfieEdim.isPath_map
#print axioms ProcgenSelfieEdim.containsH_of_iso
#print axioms ProcgenSelfieEdim.avoids_H5_iff

-- CycleCount
#print axioms ProcgenSelfieEdim.lastOr_mem
#print axioms ProcgenSelfieEdim.containsH_iff_copiesH_pos
#print axioms ProcgenSelfieEdim.closes_eq_not
#print axioms ProcgenSelfieEdim.split_closes
#print axioms ProcgenSelfieEdim.rotN_append
#print axioms ProcgenSelfieEdim.rotN_eq
#print axioms ProcgenSelfieEdim.isPath_snoc
#print axioms ProcgenSelfieEdim.rot1_closing
#print axioms ProcgenSelfieEdim.rotN_closing
#print axioms ProcgenSelfieEdim.pos_unique
#print axioms ProcgenSelfieEdim.rotN_zero_head
#print axioms ProcgenSelfieEdim.rotN_zero_split
#print axioms ProcgenSelfieEdim.rotN_nodup
#print axioms ProcgenSelfieEdim.rot_index_unique
#print axioms ProcgenSelfieEdim.rotN_inj
#print axioms ProcgenSelfieEdim.startsAtZero_iff
#print axioms ProcgenSelfieEdim.zero_mem_of_hamPath
#print axioms ProcgenSelfieEdim.closing_count
#print axioms ProcgenSelfieEdim.copiesH_add
#print axioms ProcgenSelfieEdim.c3c3_counts

-- Circulant
#print axioms ProcgenSelfieEdim.refl_spec
#print axioms ProcgenSelfieEdim.cdist_spec
#print axioms ProcgenSelfieEdim.refl_lt
#print axioms ProcgenSelfieEdim.refl_refl
#print axioms ProcgenSelfieEdim.cdist_refl
#print axioms ProcgenSelfieEdim.circ_refl
#print axioms ProcgenSelfieEdim.isPath_iff_pairs
#print axioms ProcgenSelfieEdim.isPath_reverse
#print axioms ProcgenSelfieEdim.nodup_reverse_of
#print axioms ProcgenSelfieEdim.length_filter_ne
#print axioms ProcgenSelfieEdim.even_of_involution
#print axioms ProcgenSelfieEdim.circulant_arcCount_even
#print axioms ProcgenSelfieEdim.qr7_eq_circ
#print axioms ProcgenSelfieEdim.qr7_allEven

-- HPExist
#print axioms ProcgenSelfieEdim.mem_insertV
#print axioms ProcgenSelfieEdim.nodup_insertV
#print axioms ProcgenSelfieEdim.head_insertV
#print axioms ProcgenSelfieEdim.isPath_insertV
#print axioms ProcgenSelfieEdim.exists_path_on
#print axioms ProcgenSelfieEdim.exists_hamPath
#print axioms ProcgenSelfieEdim.hpCount_pos
#print axioms ProcgenSelfieEdim.usesArc_succ_unique
#print axioms ProcgenSelfieEdim.usesArc_pred_unique
#print axioms ProcgenSelfieEdim.usesArc_last
#print axioms ProcgenSelfieEdim.no_covered_universal_arc

-- SwitchSum
#print axioms ProcgenSelfieEdim.ffact_self
#print axioms ProcgenSelfieEdim.sum_map_ite
#print axioms ProcgenSelfieEdim.count_not_mem
#print axioms ProcgenSelfieEdim.length_extR_complete
#print axioms ProcgenSelfieEdim.hpCount_complete
#print axioms ProcgenSelfieEdim.mem_allBools
#print axioms ProcgenSelfieEdim.nodup_allBools
#print axioms ProcgenSelfieEdim.length_allBools
#print axioms ProcgenSelfieEdim.nthB_ext
#print axioms ProcgenSelfieEdim.isPathB_iff
#print axioms ProcgenSelfieEdim.switch_pair
#print axioms ProcgenSelfieEdim.bool_flip
#print axioms ProcgenSelfieEdim.isPath_switch_congr
#print axioms ProcgenSelfieEdim.loops_along
#print axioms ProcgenSelfieEdim.loopsAlong_path
#print axioms ProcgenSelfieEdim.mem_of_hamPath
#print axioms ProcgenSelfieEdim.two_loop_sets
#print axioms ProcgenSelfieEdim.count_loop_sets
#print axioms ProcgenSelfieEdim.length_filter_cons
#print axioms ProcgenSelfieEdim.sum_map_add
#print axioms ProcgenSelfieEdim.sum_filter_swap
#print axioms ProcgenSelfieEdim.switch_sum

-- Redei
#print axioms ProcgenSelfieEdim.isPath_append_cons
#print axioms ProcgenSelfieEdim.isPath_prefix
#print axioms ProcgenSelfieEdim.isPath_snoc_iff
#print axioms ProcgenSelfieEdim.isPath_congr_pairs
#print axioms ProcgenSelfieEdim.lt_of_mem_ne
#print axioms ProcgenSelfieEdim.count_head
#print axioms ProcgenSelfieEdim.count_last
#print axioms ProcgenSelfieEdim.forceArc_ab
#print axioms ProcgenSelfieEdim.forceArc_ba
#print axioms ProcgenSelfieEdim.forceArc_other
#print axioms ProcgenSelfieEdim.forceArc_isTournament
#print axioms ProcgenSelfieEdim.forceArc_sameOn
#print axioms ProcgenSelfieEdim.eraseN_of_not_mem
#print axioms ProcgenSelfieEdim.eraseN_cons_self
#print axioms ProcgenSelfieEdim.eraseN_cons_ne
#print axioms ProcgenSelfieEdim.eraseN_append
#print axioms ProcgenSelfieEdim.split_around
#print axioms ProcgenSelfieEdim.nodup_insert_middle
#print axioms ProcgenSelfieEdim.mem_of_snoc_eq
#print axioms ProcgenSelfieEdim.isPath_splice
#print axioms ProcgenSelfieEdim.eraseN_splice
#print axioms ProcgenSelfieEdim.count_interior
#print axioms ProcgenSelfieEdim.lastIs_snoc
#print axioms ProcgenSelfieEdim.usesArc_into
#print axioms ProcgenSelfieEdim.usesArc_from
#print axioms ProcgenSelfieEdim.ind_and_false_left
#print axioms ProcgenSelfieEdim.ind_and_false_right
#print axioms ProcgenSelfieEdim.hp_split
#print axioms ProcgenSelfieEdim.sum_map_one'
#print axioms ProcgenSelfieEdim.hp_decomp
#print axioms ProcgenSelfieEdim.rsum_comm
#print axioms ProcgenSelfieEdim.mul_list_sum
#print axioms ProcgenSelfieEdim.list_sum_mod_two
#print axioms ProcgenSelfieEdim.rsum_rsum_weighted_point
#print axioms ProcgenSelfieEdim.weighted_arcs
#print axioms ProcgenSelfieEdim.odd_contribution
#print axioms ProcgenSelfieEdim.not_pred_and_succ
#print axioms ProcgenSelfieEdim.deleteArc_force_sameOn
#print axioms ProcgenSelfieEdim.redei
#print axioms ProcgenSelfieEdim.allOdd_mod_four_of_tournament
#print axioms ProcgenSelfieEdim.allEven_odd_of_tournament
#print axioms ProcgenSelfieEdim.oddArcs_parity_of_tournament
#print axioms ProcgenSelfieEdim.exists_odd_arc_of_even
#print axioms ProcgenSelfieEdim.shave_keeps_odd_iff
#print axioms ProcgenSelfieEdim.copiesH_odd_of_even
#print axioms ProcgenSelfieEdim.containsH_of_even

-- ConstantH
#print axioms ProcgenSelfieEdim.s2_eq
#print axioms ProcgenSelfieEdim.s2_zero
#print axioms ProcgenSelfieEdim.s2_double
#print axioms ProcgenSelfieEdim.s2_double_add_one
#print axioms ProcgenSelfieEdim.s2_eq_zero
#print axioms ProcgenSelfieEdim.pow_two_of_s2_eq_one
#print axioms ProcgenSelfieEdim.v2_mul
#print axioms ProcgenSelfieEdim.v2_odd
#print axioms ProcgenSelfieEdim.v2_double
#print axioms ProcgenSelfieEdim.s2_succ
#print axioms ProcgenSelfieEdim.fact_pos
#print axioms ProcgenSelfieEdim.legendre_two
#print axioms ProcgenSelfieEdim.constant_switching_class
#print axioms ProcgenSelfieEdim.constant_class_check
#print axioms ProcgenSelfieEdim.constant_class_example

-- ParityBreak
#print axioms ProcgenSelfieEdim.succ_or_last
#print axioms ProcgenSelfieEdim.pred_of_last
#print axioms ProcgenSelfieEdim.out_arcs_sum
#print axioms ProcgenSelfieEdim.end_sum
#print axioms ProcgenSelfieEdim.arc_split
#print axioms ProcgenSelfieEdim.delete_two
#print axioms ProcgenSelfieEdim.parity_break_two

-- CollatzDrop
#print axioms ProcgenSelfieEdim.aux_spec
#print axioms ProcgenSelfieEdim.v2_oddPart_of_eq
#print axioms ProcgenSelfieEdim.exists_decomp
#print axioms ProcgenSelfieEdim.decomp
#print axioms ProcgenSelfieEdim.oddPart_two_pow_mul
#print axioms ProcgenSelfieEdim.syr_decomp
#print axioms ProcgenSelfieEdim.drop_identity
#print axioms ProcgenSelfieEdim.two_pow_sub_three_ne_zero
#print axioms ProcgenSelfieEdim.drop_forward
#print axioms ProcgenSelfieEdim.drop_injective
#print axioms ProcgenSelfieEdim.drop_backward
#print axioms ProcgenSelfieEdim.admissible_iff
#print axioms ProcgenSelfieEdim.pow_two_ge_four
#print axioms ProcgenSelfieEdim.admissible_of_neg
#print axioms ProcgenSelfieEdim.admissible_of_nonneg
#print axioms ProcgenSelfieEdim.mod_four_iff
#print axioms ProcgenSelfieEdim.labelF_even
#print axioms ProcgenSelfieEdim.labelF_four_j_one
#print axioms ProcgenSelfieEdim.labelF_microcosm
#print axioms ProcgenSelfieEdim.labelF_eq_sub_drop
#print axioms ProcgenSelfieEdim.drop_two_regular
#print axioms ProcgenSelfieEdim.owner_K_values
#print axioms ProcgenSelfieEdim.owner_S_values
#print axioms ProcgenSelfieEdim.two_pow_mono
#print axioms ProcgenSelfieEdim.not_admissible_of_large
#print axioms ProcgenSelfieEdim.multiplicity_examples
