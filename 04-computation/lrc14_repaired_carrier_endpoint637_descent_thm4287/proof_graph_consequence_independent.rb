#!/usr/bin/env ruby
# Independent Ruby set replay of the THM-4287 typed proof-graph subtraction.

require 'digest'
require 'fileutils'
require 'set'

FNV_OFFSET = 0xcbf29ce484222325
FNV_PRIME = 0x100000001b3
MASK64 = (1 << 64) - 1

EXPECTED_INPUT_COUNT = 22_682
EXPECTED_INPUT_FNV = 0xf7563445f15efebf
EXPECTED_INPUT_SHA = '7d0044bc4c32f08b9d09dca420171a05666d26e03f38fbc48a9baa1fcb027102'
EXPECTED_CARRIER = [[100, 637], [294, 637], [520, 637]].freeze
EXPECTED_CARRIER_FNV = 0xf4b48b16a28826ad
EXPECTED_CARRIER_SHA = 'e31697893a8a4dbeb4b22ea62d7ef2bfed7fd55bbe35138a1bfc6f460638601b'
EXPECTED_SIGNATURE_COUNT = 33
EXPECTED_SIGNATURE_FNV = 0xf11425c1e1b17094
EXPECTED_SIGNATURE_SHA = '5b5de4e502802287ff2c504d34cfc4de21c10adffd1694c4bfa1c9a1d735c85b'
EXPECTED_OVERLAP = [[520, 637]].freeze
EXPECTED_OVERLAP_FNV = 0x3cefea79822890f4
EXPECTED_OVERLAP_SHA = '0d8ccfa43b9f309187ba867180e667eb9aac32341991ede0ea34a4f4b3b2fc12'
EXPECTED_UNION_COUNT = 35
EXPECTED_UNION_FNV = 0x0b5d8d9c28d7325d
EXPECTED_UNION_SHA = '1d36b9e6899496f93d44785d389523d4f3a85efe4eaadaf2e1297904116bbd78'
EXPECTED_FINAL_COUNT = 22_647
EXPECTED_FINAL_FNV = 0xdf5374d4aca67677
EXPECTED_FINAL_SHA = '14f9be0d9472bc573e582ec6f4cb92c7def6f583f6afaf0b747f2a9713330317'
EXPECTED_FINAL_MAX = 636
EXPECTED_FINAL_TOP = [
  [100, 636], [256, 636], [294, 636],
  [338, 636], [372, 636], [384, 636]
].freeze

def fail!(message)
  raise message
end

def raw_rows(rows)
  rows.map { |q, r| "#{q},#{r}\n" }.join
end

def read_pairs(path)
  data = File.binread(path)
  fail!("#{path}: ledger is not ASCII") unless data.ascii_only?
  rows = data.lines(chomp: true).map.with_index do |line, index|
    fail!("#{path}:#{index + 1}: blank row") if line.empty?
    fields = line.split(',', -1)
    fail!("#{path}:#{index + 1}: malformed pair") unless
      fields.length == 2 && fields.all? { |field| field.match?(/\A[0-9]+\z/) }
    pair = fields.map { |field| Integer(field, 10) }
    fail!("#{path}:#{index + 1}: invalid ordered pair") unless
      pair[0].positive? && pair[0] < pair[1]
    pair.freeze
  end
  fail!("#{path}: ledger not strictly ordered") unless
    rows.each_cons(2).all? { |a, b| (a <=> b) == -1 }
  fail!("#{path}: ledger serialization is not canonical") unless data == raw_rows(rows)
  rows
end

def fnv(rows)
  state = FNV_OFFSET
  rows.flatten.each do |word|
    8.times do |byte|
      state ^= (word >> (8 * byte)) & 0xff
      state = (state * FNV_PRIME) & MASK64
    end
  end
  state
end

def sha(rows)
  Digest::SHA256.hexdigest(raw_rows(rows))
end

def identity(label, rows)
  format('%s count=%d fnv=%016x sha256=%s', label, rows.length, fnv(rows), sha(rows))
end

fail!('usage: proof_graph_consequence_independent.rb U4283 CARRIER637 SIGNATURE OUTPUT') unless ARGV.length == 4

residual = read_pairs(ARGV[0])
carrier = read_pairs(ARGV[1])
signature = read_pairs(ARGV[2])

fail!('THM-4283 final residual identity changed') unless
  residual.length == EXPECTED_INPUT_COUNT && fnv(residual) == EXPECTED_INPUT_FNV &&
  sha(residual) == EXPECTED_INPUT_SHA
fail!('endpoint-637 carrier closure identity changed') unless
  carrier == EXPECTED_CARRIER && fnv(carrier) == EXPECTED_CARRIER_FNV &&
  sha(carrier) == EXPECTED_CARRIER_SHA
fail!('signature-fibre closure identity changed') unless
  signature.length == EXPECTED_SIGNATURE_COUNT &&
  fnv(signature) == EXPECTED_SIGNATURE_FNV && sha(signature) == EXPECTED_SIGNATURE_SHA

residual_set = residual.to_set
carrier_set = carrier.to_set
signature_set = signature.to_set
fail!('a closure ledger is not typed in the THM-4283 residual') unless
  carrier_set.subset?(residual_set) && signature_set.subset?(residual_set)

input_max = residual.map(&:last).max
input_top = residual.select { |_, r| r == input_max }
fail!('carrier closure is not exactly the maximal endpoint layer') unless
  input_max == 637 && input_top == carrier

overlap = (carrier_set & signature_set).to_a.sort
closure_union = (carrier_set | signature_set).to_a.sort
removed = (residual_set & closure_union.to_set).to_a.sort
final_residual = (residual_set - closure_union.to_set).to_a.sort
fail!('closure-ledger overlap changed') unless
  overlap == EXPECTED_OVERLAP && fnv(overlap) == EXPECTED_OVERLAP_FNV &&
  sha(overlap) == EXPECTED_OVERLAP_SHA
fail!('typed closure union/removal changed') unless
  closure_union.length == EXPECTED_UNION_COUNT && fnv(closure_union) == EXPECTED_UNION_FNV &&
  sha(closure_union) == EXPECTED_UNION_SHA && removed == closure_union
fail!('post-THM-4287 residual count changed') unless
  final_residual.length == EXPECTED_FINAL_COUNT && fnv(final_residual) == EXPECTED_FINAL_FNV &&
  sha(final_residual) == EXPECTED_FINAL_SHA

final_max = final_residual.map(&:last).max
final_top = final_residual.select { |_, r| r == final_max }
fail!('post-THM-4287 residual boundary changed') unless
  final_max == EXPECTED_FINAL_MAX && final_top == EXPECTED_FINAL_TOP

FileUtils.mkdir_p(File.dirname(ARGV[3]))
File.binwrite(ARGV[3], raw_rows(final_residual))

puts 'THM4287_PROOF_GRAPH_CONSEQUENCE_V1'
[
  ['INPUT_U4283', residual],
  ['CARRIER_ENDPOINT637', carrier],
  ['SIGNATURE_FIBRE_CLOSURE', signature],
  ['OVERLAP', overlap],
  ['CLOSURE_UNION', closure_union],
  ['REMOVED', removed],
  ['FINAL_RESIDUAL', final_residual]
].each { |label, rows| puts identity(label, rows) }
puts "INPUT_BOUNDARY max_r=#{input_max} top=#{input_top.map { |q, r| "#{q},#{r}" }.join(';')}"
puts "FINAL_BOUNDARY max_r=#{final_max} top=#{final_top.map { |q, r| "#{q},#{r}" }.join(';')}"
puts 'SCOPE TYPED_RESIDUAL_SUBTRACTION_NO_PHYSICAL_ENTRY_NO_LRC14_CONSEQUENCE'
puts 'VERDICT PASS EXACT_THM4287_PROOF_GRAPH_CONSEQUENCE'
