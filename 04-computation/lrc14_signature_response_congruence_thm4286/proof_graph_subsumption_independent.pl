#!/usr/bin/env perl
# Independent exact-set replay of the current THM-4286 proof-graph placement.

use strict;
use warnings;
use Digest::SHA qw(sha256_hex);
use FindBin qw($Bin);
use Cwd qw(abs_path);
use File::Spec;

my $repo = abs_path(File::Spec->catdir($Bin, "..", ".."));
my $packet = File::Spec->catdir(
    $repo, "05-knowledge", "results",
    "lrc14_signature_response_congruence_thm4286"
);
my $inherited = File::Spec->catdir(
    $repo, "05-knowledge", "results",
    "lrc14_inactive_signature_deck_surgery_thm4282"
);
my $current = File::Spec->catdir(
    $repo, "05-knowledge", "results",
    "lrc14_endpoint_carrier_signature_surgery_thm4283",
    "results", "proof_graph"
);

my %inputs = (
    post4282 => [
        File::Spec->catfile(
            $inherited, "components", "surgery256", "results",
            "final_residual23373.csv"
        ),
        "c3e5bf37887aa57af79cb166fce4a6e933e5daffc26dd8032fdfc52ce31240f3",
        23_373,
    ],
    thm4283_union => [
        File::Spec->catfile($current, "proof_union.csv"),
        "c5646e81b3815bdef5168e36bcd76174065ed21339a5d8853d9efddc8fa3efae",
        691,
    ],
    post4283 => [
        File::Spec->catfile($current, "final_residual.csv"),
        "7d0044bc4c32f08b9d09dca420171a05666d26e03f38fbc48a9baa1fcb027102",
        22_682,
    ],
    index396 => [
        File::Spec->catfile($packet, "index396_fibre36.csv"),
        "114050c0100fe58793ed32f38f5c9e6bb530854bc74cc46765e4307a877b5fc6",
        36,
    ],
    two_mask => [
        File::Spec->catfile($packet, "two_mask_net36.csv"),
        "7a8b049abe73018e650420d5773fae630733ccfb47f7b5d775023affe23220cd",
        36,
    ],
);

sub slurp_raw {
    my ($path) = @_;
    open my $handle, "<:raw", $path or die "cannot open $path: $!\n";
    local $/;
    my $bytes = <$handle>;
    close $handle or die "cannot close $path: $!\n";
    return $bytes;
}

sub read_set {
    my ($name) = @_;
    my ($path, $expected_hash, $expected_count) = @{$inputs{$name}};
    my $bytes = slurp_raw($path);
    my $actual_hash = sha256_hex($bytes);
    die "$name SHA-256 changed: $actual_hash\n"
        unless $actual_hash eq $expected_hash;
    my %set;
    my $line_number = 0;
    for my $line (split /\n/, $bytes) {
        ++$line_number;
        die "bad row at $name:$line_number\n"
            unless $line =~ /^(\d+),(\d+)$/ && 0 < $1 && $1 < $2;
        die "duplicate row $line in $name\n" if exists $set{$line};
        $set{$line} = [$1 + 0, $2 + 0];
    }
    die "$name count changed\n" unless scalar(keys %set) == $expected_count;
    return \%set;
}

sub set_difference {
    my ($left, $right) = @_;
    return { map { $_ => $left->{$_} } grep { !exists $right->{$_} } keys %{$left} };
}

sub set_intersection {
    my ($left, $right) = @_;
    return { map { $_ => $left->{$_} } grep { exists $right->{$_} } keys %{$left} };
}

sub same_keys {
    my ($left, $right) = @_;
    return join("\n", sort keys %{$left}) eq join("\n", sort keys %{$right});
}

my $post4282 = read_set("post4282");
my $proof_union = read_set("thm4283_union");
my $post4283 = read_set("post4283");
my $index396 = read_set("index396");
my $two_mask = read_set("two_mask");

die "current THM-4283 set difference changed\n"
    unless same_keys(set_difference($post4282, $proof_union), $post4283);
die "current THM-4283 split overlaps\n"
    if scalar(keys %{set_intersection($proof_union, $post4283)});

my $index_overlap = set_intersection($index396, $proof_union);
my $two_overlap = set_intersection($two_mask, $proof_union);
my $index_novel = set_difference($index396, $proof_union);
my $two_novel = set_difference($two_mask, $proof_union);
my $branch_overlap = set_intersection($index396, $two_mask);
die "index-396 branch is not fully subsumed\n"
    unless scalar(keys %{$index_overlap}) == 36 && !scalar(keys %{$index_novel});
die "two-mask branch is not fully subsumed\n"
    unless scalar(keys %{$two_overlap}) == 36 && !scalar(keys %{$two_novel});
die "alternate proof nodes overlap\n" if scalar(keys %{$branch_overlap});

my %alternate_union = (%{$index396}, %{$two_mask});
die "alternate-node union count changed\n"
    unless scalar(keys %alternate_union) == 72;
my @alternate_rows = sort {
    $alternate_union{$a}[0] <=> $alternate_union{$b}[0] ||
    $alternate_union{$a}[1] <=> $alternate_union{$b}[1]
} keys %alternate_union;
my $alternate_bytes = join("", map { "$_\n" } @alternate_rows);
my $alternate_hash = sha256_hex($alternate_bytes);
die "alternate-node union identity changed\n"
    unless $alternate_hash eq
        "253071871ede4041c658ac7705de5283794f1baa230c9df1e22d16e22ac830b3";

my $maximum = 0;
for my $pair (values %{$post4283}) {
    $maximum = $pair->[1] if $pair->[1] > $maximum;
}
my @top = sort {
    $post4283->{$a}[0] <=> $post4283->{$b}[0]
} grep { $post4283->{$_}[1] == $maximum } keys %{$post4283};
die "current top changed\n"
    unless $maximum == 637 && join(" ", @top) eq "100,637 294,637 520,637";

print "THM4286_INDEPENDENT_PROOF_GRAPH_SUBSUMPTION_V1\n";
print "POST4282 23373\n";
print "CURRENT_THM4283_UNION 691\n";
print "CURRENT_RESIDUAL 22682\n";
print "INDEX396_BRANCH 36 OVERLAP 36 NOVEL 0\n";
print "TWO_MASK_BRANCH 36 OVERLAP 36 NOVEL 0\n";
print "ALTERNATE_NODE_UNION 72 SHA256=$alternate_hash\n";
print "TOP 637 @top\n";
print "VERDICT PASS INDEPENDENT EXACT SUBSUMPTION; ZERO CURRENT LEDGER DECREMENT\n";
