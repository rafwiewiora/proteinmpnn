#!/usr/bin/env python3
"""Compare ProteinMPNN results between 1TNF and 6OOY."""

import csv
import numpy as np
from collections import Counter

def load_sequences(filename):
    """Load sequences from CSV, return only designed (non-cleaned) sequences."""
    seqs = []
    with open(filename) as f:
        reader = csv.DictReader(f)
        for row in reader:
            # Skip the "cleaned" (original) sequence
            if 'cleaned' in row['header']:
                continue
            seqs.append(row)
    return seqs


def get_monomer_sequences(sequences):
    """Extract single-chain (monomer) sequences.

    The homomer output has repeated chains separated by '/'.
    We just take the first chain since they're identical in homomer mode.
    """
    monomer_seqs = []
    for s in sequences:
        full_seq = s['sequence']
        # Split by '/' to get individual chains
        chains = full_seq.split('/')
        monomer_seqs.append({
            **s,
            'monomer_seq': chains[0],  # First chain
            'n_chains': len(chains),
        })
    return monomer_seqs


def compute_position_entropy(sequences, seq_key='monomer_seq'):
    """Compute per-position amino acid frequency and entropy."""
    if not sequences:
        return [], []

    seq_len = len(sequences[0][seq_key])
    position_counts = [Counter() for _ in range(seq_len)]

    for s in sequences:
        seq = s[seq_key]
        for i, aa in enumerate(seq):
            position_counts[i][aa] += 1

    n_seqs = len(sequences)
    entropies = []
    for i, counts in enumerate(position_counts):
        entropy = 0
        for aa, count in counts.items():
            p = count / n_seqs
            if p > 0:
                entropy -= p * np.log2(p)
        entropies.append(entropy)

    return position_counts, entropies


def compute_consensus(position_counts):
    """Get consensus sequence from position counts."""
    consensus = []
    for counts in position_counts:
        if counts:
            consensus.append(counts.most_common(1)[0][0])
        else:
            consensus.append('-')
    return ''.join(consensus)


def simple_align(seq1, seq2):
    """Simple Needleman-Wunsch alignment for comparing two sequences."""
    # Since both are TNF-alpha with similar length, do simple NW
    from functools import lru_cache

    gap_penalty = -2
    match_score = 1
    mismatch_score = -1

    n, m = len(seq1), len(seq2)

    # DP table
    dp = np.zeros((n+1, m+1))
    traceback = np.zeros((n+1, m+1), dtype=int)  # 0=diag, 1=up, 2=left

    for i in range(1, n+1):
        dp[i][0] = i * gap_penalty
        traceback[i][0] = 1
    for j in range(1, m+1):
        dp[0][j] = j * gap_penalty
        traceback[0][j] = 2

    for i in range(1, n+1):
        for j in range(1, m+1):
            match = dp[i-1][j-1] + (match_score if seq1[i-1] == seq2[j-1] else mismatch_score)
            delete = dp[i-1][j] + gap_penalty
            insert = dp[0][j-1] + gap_penalty  # Bug: should be dp[i][j-1]

            # Fix: recalculate
            insert = dp[i][j-1] + gap_penalty

            if match >= delete and match >= insert:
                dp[i][j] = match
                traceback[i][j] = 0
            elif delete >= insert:
                dp[i][j] = delete
                traceback[i][j] = 1
            else:
                dp[i][j] = insert
                traceback[i][j] = 2

    # Traceback
    aligned1, aligned2 = [], []
    i, j = n, m
    while i > 0 or j > 0:
        if i > 0 and j > 0 and traceback[i][j] == 0:
            aligned1.append(seq1[i-1])
            aligned2.append(seq2[j-1])
            i -= 1
            j -= 1
        elif i > 0 and traceback[i][j] == 1:
            aligned1.append(seq1[i-1])
            aligned2.append('-')
            i -= 1
        else:
            aligned1.append('-')
            aligned2.append(seq2[j-1])
            j -= 1

    return ''.join(reversed(aligned1)), ''.join(reversed(aligned2))


def align_position_data(counts1, counts2, consensus1, consensus2):
    """Align two sets of position data using consensus sequence alignment."""
    aln1, aln2 = simple_align(consensus1, consensus2)

    aligned_counts1 = []
    aligned_counts2 = []

    idx1, idx2 = 0, 0
    for a, b in zip(aln1, aln2):
        if a != '-':
            c1 = counts1[idx1]
            idx1 += 1
        else:
            c1 = Counter()

        if b != '-':
            c2 = counts2[idx2]
            idx2 += 1
        else:
            c2 = Counter()

        aligned_counts1.append(c1)
        aligned_counts2.append(c2)

    return aligned_counts1, aligned_counts2, aln1, aln2


def compute_divergence(counts1, counts2, n_seqs1, n_seqs2):
    """Compute Jensen-Shannon divergence between two position distributions."""
    all_aas = set(list(counts1.keys()) + list(counts2.keys()))
    if not all_aas:
        return 0

    p = np.array([counts1.get(aa, 0) / max(n_seqs1, 1) for aa in sorted(all_aas)])
    q = np.array([counts2.get(aa, 0) / max(n_seqs2, 1) for aa in sorted(all_aas)])

    # JS divergence
    m = (p + q) / 2

    def kl(a, b):
        result = 0
        for ai, bi in zip(a, b):
            if ai > 0 and bi > 0:
                result += ai * np.log2(ai / bi)
        return result

    return (kl(p, m) + kl(q, m)) / 2


if __name__ == "__main__":
    print("=" * 70)
    print("COMPARISON: 1TNF vs 6OOY ProteinMPNN Designs")
    print("=" * 70)

    # Load data
    seqs_1tnf = load_sequences("/Users/rafal/repos/proteinmpnn/1tnf_results.csv")
    seqs_6ooy = load_sequences("/Users/rafal/repos/proteinmpnn/6ooy_homomer_results.csv")

    print(f"\n1TNF: {len(seqs_1tnf)} designed sequences")
    print(f"6OOY: {len(seqs_6ooy)} designed sequences")

    # Extract monomers
    mono_1tnf = get_monomer_sequences(seqs_1tnf)
    mono_6ooy = get_monomer_sequences(seqs_6ooy)

    print(f"\n1TNF monomer length: {len(mono_1tnf[0]['monomer_seq'])}")
    print(f"6OOY monomer length: {len(mono_6ooy[0]['monomer_seq'])}")

    # Score statistics
    scores_1tnf = [float(s['score']) for s in seqs_1tnf if s['score']]
    scores_6ooy = [float(s['score']) for s in seqs_6ooy if s['score']]
    print(f"\n1TNF scores: mean={np.mean(scores_1tnf):.3f} ± {np.std(scores_1tnf):.3f} (range {min(scores_1tnf):.3f}-{max(scores_1tnf):.3f})")
    print(f"6OOY scores: mean={np.mean(scores_6ooy):.3f} ± {np.std(scores_6ooy):.3f} (range {min(scores_6ooy):.3f}-{max(scores_6ooy):.3f})")

    # Compute per-position statistics
    counts_1tnf, entropy_1tnf = compute_position_entropy(mono_1tnf)
    counts_6ooy, entropy_6ooy = compute_position_entropy(mono_6ooy)

    consensus_1tnf = compute_consensus(counts_1tnf)
    consensus_6ooy = compute_consensus(counts_6ooy)

    print(f"\n1TNF consensus: {consensus_1tnf}")
    print(f"6OOY consensus: {consensus_6ooy}")

    # Align consensus sequences
    aln_counts1, aln_counts2, aln1, aln2 = align_position_data(
        counts_1tnf, counts_6ooy, consensus_1tnf, consensus_6ooy
    )

    print(f"\nAligned 1TNF: {aln1}")
    print(f"Aligned 6OOY: {aln2}")

    # Match line
    match_line = []
    for a, b in zip(aln1, aln2):
        if a == b:
            match_line.append('|')
        elif a == '-' or b == '-':
            match_line.append(' ')
        else:
            match_line.append('.')
    print(f"             {''.join(match_line)}")

    # Compute per-position divergence
    n1 = len(mono_1tnf)
    n2 = len(mono_6ooy)

    divergences = []
    for i, (c1, c2) in enumerate(zip(aln_counts1, aln_counts2)):
        div = compute_divergence(c1, c2, n1, n2)

        # Get most common AA for each
        top1 = c1.most_common(1)[0] if c1 else ('-', 0)
        top2 = c2.most_common(1)[0] if c2 else ('-', 0)

        divergences.append({
            'aligned_pos': i + 1,
            'aa_1tnf': aln1[i],
            'aa_6ooy': aln2[i],
            'top_aa_1tnf': top1[0],
            'freq_1tnf': top1[1] / max(n1, 1),
            'top_aa_6ooy': top2[0],
            'freq_6ooy': top2[1] / max(n2, 1),
            'js_divergence': div,
            'entropy_diff': abs(
                (-sum(p/n1 * np.log2(p/n1) for p in c1.values() if p > 0) if c1 else 0) -
                (-sum(p/n2 * np.log2(p/n2) for p in c2.values() if p > 0) if c2 else 0)
            ),
        })

    # Sort by divergence
    divergences.sort(key=lambda x: x['js_divergence'], reverse=True)

    print("\n" + "=" * 70)
    print("TOP 20 POSITIONS WITH HIGHEST DIVERGENCE (Jensen-Shannon)")
    print("=" * 70)
    print(f"{'Pos':>4} {'1TNF':>5} {'6OOY':>5} {'Top1TNF':>8} {'Freq1':>6} {'Top6OOY':>8} {'Freq2':>6} {'JS-Div':>8}")
    print("-" * 70)

    for d in divergences[:20]:
        print(f"{d['aligned_pos']:>4} {d['aa_1tnf']:>5} {d['aa_6ooy']:>5} "
              f"{d['top_aa_1tnf']:>8} {d['freq_1tnf']:>6.2f} "
              f"{d['top_aa_6ooy']:>8} {d['freq_6ooy']:>6.2f} "
              f"{d['js_divergence']:>8.4f}")

    # Overall summary
    consensus_matches = sum(1 for a, b in zip(aln1, aln2) if a == b and a != '-')
    consensus_mismatches = sum(1 for a, b in zip(aln1, aln2) if a != b and a != '-' and b != '-')
    print(f"\nConsensus alignment: {consensus_matches} matches, {consensus_mismatches} mismatches out of {len(aln1)} aligned positions")
    print(f"Consensus identity: {consensus_matches / (consensus_matches + consensus_mismatches) * 100:.1f}%")

    # Save full divergence data
    outfile = "/Users/rafal/repos/proteinmpnn/divergence_analysis.csv"
    with open(outfile, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(divergences[0].keys()))
        writer.writeheader()
        # Sort by position for the CSV
        writer.writerows(sorted(divergences, key=lambda x: x['aligned_pos']))
    print(f"\nFull divergence data saved to {outfile}")

    # Distribution of divergences
    js_values = [d['js_divergence'] for d in divergences]
    print(f"\nDivergence statistics:")
    print(f"  Mean JS-divergence: {np.mean(js_values):.4f}")
    print(f"  Median JS-divergence: {np.median(js_values):.4f}")
    print(f"  Max JS-divergence: {max(js_values):.4f}")
    print(f"  Positions with JS > 0.1: {sum(1 for v in js_values if v > 0.1)}")
    print(f"  Positions with JS > 0.2: {sum(1 for v in js_values if v > 0.2)}")
    print(f"  Positions with JS > 0.5: {sum(1 for v in js_values if v > 0.5)}")
