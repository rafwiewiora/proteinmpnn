#!/usr/bin/env python3
"""Compare ProteinMPNN results between 1TNF and 6OOY.

Position matching is done via PDB residue numbers (same approach as score_mutants.py).
This avoids the error-prone NW alignment of designed consensus sequences.
"""

import csv
import os
import urllib.request
import numpy as np
from collections import Counter

AA3TO1 = {'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F',
          'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L',
          'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R',
          'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'}

WORK_DIR = os.path.dirname(os.path.abspath(__file__))


def load_sequences(filename):
    """Load sequences from CSV, return only designed (non-cleaned) sequences."""
    seqs = []
    with open(filename) as f:
        reader = csv.DictReader(f)
        for row in reader:
            if 'cleaned' in row['header']:
                continue
            seqs.append(row)
    return seqs


def get_wt_monomer(csv_path):
    """Extract the WT (cleaned) monomer sequence."""
    with open(csv_path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            if 'cleaned' in row['header']:
                return row['sequence'].split('/')[0]
    return None


def get_monomer_sequences(sequences):
    """Extract single-chain (monomer) sequences."""
    monomer_seqs = []
    for s in sequences:
        full_seq = s['sequence']
        chains = full_seq.split('/')
        monomer_seqs.append({
            **s,
            'monomer_seq': chains[0],
            'n_chains': len(chains),
        })
    return monomer_seqs


def get_pdb_chain_a_resnums(pdb_id):
    """Get CA resnums and 1-letter AAs for chain A from RCSB PDB."""
    url = f'https://files.rcsb.org/download/{pdb_id}.pdb'
    text = urllib.request.urlopen(url).read().decode()
    resnums, seq = [], []
    seen = set()
    for line in text.split('\n'):
        if line.startswith('ATOM') and line[12:16].strip() == 'CA' and line[21] == 'A':
            rn = int(line[22:26].strip())
            if rn not in seen:
                seen.add(rn)
                resnums.append(rn)
                seq.append(AA3TO1.get(line[17:20].strip(), 'X'))
    return resnums, ''.join(seq)


def map_monomer_to_pdb_resnums(monomer_seq, pdb_resnums, pdb_seq):
    """Greedy match of MPNN monomer sequence to RCSB PDB resnums."""
    mapping = []
    j = 0
    for i in range(len(pdb_seq)):
        if j < len(monomer_seq) and pdb_seq[i] == monomer_seq[j]:
            mapping.append(pdb_resnums[i])
            j += 1
    assert j == len(monomer_seq), f'Only matched {j}/{len(monomer_seq)} residues'
    return mapping


def compute_position_counts(monomer_seqs, seq_key='monomer_seq'):
    """Compute per-position amino acid frequency counts."""
    if not monomer_seqs:
        return []
    seq_len = len(monomer_seqs[0][seq_key])
    position_counts = [Counter() for _ in range(seq_len)]
    for s in monomer_seqs:
        seq = s[seq_key]
        for i, aa in enumerate(seq):
            position_counts[i][aa] += 1
    return position_counts


def compute_jsd(counts1, counts2, n1, n2):
    """Compute Jensen-Shannon divergence between two AA distributions."""
    all_aas = set(list(counts1.keys()) + list(counts2.keys()))
    if not all_aas:
        return 0.0

    p = np.array([counts1.get(aa, 0) / max(n1, 1) for aa in sorted(all_aas)])
    q = np.array([counts2.get(aa, 0) / max(n2, 1) for aa in sorted(all_aas)])
    m = (p + q) / 2

    def kl(a, b):
        result = 0.0
        for ai, bi in zip(a, b):
            if ai > 0 and bi > 0:
                result += ai * np.log2(ai / bi)
        return result

    return (kl(p, m) + kl(q, m)) / 2


if __name__ == "__main__":
    print("=" * 70)
    print("COMPARISON: 1TNF vs 6OOY ProteinMPNN Designs")
    print("(PDB resnum-based position matching)")
    print("=" * 70)

    # Load designed sequences
    seqs_1tnf = load_sequences(f"{WORK_DIR}/1tnf_results.csv")
    seqs_6ooy = load_sequences(f"{WORK_DIR}/6ooy_homomer_results.csv")
    print(f"\n1TNF: {len(seqs_1tnf)} designed sequences")
    print(f"6OOY: {len(seqs_6ooy)} designed sequences")

    # WT monomers
    wt_1tnf = get_wt_monomer(f"{WORK_DIR}/1tnf_results.csv")
    wt_6ooy = get_wt_monomer(f"{WORK_DIR}/6ooy_homomer_results.csv")
    print(f"1TNF WT monomer: {len(wt_1tnf)} residues")
    print(f"6OOY WT monomer: {len(wt_6ooy)} residues")

    # Extract monomer sequences
    mono_1tnf = get_monomer_sequences(seqs_1tnf)
    mono_6ooy = get_monomer_sequences(seqs_6ooy)

    # Score statistics
    scores_1tnf = [float(s['score']) for s in seqs_1tnf if s['score']]
    scores_6ooy = [float(s['score']) for s in seqs_6ooy if s['score']]
    print(f"\n1TNF scores: mean={np.mean(scores_1tnf):.3f} +/- {np.std(scores_1tnf):.3f}")
    print(f"6OOY scores: mean={np.mean(scores_6ooy):.3f} +/- {np.std(scores_6ooy):.3f}")

    # Per-position AA counts
    counts_1tnf = compute_position_counts(mono_1tnf)
    counts_6ooy = compute_position_counts(mono_6ooy)

    # Consensus sequences
    cons_1tnf = ''.join(c.most_common(1)[0][0] if c else '-' for c in counts_1tnf)
    cons_6ooy = ''.join(c.most_common(1)[0][0] if c else '-' for c in counts_6ooy)

    # Map monomer positions to PDB resnums
    print("\nFetching PDB structures for residue number mapping...")
    pdb_rn_1tnf, pdb_seq_1tnf = get_pdb_chain_a_resnums('1TNF')
    pdb_rn_6ooy, pdb_seq_6ooy = get_pdb_chain_a_resnums('6OOY')
    print(f"1TNF PDB chain A: {len(pdb_rn_1tnf)} residues (PDB {pdb_rn_1tnf[0]}-{pdb_rn_1tnf[-1]})")
    print(f"6OOY PDB chain A: {len(pdb_rn_6ooy)} residues (PDB {pdb_rn_6ooy[0]}-{pdb_rn_6ooy[-1]})")

    resnums_1tnf = map_monomer_to_pdb_resnums(wt_1tnf, pdb_rn_1tnf, pdb_seq_1tnf)
    resnums_6ooy = map_monomer_to_pdb_resnums(wt_6ooy, pdb_rn_6ooy, pdb_seq_6ooy)
    print(f"1TNF monomer -> PDB resnums: {resnums_1tnf[0]}-{resnums_1tnf[-1]} ({len(resnums_1tnf)} pos)")
    print(f"6OOY monomer -> PDB resnums: {resnums_6ooy[0]}-{resnums_6ooy[-1]} ({len(resnums_6ooy)} pos)")

    # Build mono_idx -> PDB resnum maps
    rn_to_mi_1 = {rn: i for i, rn in enumerate(resnums_1tnf)}
    rn_to_mi_2 = {rn: i for i, rn in enumerate(resnums_6ooy)}
    set_rn_1 = set(resnums_1tnf)
    set_rn_2 = set(resnums_6ooy)
    shared_rn = sorted(set_rn_1 & set_rn_2)
    only_1tnf = sorted(set_rn_1 - set_rn_2)
    only_6ooy = sorted(set_rn_2 - set_rn_1)
    all_rn = sorted(set_rn_1 | set_rn_2)

    print(f"\nShared PDB resnums: {len(shared_rn)}")
    print(f"Only in 1TNF: {len(only_1tnf)} ({only_1tnf})")
    print(f"Only in 6OOY: {len(only_6ooy)} ({only_6ooy})")

    # Verify WT matching
    n_wt_match = 0
    n_wt_mismatch = 0
    for rn in shared_rn:
        wt1 = wt_1tnf[rn_to_mi_1[rn]]
        wt2 = wt_6ooy[rn_to_mi_2[rn]]
        if wt1 == wt2:
            n_wt_match += 1
        else:
            n_wt_mismatch += 1
            print(f"  WT mismatch at PDB resnum {rn}: 1TNF={wt1}, 6OOY={wt2}")
    print(f"WT identical: {n_wt_match}, WT mismatch: {n_wt_mismatch}")

    # Build divergence data for all positions (using PDB resnum as key)
    n1 = len(mono_1tnf)
    n2 = len(mono_6ooy)

    divergences = []
    for rn in all_rn:
        has_1tnf = rn in rn_to_mi_1
        has_6ooy = rn in rn_to_mi_2

        if has_1tnf:
            mi1 = rn_to_mi_1[rn]
            c1 = counts_1tnf[mi1]
            top1 = c1.most_common(1)[0] if c1 else ('-', 0)
            wt1_aa = wt_1tnf[mi1]
            cons1_aa = cons_1tnf[mi1]
        else:
            c1 = Counter()
            top1 = ('-', 0)
            wt1_aa = '-'
            cons1_aa = '-'

        if has_6ooy:
            mi2 = rn_to_mi_2[rn]
            c2 = counts_6ooy[mi2]
            top2 = c2.most_common(1)[0] if c2 else ('-', 0)
            wt2_aa = wt_6ooy[mi2]
            cons2_aa = cons_6ooy[mi2]
        else:
            c2 = Counter()
            top2 = ('-', 0)
            wt2_aa = '-'
            cons2_aa = '-'

        jsd = compute_jsd(c1, c2, n1, n2) if has_1tnf and has_6ooy else 0.5

        divergences.append({
            'pdb_resnum': rn,
            'wt_1tnf': wt1_aa,
            'wt_6ooy': wt2_aa,
            'cons_1tnf': cons1_aa,
            'cons_6ooy': cons2_aa,
            'top_aa_1tnf': top1[0],
            'freq_1tnf': top1[1] / max(n1, 1),
            'top_aa_6ooy': top2[0],
            'freq_6ooy': top2[1] / max(n2, 1),
            'js_divergence': jsd,
            'in_both': has_1tnf and has_6ooy,
        })

    # Print top divergent
    by_jsd = sorted(divergences, key=lambda x: x['js_divergence'], reverse=True)

    print(f"\n{'=' * 70}")
    print("TOP 20 POSITIONS WITH HIGHEST DIVERGENCE (Jensen-Shannon)")
    print(f"{'=' * 70}")
    print(f"{'PDB#':>5} {'WT1':>4} {'WT2':>4} {'Cons1':>6} {'Cons2':>6} "
          f"{'Top1':>5} {'Freq1':>6} {'Top2':>5} {'Freq2':>6} {'JSD':>8}")
    print("-" * 70)

    for d in by_jsd[:20]:
        both = '' if d['in_both'] else '  (gap)'
        print(f"{d['pdb_resnum']:>5} {d['wt_1tnf']:>4} {d['wt_6ooy']:>4} "
              f"{d['cons_1tnf']:>6} {d['cons_6ooy']:>6} "
              f"{d['top_aa_1tnf']:>5} {d['freq_1tnf']:>6.2f} "
              f"{d['top_aa_6ooy']:>5} {d['freq_6ooy']:>6.2f} "
              f"{d['js_divergence']:>8.4f}{both}")

    # Save CSV
    outfile = f"{WORK_DIR}/divergence_analysis.csv"
    with open(outfile, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(divergences[0].keys()))
        writer.writeheader()
        writer.writerows(sorted(divergences, key=lambda x: x['pdb_resnum']))
    print(f"\nSaved to {outfile} ({len(divergences)} rows)")

    # Statistics
    shared_jsd = [d['js_divergence'] for d in divergences if d['in_both']]
    print(f"\nDivergence statistics (shared positions only):")
    print(f"  Mean JSD: {np.mean(shared_jsd):.4f}")
    print(f"  Median JSD: {np.median(shared_jsd):.4f}")
    print(f"  Max JSD: {max(shared_jsd):.4f}")
    print(f"  Positions with JSD > 0.1: {sum(1 for v in shared_jsd if v > 0.1)}")
    print(f"  Positions with JSD > 0.2: {sum(1 for v in shared_jsd if v > 0.2)}")
    print(f"  Positions with JSD > 0.5: {sum(1 for v in shared_jsd if v > 0.5)}")
