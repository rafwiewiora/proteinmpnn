#!/usr/bin/env python3
"""Score single-point mutants on both 1TNF and 6OOY backbones using ProteinMPNN.

Uses conditional_probs_only mode to efficiently compute P(aa_i | backbone, WT_context)
for all 21 amino acids at every position in a single model pass per backbone.

Every position across all 3 chains is treated independently (as if a single-chain
construct with linkers). DDscore computed for each chain×position×mutation.

Usage:
    /Users/rafal/miniconda3/envs/mpnn/bin/python score_mutants.py
"""

import os
import sys
import csv
import json
import subprocess
import urllib.request
import numpy as np

MPNN_DIR = '/Users/rafal/repos/proteinmpnn/ProteinMPNN'
WORK_DIR = '/Users/rafal/repos/proteinmpnn'
PYTHON = '/Users/rafal/miniconda3/envs/mpnn/bin/python'

ALPHABET = 'ACDEFGHIKLMNPQRSTVWYX'
AA20 = list('ACDEFGHIKLMNPQRSTVWY')  # 20 standard AAs (no X)
AA_TO_IDX = {aa: i for i, aa in enumerate(ALPHABET)}


def download_pdb(pdb_id):
    """Download PDB from RCSB if not present."""
    path = os.path.join(WORK_DIR, f'{pdb_id.lower()}.pdb')
    if os.path.exists(path):
        print(f'  {path} already exists')
        return path
    url = f'https://files.rcsb.org/download/{pdb_id}.pdb'
    print(f'  Downloading {url}...')
    urllib.request.urlretrieve(url, path)
    print(f'  Saved to {path}')
    return path


def run_conditional_probs(pdb_path, pdb_name, n_samples=4):
    """Run ProteinMPNN conditional_probs_only mode."""
    out_dir = os.path.join(WORK_DIR, f'cond_probs_{pdb_name}')
    npz_dir = os.path.join(out_dir, 'conditional_probs_only')

    if os.path.exists(npz_dir):
        npz_files = [f for f in os.listdir(npz_dir) if f.endswith('.npz')]
        if npz_files:
            npz_path = os.path.join(npz_dir, npz_files[0])
            print(f'  Loading cached results from {npz_path}')
            data = np.load(npz_path)
            if data['log_p'].shape[0] >= n_samples:
                return {k: data[k] for k in ['log_p', 'S', 'mask', 'design_mask']}

    os.makedirs(out_dir, exist_ok=True)
    cmd = [
        PYTHON, os.path.join(MPNN_DIR, 'protein_mpnn_run.py'),
        '--pdb_path', pdb_path,
        '--pdb_path_chains', 'A B C',
        '--out_folder', out_dir,
        '--conditional_probs_only', '1',
        '--num_seq_per_target', str(n_samples),
        '--batch_size', '1',
        '--seed', '42',
    ]

    print(f'  Running conditional_probs for {pdb_name} ({n_samples} samples)...')
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=1800)
    if result.stdout.strip():
        for line in result.stdout.strip().split('\n'):
            print(f'    {line}')
    if result.returncode != 0:
        print(f'  STDERR: {result.stderr}')
        raise RuntimeError(f'protein_mpnn_run.py failed for {pdb_name}')

    npz_files = [f for f in os.listdir(npz_dir) if f.endswith('.npz')]
    npz_path = os.path.join(npz_dir, npz_files[0])
    data = np.load(npz_path)
    return {k: data[k] for k in ['log_p', 'S', 'mask', 'design_mask']}


def get_wt_monomer(csv_path):
    with open(csv_path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            if 'cleaned' in row['header']:
                return row['sequence'].split('/')[0]
    return None


def decode_sequence(S_array):
    return ''.join(ALPHABET[int(s)] for s in S_array)


def build_chain_map(full_seq, wt_monomer):
    """Map monomer positions to tensor indices for all 3 chains.

    Returns: chain_len, list of (monomer_idx, [chainA_idx, chainB_idx, chainC_idx])
    for non-X positions only.
    """
    total_L = len(full_seq)
    assert total_L % 3 == 0
    chain_len = total_L // 3
    chain_a = full_seq[:chain_len]

    positions = []
    mono_pos = 0
    for i, aa in enumerate(chain_a):
        if aa == 'X':
            continue
        assert aa == wt_monomer[mono_pos], \
            f'Chain A[{i}]={aa} != monomer[{mono_pos}]={wt_monomer[mono_pos]}'
        positions.append((mono_pos, [i, i + chain_len, i + 2 * chain_len]))
        mono_pos += 1

    assert mono_pos == len(wt_monomer)
    return chain_len, positions


def load_divergence_data():
    csv_path = os.path.join(WORK_DIR, 'divergence_analysis.csv')
    data = []
    with open(csv_path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            data.append(row)
    data.sort(key=lambda r: int(r['aligned_pos']))
    return data


def build_shared_positions(div_data, wt_1tnf, wt_6ooy):
    """Build list of shared aligned positions with monomer indices."""
    idx_1tnf = 0
    idx_6ooy = 0
    shared = []
    for row in div_data:
        has_1tnf = row['aa_1tnf'] != '-'
        has_6ooy = row['aa_6ooy'] != '-'
        mi1 = idx_1tnf if has_1tnf else None
        mi2 = idx_6ooy if has_6ooy else None
        if has_1tnf:
            idx_1tnf += 1
        if has_6ooy:
            idx_6ooy += 1
        if mi1 is not None and mi2 is not None:
            shared.append({
                'aligned_pos': int(row['aligned_pos']),
                'mono_idx_1tnf': mi1,
                'mono_idx_6ooy': mi2,
                'wt_1tnf': wt_1tnf[mi1],
                'wt_6ooy': wt_6ooy[mi2],
                'js_divergence': float(row['js_divergence']),
            })
    return shared


def main():
    print('=' * 60)
    print('Single-Point Mutant Scoring: 1TNF vs 6OOY')
    print('(All positions independent, no chain averaging)')
    print('=' * 60)

    # Step 1: PDB files
    print('\n[1/6] Preparing PDB files...')
    pdb_1tnf = download_pdb('1TNF')
    pdb_6ooy = os.path.join(WORK_DIR, '6ooy_clean.pdb')
    assert os.path.exists(pdb_6ooy)
    print(f'  Using cleaned 6OOY: {pdb_6ooy}')

    # Step 2: WT sequences
    print('\n[2/6] Loading WT sequences...')
    wt_1tnf = get_wt_monomer(os.path.join(WORK_DIR, '1tnf_results.csv'))
    wt_6ooy = get_wt_monomer(os.path.join(WORK_DIR, '6ooy_homomer_results.csv'))
    L1, L2 = len(wt_1tnf), len(wt_6ooy)
    print(f'  1TNF WT monomer: {L1} residues')
    print(f'  6OOY WT monomer: {L2} residues')

    # Step 3: Run conditional probs
    print('\n[3/6] Computing conditional probabilities...')
    print('  --- 1TNF ---')
    data_1tnf = run_conditional_probs(pdb_1tnf, '1tnf')
    print('  --- 6OOY ---')
    data_6ooy = run_conditional_probs(pdb_6ooy, '6ooy_clean')

    logp_1tnf = data_1tnf['log_p'].mean(axis=0)
    logp_6ooy = data_6ooy['log_p'].mean(axis=0)
    S_1tnf = data_1tnf['S']
    S_6ooy = data_6ooy['S']

    print(f'  1TNF: {logp_1tnf.shape}')
    print(f'  6OOY: {logp_6ooy.shape}')

    # Step 4: Build chain mappings
    print('\n[4/6] Building chain mappings...')
    full_seq_1tnf = decode_sequence(S_1tnf)
    full_seq_6ooy = decode_sequence(S_6ooy)
    clen1, positions_1tnf = build_chain_map(full_seq_1tnf, wt_1tnf)
    clen2, positions_6ooy = build_chain_map(full_seq_6ooy, wt_6ooy)
    print(f'  1TNF chain len={clen1}, {len(positions_1tnf)} resolved positions/chain')
    print(f'  6OOY chain len={clen2}, {len(positions_6ooy)} resolved positions/chain')

    # Build mono_idx -> tensor indices lookup
    map1 = {mi: tis for mi, tis in positions_1tnf}
    map2 = {mi: tis for mi, tis in positions_6ooy}

    # Step 5: Build alignment and score mutations
    print('\n[5/6] Building alignment and scoring...')
    div_data = load_divergence_data()
    shared = build_shared_positions(div_data, wt_1tnf, wt_6ooy)
    print(f'  {len(shared)} shared aligned positions × 3 chains = {len(shared)*3} independent positions')

    chain_labels = ['A', 'B', 'C']
    results = []
    for sp in shared:
        mi1 = sp['mono_idx_1tnf']
        mi2 = sp['mono_idx_6ooy']
        wt1 = sp['wt_1tnf']
        wt2 = sp['wt_6ooy']
        wt1_idx = AA_TO_IDX[wt1]
        wt2_idx = AA_TO_IDX[wt2]
        tis1 = map1[mi1]
        tis2 = map2[mi2]

        for ci, chain in enumerate(chain_labels):
            ti1 = tis1[ci]
            ti2 = tis2[ci]

            for mut_aa in AA20:
                if mut_aa == wt1 and mut_aa == wt2:
                    continue

                mut_idx = AA_TO_IDX[mut_aa]
                d1 = float(logp_1tnf[ti1, wt1_idx] - logp_1tnf[ti1, mut_idx])
                d2 = float(logp_6ooy[ti2, wt2_idx] - logp_6ooy[ti2, mut_idx])
                dd = d1 - d2

                results.append({
                    'aligned_pos': sp['aligned_pos'],
                    'chain': chain,
                    'wt_1tnf': wt1,
                    'wt_6ooy': wt2,
                    'mut_aa': mut_aa,
                    'delta_1tnf': round(d1, 4),
                    'delta_6ooy': round(d2, 4),
                    'dd_score': round(dd, 4),
                    'abs_dd': round(abs(dd), 4),
                    'wt_same': wt1 == wt2,
                    'js_divergence': sp['js_divergence'],
                })

    results.sort(key=lambda r: r['abs_dd'], reverse=True)

    # Write full CSV
    out_path = os.path.join(WORK_DIR, 'mutation_ddg_scores.csv')
    with open(out_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(results[0].keys()))
        writer.writeheader()
        writer.writerows(results)
    print(f'  Full results: {out_path} ({len(results)} rows)')

    # Build per-position summary (max |DD| across all mutations, per chain)
    pos_summary = {}
    for r in results:
        key = (r['aligned_pos'], r['chain'])
        if key not in pos_summary or r['abs_dd'] > pos_summary[key]['abs_dd']:
            pos_summary[key] = r

    # Build per-position summary across all chains (max |DD| at each aligned pos)
    pos_max = {}
    for r in results:
        ap = r['aligned_pos']
        if ap not in pos_max or r['abs_dd'] > pos_max[ap]['abs_dd']:
            pos_max[ap] = r

    # Write JSON for visualization
    # Per-position data: for each aligned pos, max |DD| and best mutation per chain
    viz_data = []
    for sp in shared:
        ap = sp['aligned_pos']
        entry = {
            'ap': ap,
            'wt1': sp['wt_1tnf'],
            'wt2': sp['wt_6ooy'],
            'jsd': sp['js_divergence'],
        }
        for chain in chain_labels:
            key = (ap, chain)
            if key in pos_summary:
                r = pos_summary[key]
                entry[f'dd_{chain}'] = r['dd_score']
                entry[f'abs_dd_{chain}'] = r['abs_dd']
                entry[f'mut_{chain}'] = r['mut_aa']
                entry[f'd1_{chain}'] = r['delta_1tnf']
                entry[f'd2_{chain}'] = r['delta_6ooy']
        if ap in pos_max:
            entry['max_abs_dd'] = pos_max[ap]['abs_dd']
            entry['max_dd'] = pos_max[ap]['dd_score']
            entry['max_mut'] = pos_max[ap]['mut_aa']
            entry['max_chain'] = pos_max[ap]['chain']
        viz_data.append(entry)

    viz_path = os.path.join(WORK_DIR, 'mutation_viz_data.json')
    with open(viz_path, 'w') as f:
        json.dump(viz_data, f)
    print(f'  Viz data: {viz_path}')

    # Also save full per-position per-mutation data for heatmap
    heatmap_data = {}
    for r in results:
        ap = r['aligned_pos']
        chain = r['chain']
        if ap not in heatmap_data:
            heatmap_data[ap] = {'wt1': r['wt_1tnf'], 'wt2': r['wt_6ooy'], 'jsd': r['js_divergence'], 'muts': {}}
        key = f"{chain}_{r['mut_aa']}"
        heatmap_data[ap]['muts'][key] = {
            'd1': r['delta_1tnf'],
            'd2': r['delta_6ooy'],
            'dd': r['dd_score'],
        }

    heatmap_path = os.path.join(WORK_DIR, 'mutation_heatmap_data.json')
    with open(heatmap_path, 'w') as f:
        json.dump(heatmap_data, f)
    print(f'  Heatmap data: {heatmap_path}')

    # Print top 20
    print(f'\n{"=" * 75}')
    print(f'TOP 20 MUTATIONS BY |DDscore|')
    print(f'{"=" * 75}')
    print(f'{"Pos":>4} {"Ch":>2} {"WT1":>3} {"WT2":>3} {"Mut":>3} '
          f'{"D_1TNF":>8} {"D_6OOY":>8} {"DD":>8} {"JSD":>6}')
    print('-' * 75)
    for r in results[:20]:
        print(f'{r["aligned_pos"]:>4} {r["chain"]:>2} {r["wt_1tnf"]:>3} {r["wt_6ooy"]:>3} {r["mut_aa"]:>3} '
              f'{r["delta_1tnf"]:>8.3f} {r["delta_6ooy"]:>8.3f} '
              f'{r["dd_score"]:>8.3f} {r["js_divergence"]:>6.3f}')

    # Position-level summary
    pos_ranked = sorted(pos_max.values(), key=lambda r: r['abs_dd'], reverse=True)
    same_wt = [r for r in pos_ranked if r['wt_same']]

    print(f'\n{"=" * 75}')
    print(f'TOP 20 POSITIONS BY MAX |DDscore| (same WT only)')
    print(f'{"=" * 75}')
    print(f'{"Pos":>4} {"Ch":>2} {"WT":>3} {"Mut":>3} '
          f'{"D_1TNF":>8} {"D_6OOY":>8} {"DD":>8} {"JSD":>6}')
    print('-' * 60)
    for r in same_wt[:20]:
        print(f'{r["aligned_pos"]:>4} {r["chain"]:>2} {r["wt_1tnf"]:>3} {r["mut_aa"]:>3} '
              f'{r["delta_1tnf"]:>8.3f} {r["delta_6ooy"]:>8.3f} '
              f'{r["dd_score"]:>8.3f} {r["js_divergence"]:>6.3f}')

    # Chain asymmetry check
    print(f'\n{"=" * 75}')
    print('CHAIN ASYMMETRY CHECK (how different are A/B/C at same position)')
    print(f'{"=" * 75}')
    asym_scores = []
    for sp in shared:
        ap = sp['aligned_pos']
        for mut_aa in AA20:
            if mut_aa == sp['wt_1tnf'] and mut_aa == sp['wt_6ooy']:
                continue
            dds = []
            for chain in chain_labels:
                key = (ap, chain)
                # Find this specific mutation
                for r in results:
                    if r['aligned_pos'] == ap and r['chain'] == chain and r['mut_aa'] == mut_aa:
                        dds.append(r['dd_score'])
                        break
            if len(dds) == 3:
                spread = max(dds) - min(dds)
                if spread > 2.0:
                    asym_scores.append((ap, mut_aa, dds, spread))

    asym_scores.sort(key=lambda x: x[3], reverse=True)
    print(f'Positions with chain DD spread > 2.0: {len(asym_scores)}')
    print(f'{"Pos":>4} {"Mut":>3} {"DD_A":>8} {"DD_B":>8} {"DD_C":>8} {"Spread":>8}')
    for ap, mut, dds, spread in asym_scores[:10]:
        print(f'{ap:>4} {mut:>3} {dds[0]:>8.3f} {dds[1]:>8.3f} {dds[2]:>8.3f} {spread:>8.3f}')


if __name__ == '__main__':
    main()
