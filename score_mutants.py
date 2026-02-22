#!/usr/bin/env python3
"""Score single-point mutants on both 1TNF and 6OOY backbones using ProteinMPNN.

Uses conditional_probs_only mode to efficiently compute P(aa_i | backbone, WT_context)
for all 21 amino acids at every position in a single model pass per backbone.

Position matching is done via PDB residue numbers: we find the set of PDB resnums
present in both RCSB structures and map them to ProteinMPNN monomer indices.

Usage:
    python score_mutants.py
"""

import os
import csv
import json
import subprocess
import urllib.request
import numpy as np

WORK_DIR = os.path.dirname(os.path.abspath(__file__))
MPNN_DIR = os.path.join(WORK_DIR, 'ProteinMPNN')
PYTHON = 'python'  # assumes ProteinMPNN environment is active

ALPHABET = 'ACDEFGHIKLMNPQRSTVWYX'
AA20 = list('ACDEFGHIKLMNPQRSTVWY')  # 20 standard AAs (no X)
AA_TO_IDX = {aa: i for i, aa in enumerate(ALPHABET)}
AA3TO1 = {'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F',
          'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L',
          'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R',
          'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'}


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
    """Greedy match of MPNN monomer sequence to RCSB PDB resnums.
    Returns list of PDB resnums, one per monomer position."""
    mapping = []
    j = 0
    for i in range(len(pdb_seq)):
        if j < len(monomer_seq) and pdb_seq[i] == monomer_seq[j]:
            mapping.append(pdb_resnums[i])
            j += 1
    assert j == len(monomer_seq), f'Only matched {j}/{len(monomer_seq)} residues'
    return mapping


def build_shared_positions_by_resnum(wt_1tnf, wt_6ooy, resnums_1tnf, resnums_6ooy):
    """Build shared positions using PDB residue number matching.

    Returns list of dicts with pdb_resnum, mono_idx for each structure, and WT AAs.
    """
    # Build mono_idx -> pdb_resnum maps
    mono_to_rn_1 = {i: rn for i, rn in enumerate(resnums_1tnf)}
    mono_to_rn_2 = {i: rn for i, rn in enumerate(resnums_6ooy)}

    # Build pdb_resnum -> mono_idx maps
    rn_to_mono_1 = {rn: i for i, rn in enumerate(resnums_1tnf)}
    rn_to_mono_2 = {rn: i for i, rn in enumerate(resnums_6ooy)}

    # Find shared PDB resnums (present in both MPNN monomers)
    shared_rn = sorted(set(resnums_1tnf) & set(resnums_6ooy))

    shared = []
    for rn in shared_rn:
        mi1 = rn_to_mono_1[rn]
        mi2 = rn_to_mono_2[rn]
        shared.append({
            'pdb_resnum': rn,
            'mono_idx_1tnf': mi1,
            'mono_idx_6ooy': mi2,
            'wt_1tnf': wt_1tnf[mi1],
            'wt_6ooy': wt_6ooy[mi2],
        })
    return shared


def main():
    print('=' * 60)
    print('Single-Point Mutant Scoring: 1TNF vs 6OOY')
    print('(PDB resnum-matched, all chains independent)')
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

    print(f'  1TNF log_p shape: {data_1tnf["log_p"].shape} -> averaged {logp_1tnf.shape}')
    print(f'  6OOY log_p shape: {data_6ooy["log_p"].shape} -> averaged {logp_6ooy.shape}')

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

    # Step 5: Build PDB resnum-matched positions and score mutations
    print('\n[5/6] Building PDB resnum-matched positions and scoring...')
    print('  Fetching RCSB PDB structures for residue number mapping...')
    pdb_rn_1tnf, pdb_seq_1tnf = get_pdb_chain_a_resnums('1TNF')
    pdb_rn_6ooy, pdb_seq_6ooy = get_pdb_chain_a_resnums('6OOY')
    print(f'  1TNF PDB chain A: {len(pdb_rn_1tnf)} residues (resnums {pdb_rn_1tnf[0]}-{pdb_rn_1tnf[-1]})')
    print(f'  6OOY PDB chain A: {len(pdb_rn_6ooy)} residues (resnums {pdb_rn_6ooy[0]}-{pdb_rn_6ooy[-1]})')

    resnums_1tnf = map_monomer_to_pdb_resnums(wt_1tnf, pdb_rn_1tnf, pdb_seq_1tnf)
    resnums_6ooy = map_monomer_to_pdb_resnums(wt_6ooy, pdb_rn_6ooy, pdb_seq_6ooy)
    print(f'  1TNF monomer -> PDB resnums: {resnums_1tnf[0]}-{resnums_1tnf[-1]}')
    print(f'  6OOY monomer -> PDB resnums: {resnums_6ooy[0]}-{resnums_6ooy[-1]}')

    shared = build_shared_positions_by_resnum(wt_1tnf, wt_6ooy, resnums_1tnf, resnums_6ooy)
    n_same_wt = sum(1 for s in shared if s['wt_1tnf'] == s['wt_6ooy'])
    print(f'  {len(shared)} shared PDB resnums, {n_same_wt} with identical WT, {len(shared)-n_same_wt} with WT mismatch')
    for s in shared:
        if s['wt_1tnf'] != s['wt_6ooy']:
            print(f'    PDB resnum {s["pdb_resnum"]}: 1TNF={s["wt_1tnf"]}, 6OOY={s["wt_6ooy"]}')

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
                    'pdb_resnum': sp['pdb_resnum'],
                    'chain': chain,
                    'wt_1tnf': wt1,
                    'wt_6ooy': wt2,
                    'mut_aa': mut_aa,
                    'delta_1tnf': round(d1, 4),
                    'delta_6ooy': round(d2, 4),
                    'dd_score': round(dd, 4),
                    'abs_dd': round(abs(dd), 4),
                    'wt_same': wt1 == wt2,
                })

    results.sort(key=lambda r: r['abs_dd'], reverse=True)

    # Write full CSV
    out_path = os.path.join(WORK_DIR, 'mutation_ddg_scores.csv')
    with open(out_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(results[0].keys()))
        writer.writeheader()
        writer.writerows(results)
    print(f'\n  Full results: {out_path} ({len(results)} rows)')

    # Build per-position summary (max |DD| across all mutations, per chain)
    pos_summary = {}
    for r in results:
        key = (r['pdb_resnum'], r['chain'])
        if key not in pos_summary or r['abs_dd'] > pos_summary[key]['abs_dd']:
            pos_summary[key] = r

    # Build per-position summary across all chains (max |DD| at each position)
    pos_max = {}
    for r in results:
        rn = r['pdb_resnum']
        if rn not in pos_max or r['abs_dd'] > pos_max[rn]['abs_dd']:
            pos_max[rn] = r

    # Write JSON for visualization
    viz_data = []
    for sp in shared:
        rn = sp['pdb_resnum']
        entry = {
            'rn': rn,
            'wt1': sp['wt_1tnf'],
            'wt2': sp['wt_6ooy'],
        }
        for chain in chain_labels:
            key = (rn, chain)
            if key in pos_summary:
                r = pos_summary[key]
                entry[f'dd_{chain}'] = r['dd_score']
                entry[f'abs_dd_{chain}'] = r['abs_dd']
                entry[f'mut_{chain}'] = r['mut_aa']
                entry[f'd1_{chain}'] = r['delta_1tnf']
                entry[f'd2_{chain}'] = r['delta_6ooy']
        if rn in pos_max:
            entry['max_abs_dd'] = pos_max[rn]['abs_dd']
            entry['max_dd'] = pos_max[rn]['dd_score']
            entry['max_mut'] = pos_max[rn]['mut_aa']
            entry['max_chain'] = pos_max[rn]['chain']
        viz_data.append(entry)

    viz_path = os.path.join(WORK_DIR, 'mutation_viz_data.json')
    with open(viz_path, 'w') as f:
        json.dump(viz_data, f)
    print(f'  Viz data: {viz_path}')

    # Step 6: Print summary
    print(f'\n{"=" * 75}')
    print(f'TOP 20 MUTATIONS BY |DDscore|')
    print(f'{"=" * 75}')
    print(f'{"Res":>4} {"Ch":>2} {"WT1":>3} {"WT2":>3} {"Mut":>3} '
          f'{"D_1TNF":>8} {"D_6OOY":>8} {"DD":>8}')
    print('-' * 55)
    for r in results[:20]:
        note = '' if r['wt_same'] else '  WT-diff'
        print(f'{r["pdb_resnum"]:>4} {r["chain"]:>2} {r["wt_1tnf"]:>3} {r["wt_6ooy"]:>3} {r["mut_aa"]:>3} '
              f'{r["delta_1tnf"]:>8.3f} {r["delta_6ooy"]:>8.3f} '
              f'{r["dd_score"]:>8.3f}{note}')

    pos_ranked = sorted(pos_max.values(), key=lambda r: r['abs_dd'], reverse=True)
    same_wt = [r for r in pos_ranked if r['wt_same']]

    print(f'\n{"=" * 75}')
    print(f'TOP 20 POSITIONS BY MAX |DDscore| (same WT only)')
    print(f'{"=" * 75}')
    print(f'{"Res":>4} {"Ch":>2} {"WT":>3} {"Mut":>3} '
          f'{"D_1TNF":>8} {"D_6OOY":>8} {"DD":>8}')
    print('-' * 55)
    for r in same_wt[:20]:
        print(f'{r["pdb_resnum"]:>4} {r["chain"]:>2} {r["wt_1tnf"]:>3} {r["mut_aa"]:>3} '
              f'{r["delta_1tnf"]:>8.3f} {r["delta_6ooy"]:>8.3f} '
              f'{r["dd_score"]:>8.3f}')


if __name__ == '__main__':
    main()
