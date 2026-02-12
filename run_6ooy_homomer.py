#!/usr/bin/env python3
"""Run ProteinMPNN for 6OOY with cleaned PDB + homomer mode."""

import json
import time
import csv
import os
import re
import requests
import random

BASE_URL = "https://simonduerr-proteinmpnn.hf.space"
PDB_FILE = "/Users/rafal/repos/proteinmpnn/6ooy_clean.pdb"


def upload_file():
    """Upload PDB file to Gradio server, return file reference."""
    with open(PDB_FILE, 'rb') as f:
        resp = requests.post(f"{BASE_URL}/upload",
                             files=[('files', ('6ooy_clean.pdb', f, 'chemical/x-pdb'))])
    resp.raise_for_status()
    refs = resp.json()
    return refs[0] if isinstance(refs, list) else refs


def run_single(file_ref, num_sequences=15):
    """Run a single prediction with the uploaded file."""
    file_data = {
        'name': file_ref,
        'data': None,
        'is_file': True,
        'orig_name': '6ooy_clean.pdb'
    }
    payload = {
        "data": ["", file_data, "A,B,C", "", True, num_sequences, "0.1",
                 "vanilla\u2014v_48_020", "0", "", ""],
        "fn_index": 0,
        "session_hash": f"s{random.randint(100000,999999)}",
    }
    resp = requests.post(f"{BASE_URL}/run/predict", json=payload, timeout=300)
    resp.raise_for_status()
    return resp.json()["data"]


def parse_sequences(status_text, pdb_code, run_idx):
    sequences = []
    lines = status_text.strip().split('\n')
    current_header = None
    current_seq = []

    for line in lines:
        line = line.strip()
        if line.startswith('>'):
            if current_header and current_seq:
                seq_str = ''.join(current_seq)
                sequences.append({
                    'pdb': pdb_code,
                    'run': run_idx + 1,
                    'header': current_header,
                    'sequence': seq_str,
                    'score': extract_score(current_header),
                    'seq_recovery': extract_recovery(current_header),
                })
            current_header = line
            current_seq = []
        elif line and not line.startswith('#'):
            current_seq.append(line)

    if current_header and current_seq:
        seq_str = ''.join(current_seq)
        sequences.append({
            'pdb': pdb_code,
            'run': run_idx + 1,
            'header': current_header,
            'sequence': seq_str,
            'score': extract_score(current_header),
            'seq_recovery': extract_recovery(current_header),
        })
    return sequences


def extract_score(header):
    m = re.search(r'score=([0-9.]+)', header)
    return float(m.group(1)) if m else None


def extract_recovery(header):
    m = re.search(r'seq_recovery=([0-9.]+)', header)
    return float(m.group(1)) if m else None


if __name__ == "__main__":
    n_runs = 10
    all_seqs = []

    print("=" * 60)
    print("Running ProteinMPNN for 6OOY (cleaned, homomer=True)")
    print(f"  PDB: {PDB_FILE}")
    print(f"  Runs: {n_runs}, Sequences per run: 15")
    print("=" * 60)

    for run_idx in range(n_runs):
        print(f"  Run {run_idx + 1}/{n_runs}...", end=" ", flush=True)
        try:
            file_ref = upload_file()
            result_data = run_single(file_ref)
            status_text = result_data[0] if result_data else ""
            sequences = parse_sequences(status_text, "6ooy", run_idx)
            all_seqs.extend(sequences)
            scores = [s['score'] for s in sequences if s['score'] is not None]
            if scores:
                print(f"{len(sequences)} seqs, scores: {min(scores):.3f}-{max(scores):.3f}")
            else:
                print(f"{len(sequences)} seqs")
        except Exception as e:
            print(f"ERROR: {e}")
            time.sleep(5)
            continue
        time.sleep(1)

    # Save
    outfile = "/Users/rafal/repos/proteinmpnn/6ooy_homomer_results.csv"
    with open(outfile, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=['pdb', 'run', 'header', 'sequence', 'score', 'seq_recovery'])
        writer.writeheader()
        writer.writerows(all_seqs)

    print(f"\nSaved {len(all_seqs)} sequences to {outfile}")
