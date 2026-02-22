#!/usr/bin/env python3
"""Run ProteinMPNN via Gradio HTTP API multiple times and collect results."""

import json
import time
import csv
import os
import re
import requests
import random

BASE_URL = "https://simonduerr-proteinmpnn.hf.space"


def run_single(pdb_code, designed_chain="A,B,C", fixed_chain="", homomer=True,
               num_sequences=15, temperature="0.1"):
    """Run a single ProteinMPNN prediction via HTTP API."""
    payload = {
        "data": [
            pdb_code,
            None,
            designed_chain,
            fixed_chain,
            homomer,
            num_sequences,
            temperature,
            "vanilla\u2014v_48_020",  # em-dash
            "0",
            "",
            "",
        ],
        "fn_index": 0,
        "session_hash": f"s{random.randint(100000,999999)}",
    }

    resp = requests.post(f"{BASE_URL}/run/predict", json=payload, timeout=300)
    resp.raise_for_status()
    data = resp.json()
    return data["data"]


def parse_sequences(status_text, pdb_code, run_idx):
    """Parse sequences and scores from the FASTA-like output."""
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


def run_proteinmpnn(pdb_code, n_runs=10, **kwargs):
    """Run ProteinMPNN n_runs times and collect all sequences."""
    all_sequences = []

    for run_idx in range(n_runs):
        print(f"  Run {run_idx + 1}/{n_runs} for {pdb_code}...", end=" ", flush=True)
        try:
            result_data = run_single(pdb_code, **kwargs)
            status_text = result_data[0] if result_data else ""
            sequences = parse_sequences(status_text, pdb_code, run_idx)
            all_sequences.extend(sequences)
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

    return all_sequences


def save_results(sequences, filename):
    if not sequences:
        print(f"No sequences to save for {filename}")
        return
    with open(filename, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=['pdb', 'run', 'header', 'sequence', 'score', 'seq_recovery'])
        writer.writeheader()
        writer.writerows(sequences)
    print(f"Saved {len(sequences)} sequences to {filename}")


if __name__ == "__main__":
    output_dir = os.path.dirname(os.path.abspath(__file__))

    print("=" * 60)
    print("Running ProteinMPNN for 1TNF (10 runs, 15 seqs each)")
    print("=" * 60)
    seqs_1tnf = run_proteinmpnn("1tnf", n_runs=10)
    save_results(seqs_1tnf, os.path.join(output_dir, "1tnf_results.csv"))

    print("\n" + "=" * 60)
    print("Running ProteinMPNN for 6OOY (10 runs, 15 seqs each)")
    print("=" * 60)
    seqs_6ooy = run_proteinmpnn("6ooy", n_runs=10)
    save_results(seqs_6ooy, os.path.join(output_dir, "6ooy_results.csv"))

    print("\n" + "=" * 60)
    print(f"DONE! 1TNF: {len(seqs_1tnf)} seqs | 6OOY: {len(seqs_6ooy)} seqs")
    print("=" * 60)
