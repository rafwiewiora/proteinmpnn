#!/usr/bin/env python3
"""Generate interactive HTML visualization of ProteinMPNN divergence between 1TNF and 6OOY.

Creates an HTML file with:
- Side-by-side Mol* viewers colored by JS-divergence (B-factor replacement + uncertainty theme)
- Aligned consensus sequences with color-coded positions
- Table of top divergent positions (clickable to focus viewers)
"""

import csv
import json
import urllib.request

AA3TO1 = {'ALA':'A','CYS':'C','ASP':'D','GLU':'E','PHE':'F','GLY':'G','HIS':'H',
          'ILE':'I','LYS':'K','LEU':'L','MET':'M','ASN':'N','PRO':'P','GLN':'Q',
          'ARG':'R','SER':'S','THR':'T','VAL':'V','TRP':'W','TYR':'Y'}


def get_wt_monomer(csv_path):
    """Extract the WT (cleaned) monomer sequence from a ProteinMPNN results CSV."""
    with open(csv_path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            if 'cleaned' in row['header']:
                return row['sequence'].split('/')[0]
    return None


def get_pdb_chain_a(pdb_id):
    """Fetch PDB and extract chain A residue numbers and 1-letter sequence."""
    url = f'https://files.rcsb.org/download/{pdb_id}.pdb'
    resp = urllib.request.urlopen(url)
    text = resp.read().decode()
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


def map_cleaned_to_pdb_resnums(cleaned_seq, pdb_resnums, pdb_seq):
    """Map each position in ProteinMPNN's cleaned sequence to its PDB residue number.

    The cleaned sequence is a subsequence of the PDB sequence (ProteinMPNN removes
    residues with incomplete backbone atoms). This finds the correct mapping by
    matching characters in order.
    """
    mapping = []
    j = 0
    for i in range(len(pdb_seq)):
        if j < len(cleaned_seq) and pdb_seq[i] == cleaned_seq[j]:
            mapping.append(pdb_resnums[i])
            j += 1
    assert j == len(cleaned_seq), f'Only matched {j}/{len(cleaned_seq)} residues'
    return mapping


def read_data(filepath):
    """Read divergence analysis CSV and return sorted list of dicts."""
    data = []
    with open(filepath) as f:
        reader = csv.DictReader(f)
        for row in reader:
            data.append({
                'aligned_pos': int(row['aligned_pos']),
                'aa_1tnf': row['aa_1tnf'],
                'aa_6ooy': row['aa_6ooy'],
                'top_aa_1tnf': row['top_aa_1tnf'],
                'freq_1tnf': float(row['freq_1tnf']),
                'top_aa_6ooy': row['top_aa_6ooy'],
                'freq_6ooy': float(row['freq_6ooy']),
                'js_divergence': float(row['js_divergence']),
            })
    data.sort(key=lambda x: x['aligned_pos'])
    return data


def compute_mappings(data, wt1_seq, wt2_seq, resnums1, resnums2):
    """Add monomer index, WT residue, and PDB resnum fields to each row."""
    idx_1tnf = 0
    idx_6ooy = 0
    for row in data:
        if row['aa_1tnf'] != '-':
            row['mono_idx_1tnf'] = idx_1tnf
            row['wt_1tnf'] = wt1_seq[idx_1tnf] if idx_1tnf < len(wt1_seq) else '?'
            row['pdb_resnum_1tnf'] = resnums1[idx_1tnf] if idx_1tnf < len(resnums1) else None
            idx_1tnf += 1
        else:
            row['mono_idx_1tnf'] = None
            row['wt_1tnf'] = '-'
            row['pdb_resnum_1tnf'] = None
        if row['aa_6ooy'] != '-':
            row['mono_idx_6ooy'] = idx_6ooy
            row['wt_6ooy'] = wt2_seq[idx_6ooy] if idx_6ooy < len(wt2_seq) else '?'
            row['pdb_resnum_6ooy'] = resnums2[idx_6ooy] if idx_6ooy < len(resnums2) else None
            idx_6ooy += 1
        else:
            row['mono_idx_6ooy'] = None
            row['wt_6ooy'] = '-'
            row['pdb_resnum_6ooy'] = None
    return idx_1tnf, idx_6ooy


def generate_html(data, n_res_1tnf, n_res_6ooy):
    """Generate the complete HTML string with Mol* viewers."""
    def _row(d):
        return {
            'ap': d['aligned_pos'],
            'a1': d['aa_1tnf'],
            'a2': d['aa_6ooy'],
            't1': d['top_aa_1tnf'],
            'f1': round(d['freq_1tnf'], 3),
            't2': d['top_aa_6ooy'],
            'f2': round(d['freq_6ooy'], 3),
            'jsd': round(d['js_divergence'], 4),
            'mi1': d['mono_idx_1tnf'],
            'mi2': d['mono_idx_6ooy'],
            'w1': d['wt_1tnf'],
            'w2': d['wt_6ooy'],
            'r1': d['pdb_resnum_1tnf'],
            'r2': d['pdb_resnum_6ooy'],
        }

    positions_json = json.dumps([_row(d) for d in data])

    top_divergent = sorted(data, key=lambda x: x['js_divergence'], reverse=True)[:20]
    top_json = json.dumps([_row(d) for d in top_divergent])

    js_values = [d['js_divergence'] for d in data]
    mean_js = sum(js_values) / len(js_values)
    max_js = max(js_values)
    n_high = sum(1 for v in js_values if v > 0.5)
    n_aligned = len(data)

    # Use a raw string approach to avoid f-string brace issues with JS
    # We'll use placeholder tokens and replace them
    html_template = r"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>ProteinMPNN Divergence: 1TNF vs 6OOY (TNF-alpha)</title>
<link rel="stylesheet" type="text/css" href="https://cdn.jsdelivr.net/npm/molstar@4.5.0/build/viewer/molstar.css">
<script type="text/javascript" src="https://cdn.jsdelivr.net/npm/molstar@4.5.0/build/viewer/molstar.js"></script>
<style>
* { margin: 0; padding: 0; box-sizing: border-box; }
body { font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif; background: #f5f5f5; color: #333; }

header {
    background: linear-gradient(135deg, #1a1a2e, #16213e);
    color: white; padding: 24px 32px;
}
header h1 { font-size: 22px; font-weight: 600; margin-bottom: 6px; }
header p { font-size: 14px; opacity: 0.8; }

.stats-bar {
    display: flex; gap: 32px; padding: 12px 32px;
    background: #e8eaf6; font-size: 13px; flex-wrap: wrap;
}
.stats-bar .stat { display: flex; gap: 6px; align-items: center; }
.stats-bar .stat-label { color: #666; }
.stats-bar .stat-value { font-weight: 600; color: #1a1a2e; }

.viewers-section {
    display: flex; gap: 16px; padding: 16px 32px;
    max-width: 1800px; margin: 0 auto;
}
.viewer-container {
    flex: 1; min-width: 0;
    background: white; border-radius: 8px; box-shadow: 0 1px 4px rgba(0,0,0,0.1);
    overflow: hidden;
}
.viewer-header {
    padding: 10px 16px; background: #f8f9fa;
    font-weight: 600; font-size: 14px;
    border-bottom: 1px solid #eee;
    display: flex; justify-content: space-between; align-items: center;
}
.viewer-header .pdb-id { color: #1565c0; }

/* Mol* viewer container - needs position:relative and fixed height */
.viewer-box {
    width: 100%; height: 520px; position: relative;
}

.legend-bar {
    display: flex; align-items: center; justify-content: center;
    gap: 12px; padding: 12px 32px; font-size: 12px; color: #666;
}
.legend-gradient {
    width: 200px; height: 14px; border-radius: 3px;
    border: 1px solid #ccc;
}
.legend-gradient-uncertainty {
    /* Mol* uncertainty theme: blue → cyan → green → yellow → orange */
    background: linear-gradient(to right, #3d3dff, #30d5c8, #7cfc00, #ffd700, #ff8c00);
}

.section { padding: 16px 32px; max-width: 1800px; margin: 0 auto; }
.section h2 {
    font-size: 16px; font-weight: 600; margin-bottom: 12px;
    padding-bottom: 6px; border-bottom: 2px solid #1a1a2e;
    display: inline-block;
}

.alignment-wrap {
    background: white; border-radius: 8px; padding: 16px;
    box-shadow: 0 1px 4px rgba(0,0,0,0.1);
    overflow-x: auto; font-family: 'Courier New', monospace; font-size: 12px;
}
.aln-block { margin-bottom: 12px; white-space: nowrap; }
.aln-label { display: inline-block; width: 60px; text-align: right; margin-right: 8px; color: #888; font-size: 11px; }
.aln-pos-row { margin-bottom: 2px; }
.aln-cell {
    display: inline-block; width: 16px; height: 18px;
    text-align: center; line-height: 18px; font-size: 11px;
    border-radius: 2px; cursor: pointer; font-weight: 500;
}
.aln-cell:hover { outline: 2px solid #333; z-index: 1; position: relative; }
.aln-cell.gap { color: #bbb; background: #f5f5f5 !important; }

.table-container {
    background: white; border-radius: 8px; padding: 16px;
    box-shadow: 0 1px 4px rgba(0,0,0,0.1); overflow-x: auto;
}
table { width: 100%; border-collapse: collapse; font-size: 13px; }
thead th {
    text-align: left; padding: 8px 12px; background: #f8f9fa;
    border-bottom: 2px solid #ddd; font-weight: 600; font-size: 12px;
    white-space: nowrap; cursor: pointer; user-select: none;
}
thead th:hover { background: #e8eaf6; }
tbody td { padding: 8px 12px; border-bottom: 1px solid #f0f0f0; }
tbody tr { cursor: pointer; transition: background 0.15s; }
tbody tr:hover { background: #e3f2fd; }
tbody tr.selected { background: #bbdefb; }
.jsd-bar {
    display: inline-block; height: 12px; border-radius: 2px;
    min-width: 2px; vertical-align: middle;
}
.aa-badge {
    display: inline-block; padding: 1px 6px; border-radius: 3px;
    font-family: 'Courier New', monospace; font-weight: 600; font-size: 12px;
}
.freq-small { font-size: 11px; color: #888; margin-left: 4px; }

footer {
    text-align: center; padding: 20px 32px; font-size: 12px; color: #999;
    border-top: 1px solid #eee; margin-top: 24px;
}

.tooltip {
    position: fixed; background: #333; color: white;
    padding: 8px 12px; border-radius: 6px; font-size: 12px;
    pointer-events: none; z-index: 1000; max-width: 280px;
    box-shadow: 0 2px 8px rgba(0,0,0,0.3); display: none;
    line-height: 1.5;
}
</style>
</head>
<body>

<header>
    <h1>ProteinMPNN Divergence Analysis: 1TNF vs 6OOY</h1>
    <p>TNF-alpha homotrimers &mdash; 10 ProteinMPNN runs per structure, Jensen-Shannon divergence of designed sequence distributions</p>
</header>

<div class="stats-bar">
    <div class="stat"><span class="stat-label">Aligned positions:</span><span class="stat-value">__N_ALIGNED__</span></div>
    <div class="stat"><span class="stat-label">1TNF residues:</span><span class="stat-value">__N_RES_1TNF__</span></div>
    <div class="stat"><span class="stat-label">6OOY residues:</span><span class="stat-value">__N_RES_6OOY__</span></div>
    <div class="stat"><span class="stat-label">Mean JS-div:</span><span class="stat-value">__MEAN_JS__</span></div>
    <div class="stat"><span class="stat-label">Max JS-div:</span><span class="stat-value">__MAX_JS__</span></div>
    <div class="stat"><span class="stat-label">Positions with JS &gt; 0.5:</span><span class="stat-value">__N_HIGH__</span></div>
</div>

<div class="viewers-section">
    <div class="viewer-container">
        <div class="viewer-header">
            <span><span class="pdb-id">1TNF</span> &mdash; Human TNF-alpha</span>
            <span style="font-size:12px;color:#888">__N_RES_1TNF__ res/chain</span>
        </div>
        <div class="viewer-box" id="viewer-1tnf"></div>
    </div>
    <div class="viewer-container">
        <div class="viewer-header">
            <span><span class="pdb-id">6OOY</span> &mdash; Human TNF-alpha</span>
            <span style="font-size:12px;color:#888">__N_RES_6OOY__ res/chain</span>
        </div>
        <div class="viewer-box" id="viewer-6ooy"></div>
    </div>
</div>

<div class="legend-bar">
    <span>Low B-factor (conserved)</span>
    <div class="legend-gradient legend-gradient-uncertainty"></div>
    <span>High B-factor (divergent)</span>
    <span style="margin-left:16px;font-style:italic;color:#999">Mol* "uncertainty" color theme on replaced B-factors</span>
</div>

<div class="section">
    <h2>Aligned Consensus Sequences</h2>
    <div class="alignment-wrap" id="alignment"></div>
</div>

<div class="section">
    <h2>Top 20 Divergent Positions</h2>
    <div class="table-container">
        <table id="div-table">
            <thead>
                <tr>
                    <th>#</th>
                    <th>Aligned Pos</th>
                    <th>WT</th>
                    <th>1TNF Consensus</th>
                    <th>6OOY Consensus</th>
                    <th>Top AA (1TNF)</th>
                    <th>Top AA (6OOY)</th>
                    <th>JS Divergence</th>
                </tr>
            </thead>
            <tbody id="div-tbody"></tbody>
        </table>
    </div>
</div>

<footer>
    Generated from ProteinMPNN analysis &mdash; 1TNF (PDB) vs 6OOY (PDB) TNF-alpha homotrimers
</footer>

<div class="tooltip" id="tooltip"></div>

<script>
// ===== Embedded data =====
const POSITIONS = __POSITIONS_JSON__;
const TOP_DIVERGENT = __TOP_JSON__;

// Correct monomer_index -> PDB_resnum mapping (accounts for ProteinMPNN cleaning)
// These are derived from POSITIONS data which was precomputed in Python.
const CLEANED_RESNUMS_1TNF = POSITIONS.filter(p => p.r1 !== null).map(p => p.r1);
const CLEANED_RESNUMS_6OOY = POSITIONS.filter(p => p.r2 !== null).map(p => p.r2);

// ===== Color mapping (for alignment/table display) =====
function divToRgb(div) {
    // Blue (0) -> White (0.5) -> Red (1)
    let r, g, b;
    if (div <= 0.5) {
        const t = div * 2;
        r = Math.round(t * 255);
        g = Math.round(t * 255);
        b = 255;
    } else {
        const t = (div - 0.5) * 2;
        r = 255;
        g = Math.round((1 - t) * 255);
        b = Math.round((1 - t) * 255);
    }
    return [r, g, b];
}

function divToCss(div) {
    const [r, g, b] = divToRgb(div);
    return `rgb(${r},${g},${b})`;
}

function textColorForDiv(div) {
    const [r, g, b] = divToRgb(div);
    const lum = 0.299 * r + 0.587 * g + 0.114 * b;
    return lum > 140 ? '#222' : '#fff';
}

// ===== PDB parsing =====
function buildResnumDivMap(positions, resnumKey) {
    const m = {};
    for (const p of positions) {
        const resnum = p[resnumKey];
        if (resnum !== null) {
            m[resnum] = p.jsd;
        }
    }
    return m;
}

// Replace B-factor column in PDB with divergence values (scaled 0-99)
function replaceBfactors(pdbText, divMap) {
    return pdbText.split('\n').map(line => {
        if (line.startsWith('ATOM') || line.startsWith('HETATM')) {
            const resNum = parseInt(line.substring(22, 26).trim());
            const div = divMap[resNum];
            if (div !== undefined) {
                const bfac = (div * 99).toFixed(2).padStart(6);
                return line.substring(0, 60) + bfac + line.substring(66);
            } else {
                // Residues not in our mapping get 0
                return line.substring(0, 60) + '  0.00' + line.substring(66);
            }
        }
        return line;
    }).join('\n');
}

// ===== Mol* viewers =====
let plugin1 = null, plugin2 = null;

async function fetchPdb(pdbId) {
    const url = `https://files.rcsb.org/download/${pdbId}.pdb`;
    const resp = await fetch(url);
    if (!resp.ok) throw new Error(`Failed to fetch ${pdbId}: ${resp.status}`);
    return await resp.text();
}

// ===== Structural superposition (Kabsch algorithm) =====
// Extract CA atom coords for chain A, keyed by residue number
function getCAcoords(pdbText) {
    const cas = {};
    for (const line of pdbText.split('\n')) {
        if (!line.startsWith('ATOM')) continue;
        const atomName = line.substring(12, 16).trim();
        const chain = line[21];
        if (atomName !== 'CA' || (chain !== 'A' && chain !== ' ')) continue;
        const resNum = parseInt(line.substring(22, 26).trim());
        const x = parseFloat(line.substring(30, 38));
        const y = parseFloat(line.substring(38, 46));
        const z = parseFloat(line.substring(46, 54));
        if (!cas[resNum]) cas[resNum] = [x, y, z];
    }
    return cas;
}

// Kabsch: find rotation R and translation t so that R * mobile + t ~= target
// Returns { R: 3x3 array, t: [tx,ty,tz] }
function kabschAlign(targetCoords, mobileCoords) {
    const n = targetCoords.length;
    // Compute centroids
    let ct = [0,0,0], cm = [0,0,0];
    for (let i = 0; i < n; i++) {
        for (let j = 0; j < 3; j++) { ct[j] += targetCoords[i][j]; cm[j] += mobileCoords[i][j]; }
    }
    ct = ct.map(v => v/n);
    cm = cm.map(v => v/n);

    // Center the coordinates
    const P = mobileCoords.map(c => [c[0]-cm[0], c[1]-cm[1], c[2]-cm[2]]);
    const Q = targetCoords.map(c => [c[0]-ct[0], c[1]-ct[1], c[2]-ct[2]]);

    // Compute cross-covariance matrix H = P^T * Q (3x3)
    const H = [[0,0,0],[0,0,0],[0,0,0]];
    for (let i = 0; i < n; i++) {
        for (let r = 0; r < 3; r++) {
            for (let c = 0; c < 3; c++) {
                H[r][c] += P[i][r] * Q[i][c];
            }
        }
    }

    // SVD of H using Jacobi iterations (simple 3x3 SVD)
    const { U, V } = svd3x3(H);

    // R = V * U^T, ensuring proper rotation (det=+1)
    let R = matMul3(V, transpose3(U));
    if (det3(R) < 0) {
        // Flip last column of V
        V[0][2] *= -1; V[1][2] *= -1; V[2][2] *= -1;
        R = matMul3(V, transpose3(U));
    }

    // Translation: t = ct - R * cm
    const Rcm = matVec3(R, cm);
    const t = [ct[0]-Rcm[0], ct[1]-Rcm[1], ct[2]-Rcm[2]];

    return { R, t };
}

// 3x3 matrix helpers
function transpose3(M) { return [[M[0][0],M[1][0],M[2][0]],[M[0][1],M[1][1],M[2][1]],[M[0][2],M[1][2],M[2][2]]]; }
function matMul3(A,B) {
    const C = [[0,0,0],[0,0,0],[0,0,0]];
    for (let i=0;i<3;i++) for (let j=0;j<3;j++) for (let k=0;k<3;k++) C[i][j]+=A[i][k]*B[k][j];
    return C;
}
function matVec3(M,v) { return [M[0][0]*v[0]+M[0][1]*v[1]+M[0][2]*v[2], M[1][0]*v[0]+M[1][1]*v[1]+M[1][2]*v[2], M[2][0]*v[0]+M[2][1]*v[1]+M[2][2]*v[2]]; }
function det3(M) { return M[0][0]*(M[1][1]*M[2][2]-M[1][2]*M[2][1])-M[0][1]*(M[1][0]*M[2][2]-M[1][2]*M[2][0])+M[0][2]*(M[1][0]*M[2][1]-M[1][1]*M[2][0]); }

// Simple 3x3 SVD via eigendecomposition of H^T*H
function svd3x3(H) {
    const Ht = transpose3(H);
    const HtH = matMul3(Ht, H);
    const { eigenvalues, eigenvectors } = eigenDecomp3(HtH);

    // V = eigenvectors (columns), sorted by descending eigenvalue
    const idx = [0,1,2].sort((a,b) => eigenvalues[b]-eigenvalues[a]);
    const V = [[0,0,0],[0,0,0],[0,0,0]];
    for (let i=0;i<3;i++) for (let j=0;j<3;j++) V[i][j] = eigenvectors[i][idx[j]];

    // U = H * V * Sigma^-1
    const HV = matMul3(H, V);
    const U = [[0,0,0],[0,0,0],[0,0,0]];
    for (let j=0;j<3;j++) {
        const sigma = Math.sqrt(Math.max(eigenvalues[idx[j]], 1e-10));
        for (let i=0;i<3;i++) U[i][j] = HV[i][j] / sigma;
    }

    return { U, V };
}

// 3x3 symmetric matrix eigendecomposition via Jacobi iteration
function eigenDecomp3(A) {
    const a = [[A[0][0],A[0][1],A[0][2]],[A[1][0],A[1][1],A[1][2]],[A[2][0],A[2][1],A[2][2]]];
    const v = [[1,0,0],[0,1,0],[0,0,1]];
    for (let iter = 0; iter < 50; iter++) {
        // Find largest off-diagonal
        let p=0, q=1;
        if (Math.abs(a[0][2]) > Math.abs(a[p][q])) { p=0; q=2; }
        if (Math.abs(a[1][2]) > Math.abs(a[p][q])) { p=1; q=2; }
        if (Math.abs(a[p][q]) < 1e-12) break;

        const theta = 0.5 * Math.atan2(2*a[p][q], a[q][q]-a[p][p]);
        const c = Math.cos(theta), s = Math.sin(theta);

        // Givens rotation
        const G = [[1,0,0],[0,1,0],[0,0,1]];
        G[p][p]=c; G[q][q]=c; G[p][q]=s; G[q][p]=-s;

        // a = G^T * a * G
        const tmp = matMul3(transpose3(G), matMul3(a, G));
        for (let i=0;i<3;i++) for (let j=0;j<3;j++) a[i][j]=tmp[i][j];

        // v = v * G
        const tmpv = matMul3(v, G);
        for (let i=0;i<3;i++) for (let j=0;j<3;j++) v[i][j]=tmpv[i][j];
    }
    return { eigenvalues: [a[0][0], a[1][1], a[2][2]], eigenvectors: v };
}

// Apply rotation R and translation t to all ATOM/HETATM in PDB text
function applyTransformPdb(pdbText, R, t) {
    return pdbText.split('\n').map(line => {
        if (!line.startsWith('ATOM') && !line.startsWith('HETATM')) return line;
        const x = parseFloat(line.substring(30, 38));
        const y = parseFloat(line.substring(38, 46));
        const z = parseFloat(line.substring(46, 54));
        const nx = R[0][0]*x + R[0][1]*y + R[0][2]*z + t[0];
        const ny = R[1][0]*x + R[1][1]*y + R[1][2]*z + t[1];
        const nz = R[2][0]*x + R[2][1]*y + R[2][2]*z + t[2];
        return line.substring(0,30) + nx.toFixed(3).padStart(8) + ny.toFixed(3).padStart(8) + nz.toFixed(3).padStart(8) + line.substring(54);
    }).join('\n');
}

// Superpose mobile PDB onto target PDB using aligned residues
// Uses precomputed PDB residue numbers (r1, r2) from POSITIONS data.
function superposePdb(targetPdb, mobilePdb) {
    const targetCAs = getCAcoords(targetPdb);
    const mobileCAs = getCAcoords(mobilePdb);

    // Collect paired CA coords for aligned positions using correct PDB resnums
    const tCoords = [], mCoords = [];
    for (const p of POSITIONS) {
        if (p.r1 === null || p.r2 === null) continue;
        if (targetCAs[p.r1] && mobileCAs[p.r2]) {
            tCoords.push(targetCAs[p.r1]);
            mCoords.push(mobileCAs[p.r2]);
        }
    }

    if (tCoords.length < 3) return mobilePdb;
    const { R, t } = kabschAlign(tCoords, mCoords);
    console.log(`Superposition: ${tCoords.length} paired CAs, applied rotation+translation`);
    return applyTransformPdb(mobilePdb, R, t);
}

async function applyUncertaintyColoring(plugin) {
    // Wait briefly for the representation state tree to settle
    await new Promise(r => setTimeout(r, 500));
    const structs = plugin.managers.structure.hierarchy.current.structures;
    for (const s of structs) {
        for (const c of s.components) {
            for (const repr of c.representations) {
                const cell = repr.cell;
                const ref = cell.transform.ref;
                const oldParams = cell.transform.params;
                const newParams = {
                    ...oldParams,
                    colorTheme: { name: 'uncertainty', params: {} }
                };
                await plugin.state.data.build().to(ref).update(newParams).commit();
            }
        }
    }
}

async function loadPdbIntoViewer(elementId, pdbId, pdbText, resnumKey) {
    const viewer = await molstar.Viewer.create(elementId, {
        layoutIsExpanded: false,
        layoutShowControls: true,
        layoutShowSequence: true,
        layoutShowRemoteState: false,
        layoutShowLog: false,
        viewportShowExpand: true,
        viewportShowControls: true,
        viewportShowSettings: true,
        viewportShowSelectionMode: true,
        collapseLeftPanel: true,
        collapseRightPanel: true,
    });
    const plugin = viewer.plugin;

    try {
        const divMap = buildResnumDivMap(POSITIONS, resnumKey);
        const modifiedPdb = replaceBfactors(pdbText, divMap);

        const rawData = await plugin.builders.data.rawData({ data: modifiedPdb, label: pdbId });
        const trajectory = await plugin.builders.structure.parseTrajectory(rawData, 'pdb');
        await plugin.builders.structure.hierarchy.applyPreset(trajectory, 'default', { showUnitcell: false });
        await applyUncertaintyColoring(plugin);

        console.log(`${pdbId}: loaded with uncertainty coloring`);
        return { plugin };
    } catch (e) {
        console.error(`Failed to load ${pdbId}:`, e);
        return { plugin };
    }
}

async function initViewers() {
    // Fetch both PDBs
    const [pdb1tnf, pdb6ooy] = await Promise.all([fetchPdb('1TNF'), fetchPdb('6OOY')]);

    // Superpose 6OOY onto 1TNF coordinate frame using aligned CA atoms
    const superposed6ooy = superposePdb(pdb1tnf, pdb6ooy);

    // Load both into viewers (can be parallel since they're independent)
    // Uses 'r1'/'r2' keys which hold correct PDB residue numbers (precomputed)
    const [r1, r2] = await Promise.all([
        loadPdbIntoViewer('viewer-1tnf', '1TNF', pdb1tnf, 'r1'),
        loadPdbIntoViewer('viewer-6ooy', '6OOY', superposed6ooy, 'r2'),
    ]);
    plugin1 = r1.plugin;
    plugin2 = r2.plugin;
}

// ===== Focus on residue via Mol* =====
// Build a StructureElement.Loci for a single chain's residue (first matching unit).
// Using only one chain avoids the bounding sphere spanning the entire homotrimer.
function buildResidueLoci(plugin, resNum) {
    const structs = plugin.managers.structure.hierarchy.current.structures;
    if (!structs.length) return null;
    const struct = structs[0].cell.obj?.data;
    if (!struct) return null;

    for (const unit of struct.units) {
        if (!unit.conformation) continue;
        const elements = unit.elements;
        const residueIndex = unit.residueIndex;
        const auth_seq_id = unit.model.atomicHierarchy.residues.auth_seq_id;

        const indices = [];
        for (let eIdx = 0; eIdx < elements.length; eIdx++) {
            const eI = elements[eIdx];
            const rI = residueIndex[eI];
            if (auth_seq_id.value(rI) === resNum) {
                indices.push(eIdx);
            }
        }
        if (indices.length > 0) {
            return {
                kind: 'element-loci',
                structure: struct,
                elements: [{ unit, indices: Int32Array.from(indices) }]
            };
        }
    }
    return null;
}

function focusMolstarOnResidue(plugin, resNum) {
    try {
        const loci = buildResidueLoci(plugin, resNum);
        if (!loci) return;
        plugin.managers.structure.focus.setFromLoci(loci);
        setTimeout(() => { plugin.managers.camera.focusLoci(loci); }, 50);
    } catch (e) {
        console.warn('Focus failed:', resNum, e.message || e);
    }
}

function focusResidue(pos) {
    const p = POSITIONS.find(d => d.ap === pos);
    if (!p) return;

    if (p.r1 !== null && plugin1) {
        focusMolstarOnResidue(plugin1, p.r1);
    }
    if (p.r2 !== null && plugin2) {
        focusMolstarOnResidue(plugin2, p.r2);
    }
}

// ===== Alignment display =====
function renderAlignment() {
    const container = document.getElementById('alignment');
    const lineLen = 50;
    let html = '';

    for (let start = 0; start < POSITIONS.length; start += lineLen) {
        const chunk = POSITIONS.slice(start, start + lineLen);

        html += '<div class="aln-block">';
        html += '<div class="aln-pos-row"><span class="aln-label"></span>';
        for (let i = 0; i < chunk.length; i++) {
            const pos = start + i + 1;
            if (pos === 1 || pos % 10 === 0) {
                html += `<span class="aln-cell" style="color:#999;font-size:10px;background:none">${pos}</span>`;
            } else {
                html += '<span class="aln-cell" style="background:none"></span>';
            }
        }
        html += '</div>';

        // WT row (shared between both structures, show any differences)
        html += '<div><span class="aln-label" style="color:#555;font-weight:600">WT</span>';
        for (const p of chunk) {
            const wt = p.w1 !== '-' ? p.w1 : p.w2;
            const cls = wt === '-' ? ' gap' : '';
            const wtMismatch = p.w1 !== '-' && p.w2 !== '-' && p.w1 !== p.w2;
            const bg = wt === '-' ? '' : wtMismatch ? 'background:#fff3e0;color:#e65100;font-weight:700' : 'background:#f5f5f5;color:#555';
            const title = wtMismatch ? `WT mismatch: 1TNF=${p.w1} 6OOY=${p.w2}` : `WT: ${wt}`;
            html += `<span class="aln-cell${cls}" style="${bg}" data-pos="${p.ap}" title="${title}">${wtMismatch ? p.w1+'/'+p.w2 : wt}</span>`;
        }
        html += '</div>';

        html += '<div><span class="aln-label" style="color:#1565c0">1TNF design</span>';
        for (const p of chunk) {
            const cls = p.a1 === '-' ? ' gap' : '';
            const bg = p.a1 === '-' ? '' : `background:${divToCss(p.jsd)};color:${textColorForDiv(p.jsd)}`;
            html += `<span class="aln-cell${cls}" style="${bg}" data-pos="${p.ap}" title="Pos ${p.ap}: cons=${p.a1} WT=${p.w1} JS=${p.jsd}">${p.a1}</span>`;
        }
        html += '</div>';

        html += '<div><span class="aln-label" style="color:#c62828">6OOY design</span>';
        for (const p of chunk) {
            const cls = p.a2 === '-' ? ' gap' : '';
            const bg = p.a2 === '-' ? '' : `background:${divToCss(p.jsd)};color:${textColorForDiv(p.jsd)}`;
            html += `<span class="aln-cell${cls}" style="${bg}" data-pos="${p.ap}" title="Pos ${p.ap}: cons=${p.a2} WT=${p.w2} JS=${p.jsd}">${p.a2}</span>`;
        }
        html += '</div>';

        html += '<div><span class="aln-label"></span>';
        for (const p of chunk) {
            let sym = ' ';
            if (p.a1 === p.a2 && p.a1 !== '-') sym = '|';
            else if (p.a1 !== '-' && p.a2 !== '-') sym = '.';
            html += `<span class="aln-cell" style="background:none;color:#aaa;font-size:10px">${sym}</span>`;
        }
        html += '</div>';
        html += '</div>';
    }

    container.innerHTML = html;
    container.addEventListener('click', function(e) {
        const cell = e.target.closest('.aln-cell[data-pos]');
        if (cell) focusResidue(parseInt(cell.dataset.pos));
    });
}

// ===== Table =====
function renderTable() {
    const tbody = document.getElementById('div-tbody');
    let html = '';

    TOP_DIVERGENT.forEach((d, i) => {
        const resNum = d.r1 !== null ? d.r1 : (d.r2 !== null ? d.r2 : '?');
        const wtAA = d.w1 !== '-' ? d.w1 : d.w2;
        const wtMismatch = d.w1 !== '-' && d.w2 !== '-' && d.w1 !== d.w2;
        const wtLabel = wtMismatch ? `${d.w1}/${d.w2}` : wtAA;
        const wtStyle = wtMismatch ? 'background:#fff3e0;color:#e65100' : 'background:#f5f5f5';
        const res1 = d.r1 !== null ? `${d.a1}${d.r1}` : (d.a1 === '-' ? 'gap' : d.a1 + '?');
        const res2 = d.r2 !== null ? `${d.a2}${d.r2}` : (d.a2 === '-' ? 'gap' : d.a2 + '?');
        const barW = Math.round(d.jsd * 100);
        const barColor = divToCss(d.jsd);

        html += `<tr data-pos="${d.ap}">
            <td>${i + 1}</td>
            <td>${d.ap}</td>
            <td><span class="aa-badge" style="${wtStyle}">${wtLabel}${resNum}</span></td>
            <td><span class="aa-badge" style="background:#e3f2fd">${res1}</span></td>
            <td><span class="aa-badge" style="background:#fce4ec">${res2}</span></td>
            <td><span class="aa-badge" style="background:#e8f5e9">${d.t1}</span><span class="freq-small">${(d.f1 * 100).toFixed(0)}%</span></td>
            <td><span class="aa-badge" style="background:#e8f5e9">${d.t2}</span><span class="freq-small">${(d.f2 * 100).toFixed(0)}%</span></td>
            <td>
                <span class="jsd-bar" style="width:${barW}px;background:${barColor}"></span>
                ${d.jsd.toFixed(3)}
            </td>
        </tr>`;
    });
    tbody.innerHTML = html;

    tbody.addEventListener('click', function(e) {
        const row = e.target.closest('tr[data-pos]');
        if (row) {
            tbody.querySelectorAll('tr').forEach(r => r.classList.remove('selected'));
            row.classList.add('selected');
            focusResidue(parseInt(row.dataset.pos));
        }
    });
}

// ===== Tooltip =====
const tooltip = document.getElementById('tooltip');
document.addEventListener('mousemove', function(e) {
    const cell = e.target.closest('.aln-cell[data-pos]');
    if (cell) {
        const pos = parseInt(cell.dataset.pos);
        const p = POSITIONS.find(d => d.ap === pos);
        if (p) {
            const wt = p.w1 !== '-' ? p.w1 : p.w2;
            const wtNote = (p.w1 !== '-' && p.w2 !== '-' && p.w1 !== p.w2) ? ` (1TNF=${p.w1}, 6OOY=${p.w2})` : '';
            tooltip.innerHTML = `<b>Position ${pos}</b> &mdash; WT: ${wt}${wtNote}<br>` +
                `1TNF design: ${p.a1} (top: ${p.t1} ${(p.f1*100).toFixed(0)}%)<br>` +
                `6OOY design: ${p.a2} (top: ${p.t2} ${(p.f2*100).toFixed(0)}%)<br>` +
                `JS divergence: <b>${p.jsd.toFixed(4)}</b>`;
            tooltip.style.display = 'block';
            tooltip.style.left = (e.clientX + 12) + 'px';
            tooltip.style.top = (e.clientY + 12) + 'px';
        }
    } else {
        tooltip.style.display = 'none';
    }
});

// ===== Residue number mapping between structures =====
// Both structures are TNF-alpha with the same PDB residue numbering scheme.
// Map directly by PDB residue number (same number = same residue).
let resmap1to2 = {};  // 1TNF PDB resnum -> 6OOY PDB resnum
let resmap2to1 = {};  // 6OOY PDB resnum -> 1TNF PDB resnum

function buildResnumMaps() {
    const set1 = new Set(CLEANED_RESNUMS_1TNF);
    const set2 = new Set(CLEANED_RESNUMS_6OOY);
    for (const r of set1) {
        if (set2.has(r)) {
            resmap1to2[r] = r;
            resmap2to1[r] = r;
        }
    }
}

// ===== Camera sync =====
// Structures are superposed (same coordinate frame), so we can copy
// the full camera state (position + target + up) between viewers.
let cameraSyncPaused = false;

function setupCameraSync() {
    let syncing = false;

    function syncCamera(srcPlugin, tgtPlugin) {
        if (syncing || cameraSyncPaused) return;
        syncing = true;
        try {
            const srcSnap = srcPlugin.canvas3d.camera.getSnapshot();

            tgtPlugin.canvas3d.camera.setState({
                position: new Float64Array(srcSnap.position),
                target: new Float64Array(srcSnap.target),
                up: new Float64Array(srcSnap.up),
            });
            tgtPlugin.canvas3d.requestDraw();
        } finally {
            syncing = false;
        }
    }

    // Track which viewer the pointer is over
    let activePlugin = null;
    document.getElementById('viewer-1tnf').addEventListener('pointerenter', () => { activePlugin = plugin1; });
    document.getElementById('viewer-6ooy').addEventListener('pointerenter', () => { activePlugin = plugin2; });
    document.getElementById('viewer-1tnf').addEventListener('pointerleave', () => { if (activePlugin === plugin1) activePlugin = null; });
    document.getElementById('viewer-6ooy').addEventListener('pointerleave', () => { if (activePlugin === plugin2) activePlugin = null; });

    plugin1.canvas3d.didDraw.subscribe(() => {
        if (activePlugin === plugin1) syncCamera(plugin1, plugin2);
    });
    plugin2.canvas3d.didDraw.subscribe(() => {
        if (activePlugin === plugin2) syncCamera(plugin2, plugin1);
    });
}

// ===== Focus/selection sync =====
function setupFocusSync() {
    let syncingFocus = false;

    // Parse residue number from Mol* focus label like "GLN 67 | B"
    function getResNumFromLabel(label) {
        if (!label) return null;
        const m = label.match(/\b(\d+)\b/);
        return m ? parseInt(m[1]) : null;
    }

    // Trigger focus+camera on a single chain's residue in the given viewer
    function focusAndSelect(plugin, resNum) {
        const loci = buildResidueLoci(plugin, resNum);
        if (!loci) return;

        // Pause camera sync so focusLoci on this viewer isn't overridden
        cameraSyncPaused = true;
        // setFromLoci triggers the sidechain/nearby-residue display
        plugin.managers.structure.focus.setFromLoci(loci);
        // Small delay to let focus representation settle, then zoom camera
        setTimeout(() => {
            plugin.managers.camera.focusLoci(loci);
            // Re-enable camera sync after zoom settles
            setTimeout(() => { cameraSyncPaused = false; }, 300);
        }, 50);
    }

    // When focus changes in viewer 1, mirror to viewer 2
    plugin1.managers.structure.focus.behaviors.current.subscribe(val => {
        if (syncingFocus || !val) return;
        const resNum = getResNumFromLabel(val.label);
        if (resNum === null) return;
        const mappedRes = resmap1to2[resNum];
        if (mappedRes === undefined) return;
        syncingFocus = true;
        try { focusAndSelect(plugin2, mappedRes); } finally { syncingFocus = false; }
    });

    // When focus changes in viewer 2, mirror to viewer 1
    plugin2.managers.structure.focus.behaviors.current.subscribe(val => {
        if (syncingFocus || !val) return;
        const resNum = getResNumFromLabel(val.label);
        if (resNum === null) return;
        const mappedRes = resmap2to1[resNum];
        if (mappedRes === undefined) return;
        syncingFocus = true;
        try { focusAndSelect(plugin1, mappedRes); } finally { syncingFocus = false; }
    });
}

// ===== Init =====
async function init() {
    renderAlignment();
    await initViewers();
    buildResnumMaps();
    setupCameraSync();
    setupFocusSync();
    renderTable();
}

init();
</script>
</body>
</html>"""

    # Replace placeholder tokens
    html = html_template
    html = html.replace('__N_ALIGNED__', str(n_aligned))
    html = html.replace('__N_RES_1TNF__', str(n_res_1tnf))
    html = html.replace('__N_RES_6OOY__', str(n_res_6ooy))
    html = html.replace('__MEAN_JS__', f'{mean_js:.3f}')
    html = html.replace('__MAX_JS__', f'{max_js:.3f}')
    html = html.replace('__N_HIGH__', str(n_high))
    html = html.replace('__POSITIONS_JSON__', positions_json)
    html = html.replace('__TOP_JSON__', top_json)

    return html


def main():
    filepath = '/Users/rafal/repos/proteinmpnn/divergence_analysis.csv'
    data = read_data(filepath)

    wt1 = get_wt_monomer('/Users/rafal/repos/proteinmpnn/1tnf_results.csv')
    wt2 = get_wt_monomer('/Users/rafal/repos/proteinmpnn/6ooy_homomer_results.csv')

    # Get correct PDB residue numbers by aligning cleaned sequences to PDB
    print('Fetching PDB structures for residue number mapping...')
    pdb_resnums_1tnf, pdb_seq_1tnf = get_pdb_chain_a('1TNF')
    pdb_resnums_6ooy, pdb_seq_6ooy = get_pdb_chain_a('6OOY')

    cleaned_resnums_1tnf = map_cleaned_to_pdb_resnums(wt1, pdb_resnums_1tnf, pdb_seq_1tnf)
    cleaned_resnums_6ooy = map_cleaned_to_pdb_resnums(wt2, pdb_resnums_6ooy, pdb_seq_6ooy)

    print(f'1TNF: {len(cleaned_resnums_1tnf)} cleaned residues (PDB {cleaned_resnums_1tnf[0]}-{cleaned_resnums_1tnf[-1]})')
    print(f'6OOY: {len(cleaned_resnums_6ooy)} cleaned residues (PDB {cleaned_resnums_6ooy[0]}-{cleaned_resnums_6ooy[-1]})')

    n1, n2 = compute_mappings(data, wt1, wt2, cleaned_resnums_1tnf, cleaned_resnums_6ooy)

    html = generate_html(data, n1, n2)

    outpath = '/Users/rafal/repos/proteinmpnn/proteinmpnn_comparison.html'
    with open(outpath, 'w') as f:
        f.write(html)

    print(f'Generated {outpath}')
    print(f'1TNF: {n1} residues, 6OOY: {n2} residues')
    print(f'Aligned positions: {len(data)}')


if __name__ == '__main__':
    main()
