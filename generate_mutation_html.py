#!/usr/bin/env python3
"""Generate HTML visualization of ProteinMPNN mutation sensitivity between 1TNF and 6OOY."""

import csv
import json
import urllib.request

AA3TO1 = {'ALA':'A','CYS':'C','ASP':'D','GLU':'E','PHE':'F','GLY':'G','HIS':'H',
          'ILE':'I','LYS':'K','LEU':'L','MET':'M','ASN':'N','PRO':'P','GLN':'Q',
          'ARG':'R','SER':'S','THR':'T','VAL':'V','TRP':'W','TYR':'Y'}


def main():
    print('Loading data...')
    viz_data = json.load(open('/Users/rafal/repos/proteinmpnn/mutation_viz_data.json'))

    # Read full mutation results — one row per (pdb_resnum, chain) showing worst mutation
    table_data = []
    best_by_pos_chain = {}
    with open('/Users/rafal/repos/proteinmpnn/mutation_ddg_scores.csv') as f:
        reader = csv.DictReader(f)
        for row in reader:
            rn = int(row['pdb_resnum'])
            ch = row['chain']
            abs_dd = float(row['abs_dd'])
            key = (rn, ch)
            if key not in best_by_pos_chain or abs_dd > best_by_pos_chain[key]['abs_dd']:
                best_by_pos_chain[key] = {
                    'rn': rn,
                    'ch': ch,
                    'wt1': row['wt_1tnf'],
                    'wt2': row['wt_6ooy'],
                    'mut': row['mut_aa'],
                    'd1': round(float(row['delta_1tnf']), 2),
                    'd2': round(float(row['delta_6ooy']), 2),
                    'dd': round(float(row['dd_score']), 2),
                    'abs_dd': round(abs_dd, 2),
                }
    table_rows = sorted(best_by_pos_chain.values(), key=lambda r: r['abs_dd'], reverse=True)

    positions_json = json.dumps(viz_data)
    table_json = json.dumps(table_rows)

    # Stats
    max_dd = max(v.get('max_abs_dd', 0) for v in viz_data)
    n_high = sum(1 for v in viz_data if v.get('max_abs_dd', 0) > 3)
    n_positions = len(viz_data)

    html = HTML_TEMPLATE.replace('__POSITIONS__', positions_json)
    html = html.replace('__TABLE_DATA__', table_json)
    html = html.replace('__N_POS__', str(n_positions))
    html = html.replace('__MAX_DD__', f'{max_dd:.1f}')
    html = html.replace('__N_HIGH__', str(n_high))
    html = html.replace('__N_TABLE__', str(len(table_rows)))

    outpath = '/Users/rafal/repos/proteinmpnn/mutation_sensitivity.html'
    with open(outpath, 'w') as f:
        f.write(html)
    print(f'Written to {outpath}')
    print(f'  {n_positions} positions, {len(table_rows)} table rows ({n_positions}x3 chains)')


HTML_TEMPLATE = r"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>Mutation Sensitivity: 1TNF vs 6OOY</title>
<link rel="stylesheet" type="text/css" href="https://cdn.jsdelivr.net/npm/molstar@4.5.0/build/viewer/molstar.css">
<script type="text/javascript" src="https://cdn.jsdelivr.net/npm/molstar@4.5.0/build/viewer/molstar.js"></script>
<style>
* { margin: 0; padding: 0; box-sizing: border-box; }
body { font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif; background: #f5f5f5; color: #333; }

header {
    background: linear-gradient(135deg, #1a1a2e, #16213e);
    color: white; padding: 20px 32px;
}
header h1 { font-size: 1.5em; margin-bottom: 4px; }
header p { opacity: 0.8; font-size: 0.9em; }

.stats-bar {
    display: flex; gap: 32px; padding: 12px 32px;
    background: #e8eaf6; border-bottom: 1px solid #c5cae9;
    font-size: 0.9em;
}
.stat { display: flex; gap: 6px; align-items: center; }
.stat-val { font-weight: 700; color: #1a237e; }

.viewer-row {
    display: grid; grid-template-columns: 1fr 1fr; gap: 8px;
    padding: 12px; background: #fff;
}
.viewer-panel { position: relative; }
.viewer-label {
    position: absolute; top: 8px; left: 8px; z-index: 10;
    background: rgba(0,0,0,0.7); color: white; padding: 4px 12px;
    border-radius: 4px; font-size: 0.85em; font-weight: 600;
}
.viewer-container { width: 100%; height: 420px; position: relative; }

.color-legend {
    display: flex; align-items: center; gap: 8px;
    padding: 8px 32px; background: #fff; font-size: 0.85em;
    border-bottom: 1px solid #eee;
}
.legend-gradient {
    width: 200px; height: 14px;
    background: linear-gradient(to right, #3538A0, #2EA0C8, #41BB66, #D8BF2E, #F06633);
    border-radius: 2px; border: 1px solid #ccc;
}

.table-section { padding: 16px 32px 32px; }
.table-section h2 { margin-bottom: 12px; font-size: 1.1em; }

.controls { display: flex; gap: 12px; margin-bottom: 12px; align-items: center; }
.controls input {
    padding: 6px 12px; border: 1px solid #ccc; border-radius: 4px;
    font-size: 0.9em; width: 200px;
}
.controls label { font-size: 0.85em; cursor: pointer; }

table { width: 100%; border-collapse: collapse; font-size: 0.85em; }
thead { background: #263238; color: white; position: sticky; top: 0; z-index: 5; }
th { padding: 8px 10px; text-align: left; cursor: pointer; user-select: none; white-space: nowrap; }
th:hover { background: #37474f; }
td { padding: 6px 10px; border-bottom: 1px solid #e0e0e0; }
tr:hover { background: #e3f2fd; cursor: pointer; }
tr.selected { background: #bbdefb !important; }
.dd-pos { color: #c62828; font-weight: 600; }
.dd-neg { color: #1565c0; font-weight: 600; }
.dd-bar {
    display: inline-block; height: 12px; vertical-align: middle;
    border-radius: 1px; min-width: 1px;
}
.chain-tag {
    display: inline-block; padding: 1px 5px; border-radius: 3px;
    font-weight: 600; font-size: 0.8em;
}
.chain-A { background: #e3f2fd; color: #1565c0; }
.chain-B { background: #e8f5e9; color: #2e7d32; }
.chain-C { background: #fff3e0; color: #e65100; }

.table-wrap { max-height: 500px; overflow-y: auto; border: 1px solid #e0e0e0; border-radius: 4px; }
</style>
</head>
<body>

<header>
    <h1>Mutation Sensitivity: 1TNF vs 6OOY</h1>
    <p>ProteinMPNN conditional probability analysis &mdash; DDscore = &Delta;score<sub>1TNF</sub> &minus; &Delta;score<sub>6OOY</sub> (per chain, per position)</p>
</header>

<div class="stats-bar">
    <div class="stat"><span>Shared positions:</span><span class="stat-val">__N_POS__</span></div>
    <div class="stat"><span>Independent positions (3 chains):</span><span class="stat-val">__N_TABLE__</span></div>
    <div class="stat"><span>Max |DD|:</span><span class="stat-val">__MAX_DD__</span></div>
    <div class="stat"><span>Positions with |DD| &gt; 3:</span><span class="stat-val">__N_HIGH__</span></div>
</div>

<div class="viewer-row">
    <div class="viewer-panel">
        <div class="viewer-label">1TNF</div>
        <div id="viewer1" class="viewer-container"></div>
    </div>
    <div class="viewer-panel">
        <div class="viewer-label">6OOY</div>
        <div id="viewer2" class="viewer-container"></div>
    </div>
</div>

<div class="color-legend">
    <span>Backbone sensitivity:</span>
    <span style="color:#3538A0;font-weight:600">Low</span>
    <div class="legend-gradient"></div>
    <span style="color:#F06633;font-weight:600">High |DD|</span>
</div>

<div class="table-section">
    <h2>Backbone-Sensitive Mutations (worst per position &times; chain)</h2>
    <div class="controls">
        <input type="text" id="searchBox" placeholder="Filter by residue or AA...">
        <label><input type="checkbox" id="sameWtOnly" checked> Same WT only</label>
        <span id="rowCount" style="font-size:0.85em;color:#666;"></span>
    </div>
    <div class="table-wrap">
        <table>
            <thead>
                <tr>
                    <th data-col="rn">Res</th>
                    <th data-col="ch">Chain</th>
                    <th data-col="wt1">WT<sub>1TNF</sub></th>
                    <th data-col="wt2">WT<sub>6OOY</sub></th>
                    <th data-col="mut">Mut</th>
                    <th data-col="d1">&Delta;<sub>1TNF</sub></th>
                    <th data-col="d2">&Delta;<sub>6OOY</sub></th>
                    <th data-col="dd">DDscore</th>
                    <th data-col="abs_dd">|DD|</th>
                    <th>Bar</th>
                </tr>
            </thead>
            <tbody id="tableBody"></tbody>
        </table>
    </div>
</div>

<script>
const POSITIONS = __POSITIONS__;
const TABLE_DATA = __TABLE_DATA__;

// ========== Mol* viewers ==========
let plugin1 = null, plugin2 = null;

async function initViewers() {
    const opts = {
        layoutIsExpanded: false, layoutShowControls: false,
        layoutShowRemoteState: false, layoutShowSequence: false,
        layoutShowLog: false, viewportShowExpand: false,
        viewportShowSelectionMode: false,
    };
    const v1 = await molstar.Viewer.create('viewer1', opts);
    const v2 = await molstar.Viewer.create('viewer2', opts);
    plugin1 = v1.plugin;
    plugin2 = v2.plugin;

    // Fetch PDBs
    const [pdb1, pdb2] = await Promise.all([
        fetch('https://files.rcsb.org/download/1TNF.pdb').then(r => r.text()),
        fetch('https://files.rcsb.org/download/6OOY.pdb').then(r => r.text()),
    ]);

    // Build per-chain resnum -> DD maps
    const ddMap1 = buildDDMap(POSITIONS, true);
    const ddMap2 = buildDDMap(POSITIONS, false);

    // Replace B-factors with |DD| scores (per chain)
    const pdb1mod = replaceBfactors(pdb1, ddMap1);

    // Superpose 6OOY onto 1TNF
    const pdb2sup = superposePdb(pdb1, pdb2);
    const pdb2mod = replaceBfactors(pdb2sup, ddMap2);

    await loadPdb(v1, pdb1mod, '1TNF');
    await loadPdb(v2, pdb2mod, '6OOY');

    setupCameraSync();
    setupFocusSync();
}

function buildDDMap(positions, is1tnf) {
    // Returns {chain: {resnum: scaledValue}} for B-factor replacement
    // Square-root scaling for better visual contrast
    const maxDD = Math.max(...positions.map(p => p.max_abs_dd || 0));
    const map = {};
    for (const chain of ['A','B','C']) {
        map[chain] = {};
        for (const p of positions) {
            const rn = p.rn;
            if (rn == null) continue;
            const dd = p['abs_dd_' + chain] || 0;
            const frac = dd / Math.max(maxDD, 0.01);
            map[chain][rn] = Math.min(99, Math.sqrt(frac) * 99);
        }
    }
    return map;
}

function replaceBfactors(pdbText, ddMap) {
    return pdbText.split('\n').map(line => {
        if (!line.startsWith('ATOM') && !line.startsWith('HETATM')) return line;
        if (line.length < 66) return line;
        const chain = line[21];
        const resnum = parseInt(line.substring(22, 26));
        const val = (ddMap[chain] && ddMap[chain][resnum]) || 0;
        return line.substring(0, 60) + val.toFixed(2).padStart(6) + line.substring(66);
    }).join('\n');
}

// ========== Kabsch superposition ==========
function getCAcoords(pdb) {
    const cas = {};
    for (const line of pdb.split('\n')) {
        if (!line.startsWith('ATOM')) continue;
        if (line.substring(12, 16).trim() !== 'CA') continue;
        if (line[21] !== 'A') continue;
        const rn = parseInt(line.substring(22, 26));
        cas[rn] = [parseFloat(line.substring(30, 38)), parseFloat(line.substring(38, 46)), parseFloat(line.substring(46, 54))];
    }
    return cas;
}

function superposePdb(targetPdb, mobilePdb) {
    const tCAs = getCAcoords(targetPdb);
    const mCAs = getCAcoords(mobilePdb);
    const tCoords = [], mCoords = [];
    // Use shared PDB resnums from POSITIONS data
    for (const p of POSITIONS) {
        const rn = p.rn;
        if (rn == null) continue;
        if (tCAs[rn] && mCAs[rn]) {
            tCoords.push(tCAs[rn]);
            mCoords.push(mCAs[rn]);
        }
    }
    if (tCoords.length < 3) return mobilePdb;
    const { R, t } = kabschAlign(tCoords, mCoords);
    return applyTransformPdb(mobilePdb, R, t);
}

function kabschAlign(targetCoords, mobileCoords) {
    const n = targetCoords.length;
    let ct = [0,0,0], cm = [0,0,0];
    for (let i = 0; i < n; i++) {
        for (let j = 0; j < 3; j++) { ct[j] += targetCoords[i][j]; cm[j] += mobileCoords[i][j]; }
    }
    ct = ct.map(v => v/n);
    cm = cm.map(v => v/n);
    const P = mobileCoords.map(c => [c[0]-cm[0], c[1]-cm[1], c[2]-cm[2]]);
    const Q = targetCoords.map(c => [c[0]-ct[0], c[1]-ct[1], c[2]-ct[2]]);
    const H = [[0,0,0],[0,0,0],[0,0,0]];
    for (let i = 0; i < n; i++)
        for (let r = 0; r < 3; r++)
            for (let c = 0; c < 3; c++)
                H[r][c] += P[i][r] * Q[i][c];
    const { U, V } = svd3x3(H);
    let R = matMul3(V, transpose3(U));
    if (det3(R) < 0) {
        V[0][2] *= -1; V[1][2] *= -1; V[2][2] *= -1;
        R = matMul3(V, transpose3(U));
    }
    const Rcm = matVec3(R, cm);
    const t = [ct[0]-Rcm[0], ct[1]-Rcm[1], ct[2]-Rcm[2]];
    return { R, t };
}

function transpose3(M) { return [[M[0][0],M[1][0],M[2][0]],[M[0][1],M[1][1],M[2][1]],[M[0][2],M[1][2],M[2][2]]]; }
function matMul3(A,B) {
    const C = [[0,0,0],[0,0,0],[0,0,0]];
    for (let i=0;i<3;i++) for (let j=0;j<3;j++) for (let k=0;k<3;k++) C[i][j]+=A[i][k]*B[k][j];
    return C;
}
function matVec3(M,v) { return [M[0][0]*v[0]+M[0][1]*v[1]+M[0][2]*v[2], M[1][0]*v[0]+M[1][1]*v[1]+M[1][2]*v[2], M[2][0]*v[0]+M[2][1]*v[1]+M[2][2]*v[2]]; }
function det3(M) { return M[0][0]*(M[1][1]*M[2][2]-M[1][2]*M[2][1])-M[0][1]*(M[1][0]*M[2][2]-M[1][2]*M[2][0])+M[0][2]*(M[1][0]*M[2][1]-M[1][1]*M[2][0]); }

function svd3x3(H) {
    const Ht = transpose3(H);
    const HtH = matMul3(Ht, H);
    const { eigenvalues, eigenvectors } = eigenDecomp3(HtH);
    const idx = [0,1,2].sort((a,b) => eigenvalues[b]-eigenvalues[a]);
    const V = [[0,0,0],[0,0,0],[0,0,0]];
    for (let i=0;i<3;i++) for (let j=0;j<3;j++) V[i][j] = eigenvectors[i][idx[j]];
    const HV = matMul3(H, V);
    const U = [[0,0,0],[0,0,0],[0,0,0]];
    for (let j=0;j<3;j++) {
        const sigma = Math.sqrt(Math.max(eigenvalues[idx[j]], 1e-10));
        for (let i=0;i<3;i++) U[i][j] = HV[i][j] / sigma;
    }
    return { U, V };
}

function eigenDecomp3(A) {
    const a = [[A[0][0],A[0][1],A[0][2]],[A[1][0],A[1][1],A[1][2]],[A[2][0],A[2][1],A[2][2]]];
    const v = [[1,0,0],[0,1,0],[0,0,1]];
    for (let iter = 0; iter < 50; iter++) {
        let p=0, q=1;
        if (Math.abs(a[0][2]) > Math.abs(a[p][q])) { p=0; q=2; }
        if (Math.abs(a[1][2]) > Math.abs(a[p][q])) { p=1; q=2; }
        if (Math.abs(a[p][q]) < 1e-12) break;
        const theta = 0.5 * Math.atan2(2*a[p][q], a[q][q]-a[p][p]);
        const c = Math.cos(theta), s = Math.sin(theta);
        const G = [[1,0,0],[0,1,0],[0,0,1]];
        G[p][p]=c; G[q][q]=c; G[p][q]=s; G[q][p]=-s;
        const tmp = matMul3(transpose3(G), matMul3(a, G));
        for (let i=0;i<3;i++) for (let j=0;j<3;j++) a[i][j]=tmp[i][j];
        const tmpv = matMul3(v, G);
        for (let i=0;i<3;i++) for (let j=0;j<3;j++) v[i][j]=tmpv[i][j];
    }
    return { eigenvalues: [a[0][0], a[1][1], a[2][2]], eigenvectors: v };
}

function applyTransformPdb(pdb, R, t) {
    return pdb.split('\n').map(line => {
        if (!line.startsWith('ATOM') && !line.startsWith('HETATM')) return line;
        if (line.length < 54) return line;
        const x = parseFloat(line.substring(30, 38));
        const y = parseFloat(line.substring(38, 46));
        const z = parseFloat(line.substring(46, 54));
        const nx = R[0][0]*x + R[0][1]*y + R[0][2]*z + t[0];
        const ny = R[1][0]*x + R[1][1]*y + R[1][2]*z + t[1];
        const nz = R[2][0]*x + R[2][1]*y + R[2][2]*z + t[2];
        return line.substring(0, 30) + nx.toFixed(3).padStart(8) + ny.toFixed(3).padStart(8) + nz.toFixed(3).padStart(8) + line.substring(54);
    }).join('\n');
}

// ========== Apply uncertainty coloring ==========
async function applyUncertaintyColoring(plugin) {
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

// ========== Load PDB into viewer ==========
async function loadPdb(viewer, pdbData, label) {
    const plugin = viewer.plugin;
    try {
        const rawData = await plugin.builders.data.rawData({ data: pdbData, label: label });
        const trajectory = await plugin.builders.structure.parseTrajectory(rawData, 'pdb');
        await plugin.builders.structure.hierarchy.applyPreset(trajectory, 'default', { showUnitcell: false });
        await applyUncertaintyColoring(plugin);
        console.log(label + ': loaded with uncertainty coloring');
    } catch(e) { console.error('Failed to load ' + label + ':', e); }
}

// ========== Camera sync (matches comparison HTML exactly) ==========
let cameraSyncPaused = false;
function setupCameraSync() {
    if (!plugin1 || !plugin2) return;
    let syncing = false;
    let activePlugin = null;
    const el1 = document.getElementById('viewer1');
    const el2 = document.getElementById('viewer2');
    el1.addEventListener('pointerenter', () => { activePlugin = plugin1; });
    el2.addEventListener('pointerenter', () => { activePlugin = plugin2; });
    el1.addEventListener('pointerleave', () => { if (activePlugin === plugin1) activePlugin = null; });
    el2.addEventListener('pointerleave', () => { if (activePlugin === plugin2) activePlugin = null; });

    function syncCamera(src, tgt) {
        if (syncing || cameraSyncPaused) return;
        syncing = true;
        try {
            const snap = src.canvas3d.camera.getSnapshot();
            tgt.canvas3d.camera.setState({
                position: new Float64Array(snap.position),
                target: new Float64Array(snap.target),
                up: new Float64Array(snap.up),
            });
            tgt.canvas3d.requestDraw();
        } finally { syncing = false; }
    }

    plugin1.canvas3d.didDraw.subscribe(() => { if (activePlugin === plugin1) syncCamera(plugin1, plugin2); });
    plugin2.canvas3d.didDraw.subscribe(() => { if (activePlugin === plugin2) syncCamera(plugin2, plugin1); });
}

// ========== Focus on residue ==========
function buildResidueLoci(plugin, resNum) {
    const structs = plugin.managers.structure.hierarchy.current.structures;
    if (!structs.length) return null;
    const struct = structs[0].cell.obj?.data;
    if (!struct) return null;
    const matched = [];
    for (const unit of struct.units) {
        if (!unit.conformation) continue;
        const indices = [];
        for (let eIdx = 0; eIdx < unit.elements.length; eIdx++) {
            const eI = unit.elements[eIdx];
            const rI = unit.residueIndex[eI];
            const rn = unit.model.atomicHierarchy.residues.auth_seq_id.value(rI);
            if (rn === resNum) {
                indices.push(eIdx);
            }
        }
        if (indices.length > 0) {
            matched.push({ unit, indices: Int32Array.from(indices) });
        }
    }
    if (matched.length === 0) return null;
    return { kind: 'element-loci', structure: struct, elements: matched };
}

// ========== Focus sync (matches comparison HTML pattern) ==========
function setupFocusSync() {
    let syncingFocus = false;

    // Both structures share the same PDB resnums at shared positions
    const sharedResnums = new Set(POSITIONS.map(p => p.rn));

    function getResNumFromLabel(label) {
        if (!label) return null;
        const m = label.match(/\b(\d+)\b/);
        return m ? parseInt(m[1]) : null;
    }

    function focusAndSelect(plugin, resNum) {
        const loci = buildResidueLoci(plugin, resNum);
        if (!loci) return;
        cameraSyncPaused = true;
        plugin.managers.structure.focus.setFromLoci(loci);
        setTimeout(() => {
            plugin.managers.camera.focusLoci(loci);
            setTimeout(() => { cameraSyncPaused = false; }, 300);
        }, 50);
    }

    plugin1.managers.structure.focus.behaviors.current.subscribe(val => {
        if (syncingFocus || !val) return;
        const resNum = getResNumFromLabel(val.label);
        if (resNum === null || !sharedResnums.has(resNum)) return;
        syncingFocus = true;
        try { focusAndSelect(plugin2, resNum); } finally { syncingFocus = false; }
    });

    plugin2.managers.structure.focus.behaviors.current.subscribe(val => {
        if (syncingFocus || !val) return;
        const resNum = getResNumFromLabel(val.label);
        if (resNum === null || !sharedResnums.has(resNum)) return;
        syncingFocus = true;
        try { focusAndSelect(plugin1, resNum); } finally { syncingFocus = false; }
    });
}

// ========== Focus from table click ==========
function focusMutation(row) {
    cameraSyncPaused = true;
    const rn = row.rn;
    if (rn != null && plugin1) {
        const loci1 = buildResidueLoci(plugin1, rn);
        if (loci1) {
            plugin1.managers.structure.focus.setFromLoci(loci1);
            plugin1.managers.camera.focusLoci(loci1);
        }
    }
    if (rn != null && plugin2) {
        const loci2 = buildResidueLoci(plugin2, rn);
        if (loci2) {
            plugin2.managers.structure.focus.setFromLoci(loci2);
            plugin2.managers.camera.focusLoci(loci2);
        }
    }
    // Sync cameras after focus settles
    setTimeout(() => {
        if (plugin1 && plugin2) {
            const snap = plugin1.canvas3d.camera.getSnapshot();
            plugin2.canvas3d.camera.setState({
                position: new Float64Array(snap.position),
                target: new Float64Array(snap.target),
                up: new Float64Array(snap.up),
            });
            plugin2.canvas3d.requestDraw();
        }
        setTimeout(() => { cameraSyncPaused = false; }, 200);
    }, 300);
}

// ========== Table ==========
let sortCol = 'abs_dd', sortAsc = false;
let selectedRow = null;

function renderTable() {
    const search = document.getElementById('searchBox').value.toLowerCase();
    const sameWt = document.getElementById('sameWtOnly').checked;
    const maxDD = Math.max(...TABLE_DATA.map(r => r.abs_dd));

    let filtered = TABLE_DATA.filter(r => {
        if (sameWt && r.wt1 !== r.wt2) return false;
        if (search) {
            const s = `${r.rn} ${r.ch} ${r.wt1} ${r.wt2} ${r.mut}`.toLowerCase();
            if (!s.includes(search)) return false;
        }
        return true;
    });

    filtered.sort((a, b) => {
        let va = a[sortCol], vb = b[sortCol];
        if (typeof va === 'string') return sortAsc ? va.localeCompare(vb) : vb.localeCompare(va);
        return sortAsc ? va - vb : vb - va;
    });

    const tbody = document.getElementById('tableBody');
    tbody.innerHTML = filtered.map((r, i) => {
        const ddClass = r.dd > 0 ? 'dd-pos' : 'dd-neg';
        const barWidth = Math.min(60, (r.abs_dd / maxDD) * 60);
        const barColor = r.dd > 0 ? '#ef9a9a' : '#90caf9';
        const chainClass = 'chain-' + r.ch;
        return `<tr data-idx="${i}" onclick="onRowClick(${i})" class="${selectedRow===i?'selected':''}">
            <td>${r.rn}</td>
            <td><span class="chain-tag ${chainClass}">${r.ch}</span></td>
            <td>${r.wt1}</td>
            <td>${r.wt2}</td>
            <td><b>${r.mut}</b></td>
            <td>${r.d1.toFixed(2)}</td>
            <td>${r.d2.toFixed(2)}</td>
            <td class="${ddClass}">${r.dd > 0 ? '+' : ''}${r.dd.toFixed(2)}</td>
            <td>${r.abs_dd.toFixed(2)}</td>
            <td><div class="dd-bar" style="width:${barWidth}px;background:${barColor}"></div></td>
        </tr>`;
    }).join('');

    document.getElementById('rowCount').textContent = `Showing ${filtered.length} of ${TABLE_DATA.length}`;
    window._filteredData = filtered;
}

function onRowClick(idx) {
    selectedRow = idx;
    const row = window._filteredData[idx];
    if (row) focusMutation(row);
    renderTable();
}

// Sort on header click
document.querySelectorAll('th[data-col]').forEach(th => {
    th.addEventListener('click', () => {
        const col = th.dataset.col;
        if (sortCol === col) sortAsc = !sortAsc;
        else { sortCol = col; sortAsc = false; }
        renderTable();
    });
});

document.getElementById('searchBox').addEventListener('input', () => { selectedRow = null; renderTable(); });
document.getElementById('sameWtOnly').addEventListener('change', () => { selectedRow = null; renderTable(); });

// ========== Init ==========
renderTable();
initViewers().catch(err => console.error('Viewer init error:', err));
</script>
</body>
</html>"""


if __name__ == '__main__':
    main()
