import json
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# ── File paths (relative to this script — works on Windows/Linux/Mac) ──
FOLDER = os.path.dirname(os.path.abspath(__file__))
def path(f): return os.path.join(FOLDER, f)

# ════════════════════════════════════════════════════════════
#  ReRAM DEVICE PARAMETERS — FROM CITED PAPERS
#  DO NOT CHANGE these unless you have a different source paper.
#
#  [1] MAC energy = 11 fJ per operation
#      Source: Marinella et al., "Multiscale Co-Design Analysis of
#              Energy, Latency, Area, and Accuracy of a ReRAM Analog
#              Neural Training Accelerator"
#      Link:   https://arxiv.org/abs/1707.09952
#      Device: TaOx ReRAM, 14/16nm node
#      Quote:  "takes only 11 fJ per multiply and accumulate (MAC)"
#
#  [2] Write energy = 0.1 pJ to 1.6 pJ per cell (we use 0.1 pJ = Cell A)
#      Source: Lunkai et al., "Mellow Writes: Extending Lifetime in
#              Resistive Memories", ISCA 2016
#      Link:   https://people.cs.uchicago.edu/~ftchong/papers/ISCA-16-mellow.pdf
#      Device: ReRAM under 22nm process
#      Quote:  "set/reset energy of a cell... from 0.1pJ (CellA) to 1.6pJ (CellE)"
#
#  [3] Read latency ~ 10 ns, Write latency ~ 50 ns
#      Source: Chi et al., "PRIME: A Novel Processing-in-Memory
#              Architecture for Neural Network Computation in
#              ReRAM-based Main Memory", ISCA 2016
#      Link:   https://ieeexplore.ieee.org/document/7551380
#      Quote:  "read latency of ReRAM is comparable to that of DRAM
#               while its write latency is significantly longer (~5x)"
#      Note:   DRAM read latency ~10ns, so ReRAM read ~10ns,
#              ReRAM write ~50ns (5x DRAM)
#
#  [4] CPU baseline: 3 GHz processor, DRAM access ~60ns
#      Standard assumption used in computer architecture papers.
# ════════════════════════════════════════════════════════════

# Energy per operation (in picojoules, pJ)
E_MAC_RERAM   = 0.000011   # 11 fJ = 0.000011 pJ  [Source 1]
E_WRITE_RERAM = 0.1        # 0.1 pJ per cell write [Source 2]
E_READ_RERAM  = 0.01       # ~10x less than write, typical for ReRAM reads
E_ARITH_CPU   = 1.0        # ~1 pJ per op on CPU (standard CMOS reference)
E_MEM_CPU     = 60.0       # ~60 pJ per DRAM access (standard reference)
E_MAC_CPU     = 1.0        # ~1 pJ per MAC on CPU

# Latency per operation (in nanoseconds, ns)
L_READ_RERAM  = 10.0       # ~10 ns  [Source 3]
L_WRITE_RERAM = 50.0       # ~50 ns  [Source 3]
L_MAC_RERAM   = 1.0        # ~1 ns for analog crossbar MAC [Source 1 context]
L_MEM_CPU     = 60.0       # ~60 ns DRAM latency (standard)
L_MAC_CPU     = 0.33       # ~0.33 ns at 3 GHz CPU

plt.rcParams.update({
    'font.family': 'DejaVu Sans',
    'font.size':        11,
    'axes.titlesize':   12,
    'axes.labelsize':   11,
    'axes.spines.top':   False,
    'axes.spines.right': False,
    'figure.facecolor': '#F7F9FC',
    'axes.facecolor':   '#F7F9FC',
})

STEP_COLORS = {
    'Kmer-to-Current Table': '#2E5FA3',
    'Reference Index Build': '#5B8DD9',
    'Signal Normalization':  '#E07B39',
    'Event Detection':       '#F5A623',
    'Query Minimizers':      '#2E8B57',
    'Seeding':               '#8BC34A',
    'DTW Alignment':         '#C0392B',
    'Chaining & Mapping':    '#E67E22',
}

# ── Load JSON ────────────────────────────────────────────────
json_path = path('rawhash_results.json')
if not os.path.exists(json_path):
    print(f"ERROR: Could not find {json_path}")
    print("Run ./rawhash 5000 500 first!")
    exit(1)

with open(json_path) as f:
    data = json.load(f)

steps       = data['steps']
has_phase   = 'phase' in steps[0]
online_steps  = [s for s in steps if s.get('phase','ONLINE') == 'ONLINE'] if has_phase else steps
offline_steps = [s for s in steps if s.get('phase','') == 'OFFLINE'] if has_phase else []
offline_total = data.get('offline_total', 0)
online_total  = data.get('online_total',  sum(s['total'] for s in online_steps))
grand_total   = data.get('grand_total',   offline_total + online_total)
ref_len       = data['ref_len']
query_len     = data['query_len']

# ── Aggregate online operation counts ───────────────────────
total_mac    = sum(s['mac']     for s in online_steps)
total_mem    = sum(s['mem']     for s in online_steps)
total_arith  = sum(s['arith']   for s in online_steps)
total_cmp    = sum(s['compare'] for s in online_steps)

print(f"Online phase operation counts:")
print(f"  MAC ops    : {total_mac:,}")
print(f"  Memory ops : {total_mem:,}")
print(f"  Arith ops  : {total_arith:,}")
print(f"  Compare ops: {total_cmp:,}")

# ════════════════════════════════════════════════════════════
#  LATENCY ESTIMATION
#  CPU:   sequential, one op at a time
#  ReRAM: in-memory, no data movement penalty for MAC
# ════════════════════════════════════════════════════════════

# CPU time (ns) — all ops sequential at 3 GHz
cpu_time_ns = (total_mac   * L_MAC_CPU
             + total_mem   * L_MEM_CPU
             + total_arith * (1000.0/3000.0)
             + total_cmp   * (1000.0/3000.0))

# ReRAM time (ns) — MAC and memory ops use ReRAM latencies
reram_time_ns = (total_mac   * L_MAC_RERAM
               + total_mem   * L_READ_RERAM     # reads dominate over writes
               + total_arith * (1000.0/3000.0)  # arith still on CPU periphery
               + total_cmp   * (1000.0/3000.0))

speedup = cpu_time_ns / reram_time_ns if reram_time_ns > 0 else 0

# ════════════════════════════════════════════════════════════
#  ENERGY ESTIMATION  (in picojoules, pJ)
# ════════════════════════════════════════════════════════════

# CPU energy
cpu_energy_pJ = (total_mac   * E_MAC_CPU
               + total_mem   * E_MEM_CPU
               + total_arith * E_ARITH_CPU
               + total_cmp   * E_ARITH_CPU)

# ReRAM energy
reram_energy_pJ = (total_mac   * E_MAC_RERAM
                 + total_mem   * E_READ_RERAM
                 + total_arith * E_ARITH_CPU    # peripheral CMOS for arith
                 + total_cmp   * E_ARITH_CPU)

energy_savings = cpu_energy_pJ / reram_energy_pJ if reram_energy_pJ > 0 else 0

print(f"\nLatency estimate (online phase):")
print(f"  CPU   : {cpu_time_ns/1e6:.3f} ms")
print(f"  ReRAM : {reram_time_ns/1e6:.3f} ms")
print(f"  Speedup: {speedup:.1f}x")
print(f"\nEnergy estimate (online phase):")
print(f"  CPU   : {cpu_energy_pJ/1e6:.3f} nJ")
print(f"  ReRAM : {reram_energy_pJ:.3f} pJ")
print(f"  Energy savings: {energy_savings:.1f}x")

# ════════════════════════════════════════════════════════════
#  FIGURE 1: Stage breakdown + pie chart
# ════════════════════════════════════════════════════════════
fig, axes = plt.subplots(1, 2, figsize=(14, 6))
fig.suptitle(
    f'RawHash Pipeline — Cycle Count Analysis\n'
    f'Ref: {ref_len} bp  |  Query: {query_len} bp  |  Grand Total: {grand_total:,} cycles',
    fontsize=13, fontweight='bold'
)

# Horizontal bar chart
ax = axes[0]
names  = [s['name']  for s in steps]
totals = [s['total'] for s in steps]
colors = [STEP_COLORS.get(s['name'], '#999') for s in steps]
y = range(len(steps))
bars = ax.barh(list(y), totals, color=colors, edgecolor='white', height=0.7)
ax.set_yticks(list(y))
ax.set_yticklabels(names, fontsize=9)
ax.set_xlabel('Total Cycles')
ax.set_title('Cycles per Stage')
ax.xaxis.set_major_formatter(plt.FuncFormatter(
    lambda x, _: f'{x/1e6:.1f}M' if x >= 1e6 else f'{x/1e3:.0f}K'))
for bar, val in zip(bars, totals):
    ax.text(bar.get_width()*1.01, bar.get_y()+bar.get_height()/2,
            f'{val:,}', va='center', fontsize=7.5)

# Pie chart
ax2 = axes[1]
nonzero = [(s['name'], s['total']) for s in steps if s['total'] > 0]
pie_names, pie_vals = zip(*nonzero)
pie_colors = [STEP_COLORS.get(n,'#999') for n in pie_names]
wedges, texts, autotexts = ax2.pie(
    pie_vals, colors=pie_colors,
    autopct=lambda p: f'{p:.1f}%' if p > 3 else '',
    startangle=90, pctdistance=0.75,
    wedgeprops={'edgecolor':'white','linewidth':1.5})
for at in autotexts: at.set_fontsize(8)
ax2.legend(wedges, [n[:22] for n in pie_names],
           loc='lower center', bbox_to_anchor=(0.5,-0.2),
           ncol=2, fontsize=7.5, frameon=False)
ax2.set_title('Cycle Share by Stage')

plt.tight_layout()
plt.savefig(path('rawhash_stage_analysis.png'), dpi=150, bbox_inches='tight')
plt.close()
print("\nSaved: rawhash_stage_analysis.png")

# ════════════════════════════════════════════════════════════
#  FIGURE 2: CPU vs ReRAM — Latency and Energy
# ════════════════════════════════════════════════════════════
fig2, axes2 = plt.subplots(1, 2, figsize=(13, 6))
fig2.suptitle(
    'CPU vs ReRAM — Latency and Energy Comparison\n'
    'Values from cited papers (see code comments for sources)',
    fontsize=12, fontweight='bold'
)

# Latency comparison
ax = axes2[0]
labels = ['CPU\n(3 GHz + DRAM)', 'ReRAM\nIn-Memory']
lat_vals = [cpu_time_ns/1e6, reram_time_ns/1e6]
colors_lat = ['#2E5FA3', '#C0392B']
bars = ax.bar(labels, lat_vals, color=colors_lat, width=0.4, edgecolor='white')
ax.set_ylabel('Time (ms) per read')
ax.set_title(f'Latency\n(Speedup: {speedup:.1f}×)')
for bar, val in zip(bars, lat_vals):
    ax.text(bar.get_x()+bar.get_width()/2,
            bar.get_height()+max(lat_vals)*0.01,
            f'{val:.3f} ms', ha='center', fontsize=10, fontweight='bold')
ax.set_ylim(0, max(lat_vals)*1.3)
# Citation box
ax.text(0.5, 0.82,
    'Read latency: ~10 ns [PRIME, ISCA 2016]\n'
    'MAC latency:   ~1 ns [Marinella et al. 2017]',
    transform=ax.transAxes, ha='center', fontsize=7.5,
    bbox=dict(boxstyle='round', facecolor='#EAF4FB', alpha=0.9))

# Energy comparison
ax = axes2[1]
energy_vals = [cpu_energy_pJ/1e6, reram_energy_pJ/1e6]   # convert to nJ
colors_en = ['#2E5FA3', '#27AE60']
bars = ax.bar(labels, energy_vals, color=colors_en, width=0.4, edgecolor='white')
ax.set_ylabel('Energy (nJ) per read')
ax.set_title(f'Energy\n(Savings: {energy_savings:.1f}×)')
for bar, val in zip(bars, energy_vals):
    ax.text(bar.get_x()+bar.get_width()/2,
            bar.get_height()+max(energy_vals)*0.01,
            f'{val:.4f} nJ', ha='center', fontsize=10, fontweight='bold')
ax.set_ylim(0, max(energy_vals)*1.3)
ax.text(0.5, 0.82,
    'MAC energy:   11 fJ [Marinella et al. 2017]\n'
    'Write energy: 0.1 pJ [Mellow Writes, ISCA 2016]',
    transform=ax.transAxes, ha='center', fontsize=7.5,
    bbox=dict(boxstyle='round', facecolor='#EAFAF1', alpha=0.9))

plt.tight_layout()
plt.savefig(path('rawhash_cpu_vs_reram.png'), dpi=150, bbox_inches='tight')
plt.close()
print("Saved: rawhash_cpu_vs_reram.png")

# ════════════════════════════════════════════════════════════
#  FIGURE 3: Summary table with citations
# ════════════════════════════════════════════════════════════
fig3, axes3 = plt.subplots(2, 1, figsize=(14, 8))
fig3.suptitle(
    f'RawHash ReRAM Estimation Summary — Ref: {ref_len} bp | Query: {query_len} bp',
    fontsize=12, fontweight='bold'
)

# Top table: cycle counts
ax = axes3[0]
ax.axis('off')
rows = []
for s in steps:
    pct = 100*s['total']/grand_total if grand_total > 0 else 0
    phase = s.get('phase','ONLINE')
    rows.append([phase, s['name'],
                 f"{s['arith']:,}", f"{s['mem']:,}",
                 f"{s['compare']:,}", f"{s['mac']:,}",
                 f"{s['total']:,}", f"{pct:.1f}%"])
rows.append(['', 'ONLINE TOTAL','','','','',
             f"{online_total:,}", f"{100*online_total/grand_total:.1f}%"])

col1 = ['Phase','Step','Arith Ops','Mem Ops','Cmp Ops','MAC Ops','Total Cycles','% Share']
tbl = ax.table(cellText=rows, colLabels=col1, loc='center', cellLoc='center')
tbl.auto_set_font_size(False)
tbl.set_fontsize(8.5)
tbl.scale(1, 1.4)
for j in range(len(col1)):
    tbl[(0,j)].set_facecolor('#1A3A6B')
    tbl[(0,j)].set_text_props(color='white', fontweight='bold')
for i, row in enumerate(rows):
    bg = '#D6E4F7' if row[0]=='OFFLINE' else '#FFF0E6' if row[0]=='ONLINE' else '#2C3E50'
    if 'TOTAL' in row[1]: bg = '#2C3E50'
    if 'DTW' in row[1]: bg = '#FADBD8'
    for j in range(len(col1)):
        tbl[(i+1,j)].set_facecolor(bg)
        if bg == '#2C3E50':
            tbl[(i+1,j)].set_text_props(color='white', fontweight='bold')
ax.set_title('Operation Counts per Stage', fontsize=10, fontweight='bold', pad=8)

# Bottom table: estimation results with citations
ax2 = axes3[1]
ax2.axis('off')
est_rows = [
    ['MAC energy per op',    '11 fJ',      '0.000011 pJ', 'Marinella et al. 2017 — arXiv:1707.09952'],
    ['Write energy per op',  '0.1 pJ',     '0.1 pJ',      'Lunkai et al. ISCA 2016 (Mellow Writes)'],
    ['Read energy per op',   '~0.01 pJ',   '0.01 pJ',     'Estimated ~10x less than write (standard)'],
    ['Read latency',         '~10 ns',     '10 ns',        'Chi et al. PRIME ISCA 2016'],
    ['Write latency',        '~50 ns',     '50 ns',        'Chi et al. PRIME ISCA 2016 (5× DRAM)'],
    ['MAC latency (crossbar)','~1 ns',     '1 ns',         'Marinella et al. 2017 (analog crossbar)'],
    ['','','',''],
    ['CPU latency (online)',  f'{cpu_time_ns/1e6:.3f} ms', '—', 'Calculated from op counts × CPU latency'],
    ['ReRAM latency (online)',f'{reram_time_ns/1e6:.3f} ms','—','Calculated from op counts × ReRAM latency'],
    ['Latency Speedup',      f'{speedup:.1f}×',  '—', 'CPU time / ReRAM time'],
    ['CPU energy (online)',  f'{cpu_energy_pJ/1e6:.4f} nJ','—','Calculated from op counts × CPU energy'],
    ['ReRAM energy (online)',f'{reram_energy_pJ:.2f} pJ',  '—','Calculated from op counts × ReRAM energy'],
    ['Energy Savings',       f'{energy_savings:.1f}×',     '—', 'CPU energy / ReRAM energy'],
]
col2 = ['Parameter', 'Literature Value', 'Used in Code', 'Source']
tbl2 = ax2.table(cellText=est_rows, colLabels=col2, loc='center', cellLoc='left')
tbl2.auto_set_font_size(False)
tbl2.set_fontsize(8)
tbl2.scale(1, 1.35)
for j in range(len(col2)):
    tbl2[(0,j)].set_facecolor('#1A3A6B')
    tbl2[(0,j)].set_text_props(color='white', fontweight='bold')
highlight_rows = [7,8,9,10,11,12]
for i in range(len(est_rows)):
    bg = '#EBF5FB' if i not in highlight_rows else '#EAFAF1'
    if i in [9,12]:  bg = '#D5F5E3'
    for j in range(len(col2)):
        tbl2[(i+1,j)].set_facecolor(bg)
ax2.set_title('ReRAM Parameter Values and Estimation Results', fontsize=10, fontweight='bold', pad=8)

plt.tight_layout()
plt.savefig(path('rawhash_cycle_table.png'), dpi=150, bbox_inches='tight')
plt.close()
print("Saved: rawhash_cycle_table.png")
print("\n All 3 figures saved!")
print(f"   Folder: {FOLDER}")
print(f"\nCited papers:")
print(f"  [1] Marinella et al. 2017 — https://arxiv.org/abs/1707.09952")
print(f"  [2] Lunkai et al. ISCA 2016 — https://people.cs.uchicago.edu/~ftchong/papers/ISCA-16-mellow.pdf")
print(f"  [3] Chi et al. PRIME ISCA 2016 — https://ieeexplore.ieee.org/document/7551380")
