#!/usr/bin/env python3
# -*- coding: utf-8 -*-
###############################################################################
# Downstream Integration: MetaCHIP HGT + Virus-Carried Genes
# Identifies viral transduction events and summarizes taxonomic transfers
###############################################################################
import re, os
from collections import defaultdict

# ============================ Configuration ============================
METACHIP_HGT = "D:/onedrive/zzy/isme_revised/2.host/HGTs.txt"#output of metachip
TAXON        = "D:/onedrive/zzy/isme_revised/2.host/lakehost_taxon.tsv"
CARRIED      = "D:/onedrive/zzy/isme_revised/2.host/virus_carried_genes_final.tsv"
CLSTR        = "D:/onedrive/zzy/isme_revised/2.host/host_catalog.faa.clstr"
VHOST        = "D:/onedrive/zzy/isme_revised/2.host/votu2mag.txt"
OUTDIR       = "D:/onedrive/zzy/isme_revised/2.host/virus_HGT"
# =========================================================================

os.makedirs(OUTDIR, exist_ok=True)
LEVEL_NAME  = {'p':'cross-phylum','c':'cross-class','o':'cross-order','f':'cross-family','g':'cross-genus'}
LEVEL_ORDER = ['cross-phylum','cross-class','cross-order','cross-family','cross-genus','within-genus/unresolved']

def gene_to_mag(g):
    """Extract MAG ID from host gene name by removing the trailing number."""
    return re.sub(r'_\d+$', '', g)

# ---------- 0. Load valid vOTUs and MAGs ----------
votu_hosts = defaultdict(set)
valid_votu = set()
valid_mag  = set()

with open(VHOST) as fh:
    for line in fh:
        p = re.split(r'\t|\s+', line.rstrip('\n'), maxsplit=1)
        if len(p) < 2: continue
        v, m = p[0].strip(), p[1].strip()
        if not v or not m: continue
        votu_hosts[v].add(m)
        valid_votu.add(v)
        valid_mag.add(m)
print(f"Loaded: {len(valid_votu)} valid vOTUs, {len(valid_mag)} valid MAGs")

def votu_of_gene(carried_gene):
    """Map carried gene ID to vOTU ID and align with valid datasets."""
    cand = re.sub(r'_\d+$', '', carried_gene)
    if cand in valid_votu:
        return cand
    parts = carried_gene.split('_')
    for k in range(len(parts)-1, 0, -1):
        trial = '_'.join(parts[:k])
        if trial in valid_votu:
            return trial
    return cand

# ---------- 1. Parse GTDB taxonomy ----------
taxo = {}
with open(TAXON) as fh:
    next(fh, None)
    for line in fh:
        p = line.rstrip('\n').split('\t')
        if len(p) < 2: continue
        d = {}
        for tok in p[1].split(';'):
            tok = tok.strip()
            if len(tok) >= 3 and tok[1:3] == '__':
                d[tok[0]] = tok[3:]
        taxo[p[0]] = d

def distance_level(a, b):
    ta, tb = taxo.get(a, {}), taxo.get(b, {})
    for r in ['p','c','o','f','g']:
        va, vb = ta.get(r,''), tb.get(r,'')
        if va and vb and va != vb:
            return LEVEL_NAME[r]
    return 'within-genus/unresolved'

def phylum(m):
    return taxo.get(m, {}).get('p', 'NA')

# ---------- 2. Parse CD-HIT clusters (Member -> Representative) ----------
member2rep = {}
cur_members, cur_rep = [], None
with open(CLSTR) as fh:
    for line in fh:
        if line.startswith('>Cluster'):
            if cur_rep is not None:
                for m in cur_members: member2rep[m] = cur_rep
            cur_members, cur_rep = [], None
        else:
            mm = re.search(r'>(\S+?)\.\.\.', line)
            if not mm: continue
            gid = mm.group(1)
            cur_members.append(gid)
            if line.rstrip().endswith('*'):
                cur_rep = gid
    if cur_rep is not None:
        for m in cur_members: member2rep[m] = cur_rep

# ---------- 3. Map host cluster representatives to carrying vOTUs ----------
rep2votus = defaultdict(set)
with open(CARRIED) as fh:
    for line in fh:
        p = re.split(r'\t|\s+', line.rstrip('\n'))
        if len(p) < 2: continue
        carried_gene, host_rep = p[0], p[1]
        v = votu_of_gene(carried_gene)
        if v in valid_votu:
            rep2votus[host_rep].add(v)

# ---------- 4. Parse MetaCHIP HGT events and identify transductions ----------
with open(METACHIP_HGT) as fh:
    lines = fh.read().splitlines()
header = lines[0].split('\t')

def find_col(prefix):
    for i, h in enumerate(header):
        if h.lower().startswith(prefix.lower()): return i
    return None

i_g1, i_g2 = find_col('Gene_1'), find_col('Gene_2')
i_dir      = find_col('direction')
if i_g1 is None or i_g2 is None:
    raise SystemExit("ERROR: Gene_1 or Gene_2 columns not found")

def parse_direction(s):
    if not s: return ''
    return re.sub(r'\([^)]*\)', '', s).strip()

events = []
unmapped = 0
skipped_mag = 0

for line in lines[1:]:
    if not line.strip(): continue
    p = line.split('\t')
    g1, g2 = p[i_g1], p[i_g2]
    A, B = gene_to_mag(g1), gene_to_mag(g2)

    if A not in valid_mag or B not in valid_mag:
        skipped_mag += 1
        continue

    direction = parse_direction(p[i_dir] if (i_dir is not None and i_dir < len(p)) else '')
    level = distance_level(A, B)

    reps = set()
    if g1 in member2rep: reps.add(member2rep[g1])
    if g2 in member2rep: reps.add(member2rep[g2])
    if not reps: unmapped += 1

    carrying = set()
    for r in reps: carrying |= rep2votus.get(r, set())
    support = sorted(v for v in carrying if {A, B} <= votu_hosts.get(v, set()))
    transduced = len(support) > 0

    events.append(dict(gene_1=g1, gene_2=g2, MAG_A=A, MAG_B=B, direction=direction,
                       level=level, transduced=transduced, vOTUs=";".join(support)))

print(f"Read {len(events)} HGT events (Skipped MAGs: {skipped_mag}; Unmapped clusters: {unmapped})")

# ---------- 5. Output event-level table ----------
with open(os.path.join(OUTDIR, "hgt_events_annotated.tsv"), "w") as out:
    out.write("Gene_1\tGene_2\tMAG_A\tMAG_B\tDirection\tDistance_level\tTransduced\tSupporting_vOTUs\n")
    for e in events:
        out.write("\t".join([e['gene_1'],e['gene_2'],e['MAG_A'],e['MAG_B'],e['direction'],
                              e['level'],"yes" if e['transduced'] else "no",e['vOTUs']])+"\n")

# ---------- 6. Aggregate to unique MAG pairs ----------
pairs = {}
for e in events:
    key = frozenset((e['MAG_A'], e['MAG_B']))
    if key not in pairs:
        a, b = e['MAG_A'], e['MAG_B']
        pairs[key] = dict(MAG_A=a, MAG_B=b, level=e['level'],
                          phylum_A=phylum(a), phylum_B=phylum(b),
                          n_events=0, n_transduced_events=0, vOTUs=set())
    pairs[key]['n_events'] += 1
    if e['transduced']:
        pairs[key]['n_transduced_events'] += 1
        if e['vOTUs']: pairs[key]['vOTUs'].update(e['vOTUs'].split(';'))

with open(os.path.join(OUTDIR, "mag_pairs.tsv"), "w") as out:
    out.write("MAG_A\tMAG_B\tDistance_level\tPhylum_A\tPhylum_B\tN_transfer_events\t"
              "N_transduced_events\tTransduced_any\tSupporting_vOTUs\n")
    for v in pairs.values():
        out.write("\t".join([v['MAG_A'],v['MAG_B'],v['level'],v['phylum_A'],v['phylum_B'],
                             str(v['n_events']),str(v['n_transduced_events']),
                             "yes" if v['n_transduced_events']>0 else "no",
                             ";".join(sorted(v['vOTUs']))])+"\n")

# ---------- 7. Summarize by taxonomic distance level ----------
total_pairs = len(pairs)
by_level = defaultdict(lambda: [0, 0])
for v in pairs.values():
    by_level[v['level']][0] += 1
    if v['n_transduced_events'] > 0: by_level[v['level']][1] += 1

with open(os.path.join(OUTDIR, "level_summary.tsv"), "w") as out:
    out.write("Distance_level\tN_pairs\tPct_pairs(%)\tN_transduced_pairs\tPct_transduced_within_level(%)\n")
    for lv in LEVEL_ORDER:
        if lv not in by_level: continue
        n, nt = by_level[lv]
        pct = 100.0*n/total_pairs if total_pairs else 0
        ptv = 100.0*nt/n if n else 0
        out.write(f"{lv}\t{n}\t{pct:.2f}\t{nt}\t{ptv:.2f}\n")

# ---------- 8. Output transduced MAG pairs ----------
with open(os.path.join(OUTDIR, "transduction_mag_pairs.tsv"), "w") as out:
    out.write("MAG_A\tMAG_B\tDistance_level\tPhylum_A\tPhylum_B\tN_transduced_events\tSupporting_vOTUs\n")
    for v in pairs.values():
        if v['n_transduced_events'] > 0:
            out.write("\t".join([v['MAG_A'],v['MAG_B'],v['level'],v['phylum_A'],v['phylum_B'],
                                 str(v['n_transduced_events']),";".join(sorted(v['vOTUs']))])+"\n")

# ---------- Summary ----------
n_events = len(events)
n_trans_events = sum(1 for e in events if e['transduced'])
n_trans_pairs = sum(1 for v in pairs.values() if v['n_transduced_events']>0)

print("============== Summary ==============")
print(f"Valid HGT events          : {n_events}")
if n_events:
    print(f"Viral transduction events : {n_trans_events} ({100.0*n_trans_events/n_events:.1f}%)")
print(f"Transferred MAG pairs     : {total_pairs}")
if total_pairs:
    print(f"Transduced MAG pairs      : {n_trans_pairs} ({100.0*n_trans_pairs/total_pairs:.1f}%)")
print("By distance level (MAG pairs):")

for lv in LEVEL_ORDER:
    if lv in by_level:
        n, nt = by_level[lv]
        print(f"  {lv:<26} {n:>6} pairs ({100.0*n/total_pairs:.1f}%), transduced: {nt}")
print("Outputs generated: hgt_events_annotated.tsv / mag_pairs.tsv / level_summary.tsv / transduction_mag_pairs.tsv")