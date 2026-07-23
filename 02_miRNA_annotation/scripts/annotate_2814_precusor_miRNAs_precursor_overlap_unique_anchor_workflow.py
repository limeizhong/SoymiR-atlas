#!/usr/bin/env python3
from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
import re


SCRIPT_DIR = Path(__file__).resolve().parent
MODULE_DIR = SCRIPT_DIR.parent
INPUT_DIR = MODULE_DIR / "input"
RESULTS_DIR = MODULE_DIR / "results"
INTERMEDIATE_DIR = RESULTS_DIR / "intermediate"
MIN_PRECURSOR_OVERLAP_BP = 1
OPPOSITE_ARM_MIN_PRECURSOR_OVERLAP = 0.8
OPPOSITE_ARM_PRECURSOR_TOLERANCE_NT = 3
OPPOSITE_ARM_MIN_RELATIVE_SEPARATION = 0.45


@dataclass
class Row:
    idx: int
    file_line: int
    seq_id: str
    seq: str
    chr: str
    m_start: int
    m_end: int
    h_start: int
    h_end: int
    strand: str
    status: str = ""
    annotation: str = ""
    source: str = ""
    evidence: str = ""
    matched_mature: str = ""
    matched_precursor: str = ""
    matched_chr: str = ""
    matched_start: str = ""
    matched_end: str = ""
    matched_strand: str = ""
    distance: str = ""

    @property
    def key(self):
        return (
            self.seq_id,
            self.seq,
            self.chr,
            str(self.m_start),
            str(self.m_end),
            str(self.h_start),
            str(self.h_end),
            self.strand,
        )

    @property
    def precursor_mid(self) -> float:
        return (self.h_start + self.h_end) / 2


@dataclass
class Precursor:
    name: str
    accession: str
    family: str
    chr: str
    start: int
    end: int
    strand: str
    anchorable: bool = False

    @property
    def mid(self) -> float:
        return (self.start + self.end) / 2


@dataclass
class Hit:
    qseqid: str
    sseqid: str
    pident: float
    length: int
    qstart: int
    qend: int
    sstart: int
    send: int
    evalue: str
    bitscore: float
    qseq: str
    sseq: str
    stitle: str
    qlen: int
    slen: int

    @property
    def is_pident100(self) -> bool:
        return abs(self.pident - 100.0) < 1e-9

    @property
    def is_full_exact(self) -> bool:
        return self.is_pident100 and self.length == self.qlen == self.slen


def norm_seq(seq: str) -> str:
    return seq.upper().replace("U", "T")


def read_fasta(path: Path) -> dict[str, str]:
    seqs = {}
    current = None
    parts = []
    with path.open() as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if current is not None:
                    seqs[current] = norm_seq("".join(parts))
                current = line[1:].split()[0]
                parts = []
            else:
                parts.append(line.strip())
    if current is not None:
        seqs[current] = norm_seq("".join(parts))
    return seqs


def clean_chr(value: str) -> str:
    value = value.strip()
    value = re.sub(r"^chr", "", value, flags=re.I)
    value = re.sub(r"^Gma-Gm", "", value, flags=re.I)
    if value.isdigit():
        value = str(int(value))
    return value


def family_of_name(name: str) -> str:
    m = re.search(r"(miRNC|miRN|miR)(\d+)", name, flags=re.I)
    if not m:
        return ""
    kind = m.group(1).lower()
    if kind == "mirnc":
        prefix = "miRNC"
    elif kind == "mirn":
        prefix = "miRN"
    else:
        prefix = "miR"
    return f"{prefix}{m.group(2)}"


def family_label(annotation: str) -> str:
    name = annotation.split()[0]
    name = re.sub(r"\.v\d+$", "", name)
    m = re.match(r"^(.*?(?:miRNC|miRN|miR)\d+)", name, flags=re.I)
    return m.group(1) if m else name


def normalize_mature_annotation_name(name: str) -> str:
    """Normalize mature miRNA-style output names without changing precursor IDs."""
    if not name:
        return name
    name = re.sub(r"^Gma-MIRN", "Gma-miRN", name)
    name = re.sub(r"^Gma-MIR", "Gma-miR", name)
    name = re.sub(r"^Gma-miRN", "Gma-miRN", name)
    name = re.sub(r"^Gma-miR", "Gma-miR", name)
    name = re.sub(r"^gma-MIRN", "gma-miRN", name)
    name = re.sub(r"^gma-MIR", "gma-miR", name)
    name = re.sub(r"^gma-mirn", "gma-miRN", name, flags=re.I)
    name = re.sub(r"^gma-mir", "gma-miR", name, flags=re.I)
    return name


def normalize_mature_annotation_text(text: str) -> str:
    """Normalize mature miRNA-style names embedded in free-text evidence."""
    if not text:
        return text
    text = re.sub(r"\bGma-MIRN", "Gma-miRN", text)
    text = re.sub(r"\bGma-MIR", "Gma-miR", text)
    text = re.sub(r"\bgma-MIRN", "gma-miRN", text)
    text = re.sub(r"\bgma-MIR", "gma-miR", text)
    return text


def output_family_label(annotation: str) -> str:
    label = family_label(normalize_mature_annotation_name(annotation))
    if re.match(r"^gma_mirnc\d+$", label, flags=re.I):
        return re.sub(r"^gma_mirnc", "gma_miRNC", label, flags=re.I)
    label = re.sub(r"^Gma-", "gma-", label)
    label = re.sub(r"^gma-MIRN", "gma-miRN", label)
    label = re.sub(r"^gma-MIR", "gma-miR", label)
    label = re.sub(r"^gma-mirn", "gma-miRN", label)
    label = re.sub(r"^gma-mir", "gma-miR", label)
    return label


def read_count(seq_id: str) -> int:
    m = re.search(r"_R(\d+)_", seq_id)
    return int(m.group(1)) if m else 0


def precursor_from_mature(mature: str) -> str:
    base = re.sub(r"-(?:5p|3p)$", "", mature, flags=re.I)
    return re.sub(r"miR", "MIR", base, flags=re.I)


def mature_from_precursor(precursor: str) -> str:
    return re.sub(r"MIR", "miR", precursor, flags=re.I)


def suffix_of_member(name: str) -> str:
    name = re.sub(r"-(?:5p|3p)$", "", name, flags=re.I)
    m = re.search(r"(?:miRNC|miRN|miR)\d+([a-z]+)$", name, flags=re.I)
    return m.group(1).lower() if m else ""


def suffix_to_int(s: str) -> int:
    if not s:
        return 0
    n = 0
    for ch in s:
        n = n * 26 + (ord(ch) - 96)
    return n


def int_to_suffix(n: int) -> str:
    out = ""
    while n:
        n, r = divmod(n - 1, 26)
        out = chr(97 + r) + out
    return out or "a"


def next_family_name(prefix: str, family: str, counters: dict[str, int]) -> str:
    key = f"{prefix}:{family}"
    counters[key] += 1
    m = re.search(r"(miRNC|miRN|miR)(\d+)", family, flags=re.I)
    return f"{prefix}-{m.group(1)}{m.group(2)}{int_to_suffix(counters[key])}"


def load_rows(path: Path) -> list[Row]:
    rows = []
    with path.open() as f:
        header = f.readline().rstrip("\n").split("\t")
        idx = {name: i for i, name in enumerate(header)}
        for file_line, line in enumerate(f, 2):
            if not line.strip():
                continue
            r = line.rstrip("\n").split("\t")
            rows.append(
                Row(
                    idx=len(rows),
                    file_line=file_line,
                    seq_id=r[idx["Seq-ID"]],
                    seq=norm_seq(r[idx["Sequences_Mature"]]),
                    chr=clean_chr(r[idx["Chr"]]),
                    m_start=int(r[idx["M_start"]]),
                    m_end=int(r[idx["M_end"]]),
                    h_start=int(r[idx["H_start"]]),
                    h_end=int(r[idx["H_end"]]),
                    strand=r[idx["Strand"]],
                )
            )
    return rows


def load_variant_clusters(path: Path, rows_by_key: dict[tuple[str, ...], int]):
    by_row = {}
    members = defaultdict(list)
    if not path.exists():
        return by_row, members
    with path.open() as f:
        header = f.readline().rstrip("\n").split("\t")
        idx = {name: i for i, name in enumerate(header)}
        for line in f:
            if not line.strip():
                continue
            r = line.rstrip("\n").split("\t")
            key = (
                r[idx["Seq-ID"]],
                norm_seq(r[idx["Sequences_Mature"]]),
                clean_chr(r[idx["Chr"]]),
                r[idx["M_start"]],
                r[idx["M_end"]],
                r[idx["H_start"]],
                r[idx["H_end"]],
                r[idx["Strand"]],
            )
            if key in rows_by_key:
                i = rows_by_key[key]
                cid = r[idx["Variant_Cluster"]]
                by_row[i] = cid
                members[cid].append(i)
    return by_row, members


def load_cdhit_clstr(path: Path) -> dict[str, str]:
    mapping = {}
    current = None
    with path.open() as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">Cluster"):
                current = f"cluster_{int(line.split()[1]) + 1}"
                continue
            m = re.search(r">([^.\s]+)\.\.\.", line)
            if m and current:
                mapping[m.group(1)] = current
    return mapping


def load_mirbase_precursors(path: Path) -> tuple[list[Precursor], dict[str, Precursor]]:
    precursors = []
    by_name = {}
    with path.open() as f:
        next(f)
        for line in f:
            if not line.strip():
                continue
            parts = line.split()
            if len(parts) < 6:
                continue
            name, acc, chrom, start, end, strand, *_ = parts
            if not start.isdigit() or not end.isdigit():
                continue
            fam = family_of_name(name)
            p = Precursor(name, acc, fam, clean_chr(chrom), int(start), int(end), strand)
            precursors.append(p)
            by_name[name] = p
    return precursors, by_name


def load_pmiren_precursors(path: Path) -> tuple[list[Precursor], dict[str, Precursor]]:
    precursors = []
    by_name = {}
    with path.open() as f:
        header = f.readline().rstrip("\n").split("\t")
        idx = {name: i for i, name in enumerate(header)}
        for line in f:
            if not line.strip():
                continue
            r = line.rstrip("\n").split("\t")
            fam = family_of_name(r[idx["miRNA_family"]])
            p = Precursor(
                r[idx["miRNA_locus_ID"]],
                r[idx["miRNA_locus_accession"]],
                fam,
                clean_chr(r[idx["Chromosome"]]),
                int(r[idx["Start"]]),
                int(r[idx["End"]]),
                r[idx["Strand"]],
            )
            precursors.append(p)
            by_name[p.name] = p
    return precursors, by_name


def load_zh13_hairpin_positions(path: Path) -> dict[str, tuple[str, int, int, str, str]]:
    """Load database hairpin best-hit coordinates on ZH13.

    Query_ID is either a plain precursor ID (miRBase) or a pipe-delimited
    pmiREN header. The first pipe-delimited field is the precursor locus ID.
    """
    positions = {}
    with path.open() as f:
        header = f.readline().rstrip("\n").split("\t")
        idx = {name: i for i, name in enumerate(header)}
        for line in f:
            if not line.strip():
                continue
            r = line.rstrip("\n").split("\t")
            qparts = r[idx["Query_ID"]].split("|")
            name = qparts[0]
            accession = qparts[1] if len(qparts) > 1 else ""
            strand = "+" if r[idx["S_strand"]].lower() == "plus" else "-"
            positions[name] = (
                clean_chr(r[idx["Chr"]]),
                int(r[idx["Genome_start"]]),
                int(r[idx["Genome_end"]]),
                strand,
                accession,
            )
    return positions


def apply_zh13_hairpin_positions(
    precursors: list[Precursor],
    by_name: dict[str, Precursor],
    positions: dict[str, tuple[str, int, int, str, str]],
) -> int:
    updated = 0
    for name, pos in positions.items():
        chrom, start, end, strand, accession = pos
        p = by_name.get(name)
        if p is None:
            fam = family_of_name(name)
            if not fam:
                continue
            p = Precursor(name, accession, fam, chrom, start, end, strand, True)
            precursors.append(p)
            by_name[name] = p
            updated += 1
            continue
        p.chr, p.start, p.end, p.strand = chrom, start, end, strand
        if accession and not p.accession:
            p.accession = accession
        p.anchorable = True
        updated += 1
    for p in precursors:
        if p.name not in positions:
            p.anchorable = False
    return updated


def precursor_overlap_bp(row: Row, p: Precursor) -> int:
    if row.chr != p.chr or row.strand != p.strand:
        return 0
    return max(0, min(row.h_end, p.end) - max(row.h_start, p.start) + 1)


def load_blast(path: Path, qseqs: dict[str, str], sseqs: dict[str, str]) -> dict[str, list[Hit]]:
    hits = defaultdict(list)
    with path.open() as f:
        for line in f:
            if not line.strip():
                continue
            r = line.rstrip("\n").split("\t")
            if len(r) < 15:
                continue
            q, s = r[0], r[1]
            hit = Hit(
                qseqid=q,
                sseqid=s,
                pident=float(r[2]),
                length=int(r[3]),
                qstart=int(r[6]),
                qend=int(r[7]),
                sstart=int(r[8]),
                send=int(r[9]),
                evalue=r[10],
                bitscore=float(r[11]),
                qseq=norm_seq(r[12]),
                sseq=norm_seq(r[13]),
                stitle=r[14],
                qlen=len(qseqs.get(q, r[12])),
                slen=len(sseqs.get(s, r[13])),
            )
            hits[q].append(hit)
    for q in hits:
        hits[q].sort(key=lambda h: (-h.pident, -h.length, -h.bitscore, h.sseqid))
    return hits


def same_anchor_candidates(row: Row, family: str, precursors: list[Precursor]) -> list[Precursor]:
    cands = []
    for p in precursors:
        if not p.anchorable:
            continue
        if p.family != family:
            continue
        if precursor_overlap_bp(row, p) >= MIN_PRECURSOR_OVERLAP_BP:
            cands.append(p)
    cands.sort(key=lambda p: (-precursor_overlap_bp(row, p), p.start, p.end, p.name))
    return cands


def exact_reported_candidate(row: Row, hit: Hit, by_precursor: dict[str, Precursor]) -> Precursor | None:
    pname = precursor_from_mature(hit.sseqid)
    p = by_precursor.get(pname)
    if not p:
        return None
    if p.anchorable and precursor_overlap_bp(row, p) >= MIN_PRECURSOR_OVERLAP_BP:
        return p
    return None


def status_variant(status: str) -> str:
    if status == "reported":
        return "known_member_variant"
    if status.endswith("_variant"):
        return status
    if status == "unreported_mature_arm":
        return "unreported_mature_arm_variant"
    if status == "known_family_new_member":
        return "known_family_new_member_variant"
    if status == "new_family_new_member":
        return "new_family_new_member_variant"
    if status == "new_family_member":
        return "new_family_member_variant"
    return status


def base_annotation(annotation: str) -> str:
    return re.sub(r"\.v\d+$", "", annotation)


def mirna_locus(annotation: str) -> str:
    """Annotation root without variant or mature-arm suffix."""
    return normalize_mature_annotation_name(re.sub(r"-(?:5p|3p)$", "", base_annotation(annotation), flags=re.I))


def output_status(status: str) -> str:
    if status == "new_family_new_member":
        return "new_family_member"
    if status == "new_family_new_member_variant":
        return "new_family_member_variant"
    return status


def variant_name(base: str, n: int) -> str:
    return f"{base}.v{n}"


def with_arm(base: str, arm: str) -> str:
    return f"{re.sub(r'-(?:5p|3p)$', '', base, flags=re.I)}-{arm}"


def mature_region_similar(query: str, mature: str) -> bool:
    """Detect same mature-arm isomiRs despite 5'/3' length shifts."""
    if not query or not mature:
        return False
    q = norm_seq(query)
    m = norm_seq(mature)
    best = mature_region_overlap_score(q, m)
    return best >= 17 and best / min(len(q), len(m)) >= 0.8


def mature_region_overlap_score(query: str, mature: str) -> int:
    """Longest exact shared segment between two mature sequences."""
    if not query or not mature:
        return 0
    q = norm_seq(query)
    m = norm_seq(mature)
    best = 0
    for qi in range(len(q)):
        for mi in range(len(m)):
            n = 0
            while qi + n < len(q) and mi + n < len(m) and q[qi + n] == m[mi + n]:
                n += 1
            best = max(best, n)
    return best


def assign(row: Row, status: str, annotation: str, source: str, evidence: str, hit: Hit | None = None, p: Precursor | None = None):
    row.status = status
    row.annotation = annotation
    row.source = source
    row.evidence = evidence
    if hit:
        row.matched_mature = hit.sseqid
    if p:
        row.matched_precursor = p.name
        row.matched_chr = p.chr
        row.matched_start = str(p.start)
        row.matched_end = str(p.end)
        row.matched_strand = p.strand
        row.distance = str(precursor_overlap_bp(row, p))


def infer_arm(row: Row, anchor_row: Row) -> str:
    # Only called for the same or near-identical precursor interval.
    if row.strand == "+":
        return "5p" if row.m_start < anchor_row.m_start else "3p"
    return "5p" if row.m_start > anchor_row.m_start else "3p"


def row_sort_key(rows: list[Row], i: int):
    r = rows[i]
    try:
        c = (0, int(r.chr))
    except ValueError:
        c = (1, r.chr)
    return (c, r.strand, r.m_start, r.m_end, r.file_line)


def propagate_clusters(
    rows: list[Row],
    variant_members: dict[str, list[int]],
    just_indices: set[int] | None = None,
    precursors_by_source: dict[str, list[Precursor]] | None = None,
    only_statuses: set[str] | None = None,
):
    for cid, members in variant_members.items():
        annotated = [i for i in members if rows[i].status]
        if not annotated:
            continue
        if just_indices is not None and not any(i in just_indices for i in annotated):
            continue
        # Prefer highest evidence via fixed order, then read count.
        priority = {
            "reported": 0,
            "known_member_variant": 1,
            "unreported_mature_arm": 2,
            "unreported_mature_arm_variant": 3,
            "known_family_new_member": 4,
            "known_family_new_member_variant": 5,
            "new_family_new_member": 6,
            "new_family_new_member_variant": 7,
            "new_family_member": 8,
            "new_family_member_variant": 9,
        }
        rep_i = sorted(
            annotated,
            key=lambda i: (priority.get(rows[i].status, 99), -read_count(rows[i].seq_id), row_sort_key(rows, i)),
        )[0]
        rep = rows[rep_i]
        b = base_annotation(rep.annotation)
        target_status = status_variant(rep.status)
        # Renumber every non-representative member in this cluster by read count. A position
        # variant cluster is a single mature locus, so lower/equal evidence names in the same
        # cluster must be harmonized to the representative base instead of keeping separate names.
        others = [i for i in members if i != rep_i]
        others.sort(key=lambda i: (-read_count(rows[i].seq_id), -len(rows[i].seq), row_sort_key(rows, i)))
        for rank, i in enumerate(others, 1):
            if only_statuses is not None and rows[i].status not in only_statuses:
                continue
            target_precursor = None
            if precursors_by_source is not None and rep.source in precursors_by_source and rep.matched_precursor:
                for p in precursors_by_source[rep.source]:
                    if p.name == rep.matched_precursor:
                        target_precursor = p
                        break
                if target_precursor is not None and precursor_overlap_bp(rows[i], target_precursor) < MIN_PRECURSOR_OVERLAP_BP:
                    continue
            assign(
                rows[i],
                target_status,
                variant_name(b, rank),
                rep.source,
                f"Propagated from {rep.seq_id} through position variant cluster {cid}",
                p=target_precursor,
            )


def annotate_with_hits(
    rows: list[Row],
    hits_by_query: dict[str, list[Hit]],
    precursors: list[Precursor],
    by_precursor: dict[str, Precursor],
    mature_seqs: dict[str, str],
    source: str,
    prefix: str,
    counters: dict[str, int],
    layer: str,
) -> set[int]:
    mature_ids = set(mature_seqs)
    changed = set()
    for row in rows:
        if row.status:
            continue
        hits = hits_by_query.get(row.seq_id, [])
        if source == "miRbase":
            hits = [h for h in hits if h.sseqid.startswith("gma-")]
        else:
            hits = [h for h in hits if h.sseqid.startswith("Gma-")]
        if layer == "full":
            hits = [h for h in hits if h.is_full_exact]
        elif layer == "pident100_nonfull":
            hits = [h for h in hits if h.is_pident100 and not h.is_full_exact]
        elif layer == "non100":
            hits = [h for h in hits if not h.is_pident100]
        if not hits:
            continue
        hit = hits[0]
        fam = family_of_name(hit.sseqid)
        if not fam:
            continue
        if layer == "full":
            reported = []
            for h in hits:
                p_report = exact_reported_candidate(row, h, by_precursor)
                if p_report:
                    reported.append((-precursor_overlap_bp(row, p_report), h.sseqid, h, p_report))
            if reported:
                _overlap, _name, hit, p_report = sorted(reported, key=lambda x: (x[0], x[1]))[0]
                assign(row, "reported", hit.sseqid, source, "full-length 100% mature hit; ZH13 precursor intervals overlap", hit, p_report)
                changed.add(row.idx)
                continue
        cands = same_anchor_candidates(row, fam, precursors)
        if cands:
            p = cands[0]
            ann = mature_from_precursor(p.name)
            arm_match = None
            for h in hits:
                if family_of_name(h.sseqid) != fam:
                    continue
                arm_match = re.search(r"-(5p|3p)$", h.sseqid, flags=re.I)
                if arm_match:
                    break
            target_similar = mature_region_similar(row.seq, mature_seqs.get(ann, ""))
            if arm_match:
                arm_ann = f"{ann}-{arm_match.group(1).lower()}"
                arm_seq = mature_seqs.get(arm_ann, "")
                arm_similar = mature_region_similar(row.seq, arm_seq)
                if arm_ann in mature_ids and arm_similar:
                    ann = arm_ann
                    assign(row, "known_member_variant", variant_name(ann, 1), source, f"{layer} mature evidence; anchored by overlapping same-family precursor on ZH13", hit, p)
                elif target_similar:
                    assign(row, "known_member_variant", variant_name(ann, 1), source, f"{layer} arm-named hit but query is similar to target mature; anchored by overlapping same-family precursor on ZH13", hit, p)
                elif arm_ann not in mature_ids:
                    assign(row, "unreported_mature_arm", arm_ann, source, f"{layer} arm-specific mature evidence; anchored by overlapping same-family precursor on ZH13", hit, p)
                else:
                    ann = next_family_name(prefix, fam, counters)
                    assign(
                        row,
                        "known_family_new_member",
                        ann,
                        source,
                        f"{layer} arm-specific mature evidence but query is dissimilar to the overlapping reference arm; not assigned as a .v variant",
                        hit,
                        p,
                    )
            else:
                if target_similar:
                    assign(row, "known_member_variant", variant_name(ann, 1), source, f"{layer} mature evidence; anchored by overlapping same-family precursor on ZH13", hit, p)
                else:
                    ann = next_family_name(prefix, fam, counters)
                    assign(row, "known_family_new_member", ann, source, f"{layer} mature family evidence but query is dissimilar to overlapping target mature; no explicit same-species arm evidence", hit, p)
        else:
            ann = next_family_name(prefix, fam, counters)
            assign(row, "known_family_new_member", ann, source, f"{layer} mature evidence; no suitable existing same-family anchor", hit, None)
        changed.add(row.idx)
    return changed


def annotate_other_species(
    rows: list[Row],
    hits_by_query: dict[str, list[Hit]],
    precursors: list[Precursor],
    counters: dict[str, int],
    variant_members: dict[str, list[int]],
    source: str,
    prefix: str,
) -> set[int]:
    changed = set()
    # Existing annotated rows by exact precursor interval to infer unreported arms.
    by_precursor_interval = defaultdict(list)
    for r in rows:
        if r.status and r.source == source:
            by_precursor_interval[(r.chr, r.h_start, r.h_end, r.strand)].append(r)

    for row in rows:
        if row.status:
            continue
        hits = [h for h in hits_by_query.get(row.seq_id, []) if not h.sseqid.lower().startswith(f"{prefix.lower()}-")]
        if not hits:
            continue
        hit = hits[0]
        fam = family_of_name(hit.sseqid)
        if not fam:
            continue
        same_prec = [
            r for r in by_precursor_interval.get((row.chr, row.h_start, row.h_end, row.strand), [])
            if family_of_name(r.annotation) == fam and opposite_arm_row_groups_supported([row], [r])
        ]
        if same_prec:
            anchor = sorted(same_prec, key=lambda r: (-read_count(r.seq_id), r.file_line))[0]
            arm = infer_arm(row, anchor)
            ann = with_arm(base_annotation(anchor.annotation), arm)
            assign(row, "unreported_mature_arm", ann, source, "other-species family/arm evidence on the same predicted precursor", hit, None)
        else:
            cands = same_anchor_candidates(row, fam, precursors)
            if cands:
                # Conservative: do not invent 5p/3p across different predicted precursors.
                ann = next_family_name(prefix, fam, counters)
                assign(row, "known_family_new_member", ann, source, "other-species family evidence; no same-predicted-precursor arm", hit, cands[0])
            else:
                ann = next_family_name(prefix, fam, counters)
                if any(p.family == fam for p in precursors):
                    assign(row, "known_family_new_member", ann, source, f"known {prefix} {source} family; no suitable existing member anchor", hit, None)
                else:
                    assign(row, "new_family_new_member", ann, source, f"family absent from {prefix} {source}", hit, None)
        changed.add(row.idx)
    return changed


def init_family_counters(mature_ids: list[str], precursor_names: list[str], prefix: str) -> dict[str, int]:
    counters = defaultdict(int)
    for name in list(mature_ids) + [mature_from_precursor(n) for n in precursor_names]:
        fam = family_of_name(name)
        if not fam:
            continue
        if prefix == "gma" and not name.startswith("gma-"):
            continue
        if prefix == "Gma" and not name.startswith("Gma-"):
            continue
        counters[f"{prefix}:{fam}"] = max(counters[f"{prefix}:{fam}"], suffix_to_int(suffix_of_member(name)))
    return counters


def assign_soymir(rows: list[Row], variant_by_row: dict[int, str], variant_members: dict[str, list[int]], cdhit: dict[str, str]):
    nohit = [i for i, r in enumerate(rows) if not r.status]
    if not nohit:
        return

    # Build loci: position variant cluster first, singleton otherwise.
    loci = defaultdict(list)
    for i in nohit:
        cid = variant_by_row.get(i)
        loci[("variant", cid) if cid else ("singleton", i)].append(i)

    def primary_member(members):
        return sorted(members, key=lambda i: (-read_count(rows[i].seq_id), -len(rows[i].seq), row_sort_key(rows, i)))[0]

    by_cdhit = defaultdict(list)
    for key, members in loci.items():
        p = primary_member(members)
        by_cdhit[cdhit.get(rows[p].seq_id, f"singleton_{rows[p].seq_id}")].append((key, members, p))

    def loc_bounds(members):
        return {
            "chr": rows[members[0]].chr,
            "strand": rows[members[0]].strand,
            "h_start": min(rows[i].h_start for i in members),
            "h_end": max(rows[i].h_end for i in members),
            "m_start": min(rows[i].m_start for i in members),
            "m_end": max(rows[i].m_end for i in members),
        }

    def soymir_arm_pair_units(locs):
        """Pair no-hit loci on the same predicted precursor before assigning soymir names."""
        items = []
        for key, members, p in locs:
            b = loc_bounds(members)
            items.append({"key": key, "members": members, "p": p, **b})
        used = set()
        units = []
        ordered = sorted(range(len(items)), key=lambda n: row_sort_key(rows, items[n]["p"]))
        for i in ordered:
            if i in used:
                continue
            left = items[i]
            candidates = []
            for j in ordered:
                if j == i or j in used:
                    continue
                right = items[j]
                if left["chr"] != right["chr"] or left["strand"] != right["strand"]:
                    continue
                if not opposite_arm_loci_supported(left, right):
                    continue
                h_ov = interval_overlap_ratio(left["h_start"], left["h_end"], right["h_start"], right["h_end"])
                candidates.append((j, h_ov, abs(((left["h_start"] + left["h_end"]) / 2) - ((right["h_start"] + right["h_end"]) / 2))))
            if candidates:
                j = sorted(candidates, key=lambda x: (-x[1], x[2], row_sort_key(rows, items[x[0]]["p"])))[0][0]
                units.append([left, items[j]])
                used.update({i, j})
            else:
                units.append([left])
                used.add(i)
        return units

    cluster_order = sorted(
        by_cdhit,
        key=lambda c: (
            -sum(len(members) for _key, members, _p in by_cdhit[c]),
            min(row_sort_key(rows, loc[2]) for loc in by_cdhit[c]),
            c,
        ),
    )
    base_names = {cid: f"gma_miRNC{n:04d}" for n, cid in enumerate(cluster_order, 1)}

    for cid in cluster_order:
        units = sorted(soymir_arm_pair_units(by_cdhit[cid]), key=lambda unit: min(row_sort_key(rows, item["p"]) for item in unit))
        use_letters = len(units) > 1
        for loc_idx, unit in enumerate(units):
            suffix = int_to_suffix(loc_idx + 1) if use_letters else ""
            locus_name = base_names[cid] + suffix
            if len(unit) == 2:
                fivep_item = min(unit, key=lambda x: x["m_start"]) if unit[0]["strand"] == "+" else max(unit, key=lambda x: x["m_start"])
                anchor_item = sorted(unit, key=lambda x: (-read_count(rows[x["p"]].seq_id), -len(rows[x["p"]].seq), row_sort_key(rows, x["p"])))[0]
                for item in unit:
                    arm = "5p" if item is fivep_item else "3p"
                    arm_name = with_arm(locus_name, arm)
                    ordered = sorted(item["members"], key=lambda i: (-read_count(rows[i].seq_id), -len(rows[i].seq), row_sort_key(rows, i)))
                    for rank, i in enumerate(ordered):
                        if item is anchor_item:
                            if rank == 0:
                                assign(rows[i], "new_family_member", locus_name, "soymir", "no miRBase/pmiREN mature evidence; CD-HIT-assisted same-precursor candidate locus")
                            else:
                                assign(rows[i], "new_family_member_variant", variant_name(locus_name, rank), "soymir", "position variant of a no-hit soymir candidate")
                        else:
                            if rank == 0:
                                assign(rows[i], "unreported_mature_arm", arm_name, "soymir", "no miRBase/pmiREN mature evidence; paired as the other arm of a soymir no-hit locus")
                            else:
                                assign(rows[i], "unreported_mature_arm_variant", variant_name(arm_name, rank), "soymir", "position variant of a no-hit soymir arm")
            else:
                members = unit[0]["members"]
                ordered = sorted(members, key=lambda i: (-read_count(rows[i].seq_id), -len(rows[i].seq), row_sort_key(rows, i)))
                for rank, i in enumerate(ordered):
                    if rank == 0:
                        assign(rows[i], "new_family_member", locus_name, "soymir", "no miRBase/pmiREN mature evidence; CD-HIT-assisted candidate family")
                    else:
                        assign(rows[i], "new_family_member_variant", variant_name(locus_name, rank), "soymir", "position variant of a no-hit soymir candidate")


def renumber_soymir_families_by_size(rows: list[Row]):
    groups = defaultdict(list)
    for r in rows:
        if r.source != "soymir" or not r.annotation:
            continue
        fam = family_label(r.annotation)
        if fam.startswith("gma_miRNC"):
            groups[fam].append(r)
    if not groups:
        return

    def natural_mirnc_key(fam: str):
        m = re.search(r"miRNC(\d+)", fam)
        return int(m.group(1)) if m else 10**9

    ordered = sorted(
        groups.items(),
        key=lambda item: (
            -len(item[1]),
            min(row_sort_key(rows, r.idx) for r in item[1]),
            natural_mirnc_key(item[0]),
            item[0],
        ),
    )
    fam_map = {old: f"gma_miRNC{n:04d}" for n, (old, _members) in enumerate(ordered, 1)}

    for old, new in fam_map.items():
        if old == new:
            continue
        for r in groups[old]:
            r.annotation = r.annotation.replace(old, new, 1)


def renumber_variants_by_base(rows: list[Row]):
    groups = defaultdict(list)
    for r in rows:
        if re.search(r"\.v\d+$", r.annotation):
            groups[base_annotation(r.annotation)].append(r)
    for base, members in groups.items():
        members.sort(key=lambda r: (-read_count(r.seq_id), -len(r.seq), r.file_line))
        for n, r in enumerate(members, 1):
            r.annotation = variant_name(base, n)


def split_known_member_variants_across_loci(rows: list[Row], variant_by_row: dict[int, str], counters: dict[str, int]):
    """Prevent one known member name from spanning non-overlapping mature loci."""
    by_base = defaultdict(list)
    for r in rows:
        if r.status == "known_member_variant" and re.search(r"\.v\d+$", r.annotation):
            by_base[base_annotation(r.annotation)].append(r)
        elif r.status == "reported" and r.annotation:
            by_base[base_annotation(r.annotation)].append(r)

    for base, members in by_base.items():
        loci = defaultdict(list)
        for r in members:
            loci[variant_by_row.get(r.idx, f"singleton:{r.idx}")].append(r)
        if len(loci) <= 1:
            continue

        # If the base annotation itself exists (e.g. a reported row), keep variants
        # in that same position locus with the base instead of splitting them away.
        anchor_loci = {
            variant_by_row.get(r.idx, f"singleton:{r.idx}")
            for r in rows
            if r.annotation == base
        }
        if anchor_loci:
            ordered_loci = sorted(
                loci.items(),
                key=lambda kv: (0 if kv[0] in anchor_loci else 1, -len(kv[1]), -max(read_count(r.seq_id) for r in kv[1]), min(r.file_line for r in kv[1])),
            )
        else:
            ordered_loci = sorted(
                loci.items(),
                key=lambda kv: (-len(kv[1]), -max(read_count(r.seq_id) for r in kv[1]), min(r.file_line for r in kv[1])),
            )

        # Keep the anchor/largest/highest-read locus as the actual known-member variant.
        keep_key, keep_members = ordered_loci[0]
        keep_members.sort(key=lambda r: (-read_count(r.seq_id), -len(r.seq), r.file_line))
        variant_rank = 1
        for r in keep_members:
            if r.status == "reported" and r.annotation == base:
                r.evidence += "; kept as the representative reported locus for this known member"
                continue
            r.status = "known_member_variant"
            r.annotation = variant_name(base, variant_rank)
            variant_rank += 1
            r.evidence += "; kept as the representative non-overlapping locus for this known member"

        fam = family_of_name(base)
        prefix = "Gma" if base.startswith("Gma-") else "gma"
        for _, split_members in ordered_loci[1:]:
            anchor_rows = [r for r in rows if r.annotation == base]
            arm_anchor = None
            if anchor_rows and not re.search(r"-(?:5p|3p)$", base, flags=re.I):
                split_h_start = min(r.h_start for r in split_members)
                split_h_end = max(r.h_end for r in split_members)
                split_m_start = min(r.m_start for r in split_members)
                split_m_end = max(r.m_end for r in split_members)
                for a in anchor_rows:
                    split_bounds = {
                        "chr": split_members[0].chr,
                        "strand": split_members[0].strand,
                        "h_start": split_h_start,
                        "h_end": split_h_end,
                        "m_start": split_m_start,
                        "m_end": split_m_end,
                    }
                    if not opposite_arm_loci_supported(row_group_bounds([a]), split_bounds):
                        continue
                    arm_anchor = a
                    break
            if arm_anchor:
                split_members.sort(key=lambda r: (-read_count(r.seq_id), -len(r.seq), r.file_line))
                fivep_is_split = split_members[0].m_start < arm_anchor.m_start if split_members[0].strand == "+" else split_members[0].m_start > arm_anchor.m_start
                arm_base = with_arm(base, "5p" if fivep_is_split else "3p")
                for n, r in enumerate(split_members):
                    if n == 0:
                        r.status = "unreported_mature_arm"
                        r.annotation = arm_base
                    else:
                        r.status = "unreported_mature_arm_variant"
                        r.annotation = variant_name(arm_base, n)
                    r.evidence += f"; retained as unreported mature arm of {base} on the same predicted precursor"
                continue
            new_base = next_family_name(prefix, fam, counters)
            split_members.sort(key=lambda r: (-read_count(r.seq_id), -len(r.seq), r.file_line))
            for n, r in enumerate(split_members):
                if n == 0:
                    r.status = "known_family_new_member"
                    r.annotation = new_base
                else:
                    r.status = "known_family_new_member_variant"
                    r.annotation = variant_name(new_base, n)
                r.evidence += f"; split from {base} because this mature-position locus does not overlap the retained known-member variant locus"


def split_unreported_arms_across_loci(rows: list[Row], variant_by_row: dict[int, str], counters: dict[str, int]):
    """An inferred unreported arm name cannot represent multiple independent precursors."""
    by_ann = defaultdict(list)
    for r in rows:
        if r.status in {"unreported_mature_arm", "unreported_mature_arm_variant"}:
            by_ann[base_annotation(r.annotation)].append(r)

    for ann, members in by_ann.items():
        loci = defaultdict(list)
        for r in members:
            vc = variant_by_row.get(r.idx)
            locus_key = ("variant", vc) if vc else ("precursor", r.chr, r.h_start, r.h_end, r.strand, r.idx)
            loci[locus_key].append(r)
        if len(loci) <= 1:
            continue

        # Keep the locus that is closest to an already reported/known anchor in evidence order;
        # if tied, use highest reads. Split the other independent precursors to family-new names.
        ordered_loci = sorted(
            loci.items(),
            key=lambda kv: (-max(read_count(r.seq_id) for r in kv[1]), min(r.file_line for r in kv[1])),
        )
        keep_key, keep_members = ordered_loci[0]
        keep_members.sort(key=lambda r: (-read_count(r.seq_id), -len(r.seq), r.file_line))
        for n, r in enumerate(keep_members):
            r.status = "unreported_mature_arm" if n == 0 else "unreported_mature_arm_variant"
            r.annotation = ann if n == 0 else variant_name(ann, n)
            r.evidence += "; retained as the representative locus for this inferred unreported arm"

        fam = family_of_name(ann)
        prefix = "Gma" if ann.startswith("Gma-") else "gma"
        for _, split_members in ordered_loci[1:]:
            new_base = next_family_name(prefix, fam, counters)
            split_members.sort(key=lambda r: (-read_count(r.seq_id), -len(r.seq), r.file_line))
            for n, r in enumerate(split_members):
                if n == 0:
                    r.status = "known_family_new_member"
                    r.annotation = new_base
                else:
                    r.status = "known_family_new_member_variant"
                    r.annotation = variant_name(new_base, n)
                r.evidence += f"; split from inferred arm {ann} because that arm name is already used by an independent precursor"


def split_reported_precursors_across_loci(rows: list[Row], variant_by_row: dict[int, str], counters: dict[str, int]):
    """A single database precursor anchor should not be consumed by distant loci."""
    by_precursor = defaultdict(list)
    for r in rows:
        if r.status == "reported" and r.matched_precursor:
            # A precursor can legitimately contribute two mature arms. Split only
            # when the same reported mature annotation is reused by independent loci.
            by_precursor[(r.matched_precursor, r.annotation)].append(r)

    for (precursor, _ann), members in by_precursor.items():
        loci = defaultdict(list)
        for r in members:
            loci[variant_by_row.get(r.idx, f"singleton:{r.idx}")].append(r)
        if len(loci) <= 1:
            continue

        def locus_overlap(rs):
            vals = []
            for r in rs:
                try:
                    vals.append(float(r.distance))
                except ValueError:
                    pass
            return max(vals) if vals else -1

        ordered_loci = sorted(
            loci.items(),
            key=lambda kv: (-locus_overlap(kv[1]), -max(read_count(r.seq_id) for r in kv[1]), min(r.file_line for r in kv[1])),
        )
        keep_key, keep_members = ordered_loci[0]
        for r in keep_members:
            r.evidence += "; retained as the highest-overlap locus for this reported precursor anchor"

        fam = family_of_name(members[0].annotation)
        prefix = "Gma" if members[0].annotation.startswith("Gma-") else "gma"
        for _, split_members in ordered_loci[1:]:
            new_base = next_family_name(prefix, fam, counters)
            split_members.sort(key=lambda r: (-read_count(r.seq_id), -len(r.seq), r.file_line))
            for n, r in enumerate(split_members):
                if n == 0:
                    r.status = "known_family_new_member"
                    r.annotation = new_base
                else:
                    r.status = "known_family_new_member_variant"
                    r.annotation = variant_name(new_base, n)
                r.evidence += f"; split from reported {precursor} because that database precursor is already occupied by a higher-overlap non-overlapping locus"


def rescue_tandem_ordered_known_members(
    rows: list[Row],
    variant_by_row: dict[int, str],
    precursors: list[Precursor],
    mature_seqs: dict[str, str],
    source: str,
    counters: dict[str, int],
):
    """Map local tandem query loci to local tandem database precursors by order.

    This handles conserved mature families where many members have equivalent
    mature BLAST hits and the genome versions are offset as a block.
    """
    db_groups = defaultdict(list)
    for p in precursors:
        if p.family:
            db_groups[(p.family, p.chr, p.strand)].append(p)
    for key in db_groups:
        db_groups[key].sort(key=lambda p: (p.start, p.end, p.name))

    existing_base_loci = defaultdict(set)
    for r in rows:
        if r.source == source and r.annotation:
            existing_base_loci[base_annotation(r.annotation)].add(variant_by_row.get(r.idx, f"singleton:{r.idx}"))

    eligible = {
        "reported",
        "known_member_variant",
        "known_family_new_member",
        "known_family_new_member_variant",
        "unreported_mature_arm",
        "unreported_mature_arm_variant",
    }
    grouped_rows = defaultdict(list)
    for r in rows:
        if r.source != source or r.status not in eligible:
            continue
        fam = family_of_name(r.annotation) or family_of_name(r.matched_mature)
        if not fam:
            continue
        grouped_rows[(fam, r.chr, r.strand)].append(r)

    def precursor_loci(members):
        loci = []
        for r in sorted(members, key=lambda x: (x.h_start, x.h_end, x.file_line)):
            placed = False
            for locus in loci:
                if interval_overlap_ratio(
                    r.h_start,
                    r.h_end,
                    min(x.h_start for x in locus),
                    max(x.h_end for x in locus),
                ) >= 0.8:
                    locus.append(r)
                    placed = True
                    break
            if not placed:
                loci.append([r])
        items = []
        for members in loci:
            rep = sorted(members, key=lambda r: (-read_count(r.seq_id), -len(r.seq), r.file_line))[0]
            locus = "|".join(sorted({variant_by_row.get(r.idx, f"singleton:{r.idx}") for r in members}))
            items.append((locus, members, min(r.h_start for r in members), max(r.h_end for r in members), rep))
        items.sort(key=lambda x: (x[2], x[3], x[4].file_line))
        return items

    def query_clusters(items):
        clusters = []
        current = []
        last_end = None
        for item in items:
            if current and item[2] - last_end > 100_000:
                clusters.append(current)
                current = []
            current.append(item)
            last_end = max(last_end or item[3], item[3])
        if current:
            clusters.append(current)
        return clusters

    def db_clusters(db_items):
        clusters = []
        current = []
        last_end = None
        for p in db_items:
            if current and p.start - last_end > 100_000:
                clusters.append(current)
                current = []
            current.append(p)
            last_end = max(last_end or p.end, p.end)
        if current:
            clusters.append(current)
        return clusters

    def assign_precursor_locus(members: list[Row], precursor: Precursor):
        base = mature_from_precursor(precursor.name)
        mature_seq = mature_seqs.get(base, "")
        ordered = sorted(members, key=lambda r: (-read_count(r.seq_id), -len(r.seq), r.file_line))
        exact = [r for r in ordered if mature_seq and r.seq == mature_seq]
        anchor_rows = exact or [
            r for r in ordered
            if interval_overlap_ratio(r.m_start, r.m_end, ordered[0].m_start, ordered[0].m_end) > 0
        ]
        anchor_start = min(r.m_start for r in anchor_rows)
        anchor_end = max(r.m_end for r in anchor_rows)

        variants = []
        arm_groups = defaultdict(list)
        for r in ordered:
            if interval_overlap_ratio(r.m_start, r.m_end, anchor_start, anchor_end) > 0:
                variants.append(r)
            else:
                fivep_is_row = r.m_start < anchor_start if r.strand == "+" else r.m_start > anchor_start
                arm_groups["5p" if fivep_is_row else "3p"].append(r)

        exact_ids = {r.idx for r in exact}
        variant_base = base
        for arm in ("5p", "3p"):
            arm_base = f"{base}-{arm}"
            if arm_base in mature_seqs and any(base_annotation(r.matched_mature) == arm_base for r in variants):
                variant_base = arm_base
                break
        variant_rank = 1
        dissimilar = []
        for r in variants:
            if r.idx in exact_ids:
                r.status = "reported"
                r.annotation = variant_base
            elif mature_region_similar(r.seq, mature_seqs.get(variant_base, "")):
                r.status = "known_member_variant"
                r.annotation = variant_name(variant_base, variant_rank)
                variant_rank += 1
            else:
                dissimilar.append(r)
                continue
            r.source = source
            r.matched_mature = variant_base if variant_base in mature_seqs else r.matched_mature
            r.matched_precursor = precursor.name
            r.matched_chr = precursor.chr
            r.matched_start = str(precursor.start)
            r.matched_end = str(precursor.end)
            r.matched_strand = precursor.strand
            r.distance = str(precursor_overlap_bp(r, precursor))
            r.evidence += f"; tandem-order rescue mapped this predicted precursor locus to {precursor.name}"

        if dissimilar:
            fam = family_of_name(variant_base)
            prefix = "Gma" if source == "pmiren" else "gma"
            new_base = next_family_name(prefix, fam, counters)
            dissimilar.sort(key=lambda r: (-read_count(r.seq_id), -len(r.seq), r.file_line))
            for n, r in enumerate(dissimilar):
                r.status = "known_family_new_member" if n == 0 else "known_family_new_member_variant"
                r.annotation = new_base if n == 0 else variant_name(new_base, n)
                r.source = source
                r.matched_mature = ""
                r.matched_precursor = precursor.name
                r.matched_chr = precursor.chr
                r.matched_start = str(precursor.start)
                r.matched_end = str(precursor.end)
                r.matched_strand = precursor.strand
                r.distance = str(precursor_overlap_bp(r, precursor))
                r.evidence += (
                    f"; tandem-order rescue found overlapping precursor {precursor.name}, "
                    "but mature sequence is dissimilar to the mapped reference member; kept as same-family new locus"
                )

        for arm, arm_members in arm_groups.items():
            arm_base = with_arm(base, arm)
            for n, r in enumerate(sorted(arm_members, key=lambda x: (-read_count(x.seq_id), -len(x.seq), x.file_line))):
                r.status = "unreported_mature_arm" if n == 0 else "unreported_mature_arm_variant"
                r.annotation = arm_base if n == 0 else variant_name(arm_base, n)
                r.source = source
                r.matched_precursor = precursor.name
                r.matched_chr = precursor.chr
                r.matched_start = str(precursor.start)
                r.matched_end = str(precursor.end)
                r.matched_strand = precursor.strand
                r.distance = str(precursor_overlap_bp(r, precursor))
                r.evidence += f"; tandem-order rescue mapped this predicted precursor locus to {precursor.name} and kept this as unreported {arm}"

    for key, members in grouped_rows.items():
        db_items = db_groups.get(key, [])
        if not db_items:
            continue
        for qcluster in query_clusters(precursor_loci(members)):
            if len(qcluster) < 2:
                continue
            q_start = min(x[2] for x in qcluster)
            q_end = max(x[3] for x in qcluster)
            candidates = []
            for dcluster in db_clusters(db_items):
                if len(dcluster) != len(qcluster):
                    continue
                if any(not p.anchorable for p in dcluster):
                    continue
                pair_overlaps = [
                    max(precursor_overlap_bp(r, p) for r in qitem[1])
                    for qitem, p in zip(qcluster, dcluster)
                ]
                if any(ov < MIN_PRECURSOR_OVERLAP_BP for ov in pair_overlaps):
                    continue
                d_start = min(p.start for p in dcluster)
                d_end = max(p.end for p in dcluster)
                overlap = max(0, min(q_end, d_end) - max(q_start, d_start) + 1)
                if overlap >= MIN_PRECURSOR_OVERLAP_BP:
                    candidates.append((overlap, dcluster))
            if not candidates:
                continue
            _, dcluster = sorted(candidates, key=lambda x: (-x[0], x[1][0].start, x[1][0].name))[0]
            q_loci = set()
            for x in qcluster:
                q_loci.update(x[0].split("|"))
            target_bases = [mature_from_precursor(p.name) for p in dcluster]
            if any(existing_base_loci.get(base, set()) - q_loci for base in target_bases):
                continue
            for (_, members, _qs, _qe, _rep), precursor in zip(qcluster, dcluster):
                assign_precursor_locus(members, precursor)


def interval_overlap_ratio(a_start: int, a_end: int, b_start: int, b_end: int) -> float:
    overlap = max(0, min(a_end, b_end) - max(a_start, b_start) + 1)
    if overlap <= 0:
        return 0.0
    return overlap / min(a_end - a_start + 1, b_end - b_start + 1)


def row_group_bounds(members: list[Row]) -> dict[str, int | str]:
    return {
        "chr": members[0].chr,
        "strand": members[0].strand,
        "h_start": min(r.h_start for r in members),
        "h_end": max(r.h_end for r in members),
        "m_start": min(r.m_start for r in members),
        "m_end": max(r.m_end for r in members),
    }


def mature_interval_inside_precursor(m_start: int, m_end: int, h_start: int, h_end: int, tolerance: int = OPPOSITE_ARM_PRECURSOR_TOLERANCE_NT) -> bool:
    return m_start >= h_start - tolerance and m_end <= h_end + tolerance


def mature_relative_midpoint(item: dict[str, int | str]) -> float:
    h_len = int(item["h_end"]) - int(item["h_start"]) + 1
    if h_len <= 0:
        return 0.5
    return (((int(item["m_start"]) + int(item["m_end"])) / 2) - int(item["h_start"])) / h_len


def opposite_arm_loci_supported(left: dict[str, int | str], right: dict[str, int | str]) -> bool:
    """Return True only for two mature regions that plausibly occupy opposite precursor arms."""
    if left["chr"] != right["chr"] or left["strand"] != right["strand"]:
        return False
    if interval_overlap_ratio(int(left["h_start"]), int(left["h_end"]), int(right["h_start"]), int(right["h_end"])) < OPPOSITE_ARM_MIN_PRECURSOR_OVERLAP:
        return False
    if interval_overlap_ratio(int(left["m_start"]), int(left["m_end"]), int(right["m_start"]), int(right["m_end"])) > 0:
        return False
    if not mature_interval_inside_precursor(int(left["m_start"]), int(left["m_end"]), int(right["h_start"]), int(right["h_end"])):
        return False
    if not mature_interval_inside_precursor(int(right["m_start"]), int(right["m_end"]), int(left["h_start"]), int(left["h_end"])):
        return False
    rel_left = mature_relative_midpoint(left)
    rel_right = mature_relative_midpoint(right)
    return abs(rel_left - rel_right) >= OPPOSITE_ARM_MIN_RELATIVE_SEPARATION


def opposite_arm_row_groups_supported(left_members: list[Row], right_members: list[Row]) -> bool:
    return opposite_arm_loci_supported(row_group_bounds(left_members), row_group_bounds(right_members))


def merge_new_member_precursor_arms(rows: list[Row], variant_by_row: dict[int, str]):
    """Merge two new-family-member mature loci on the same precursor as 5p/3p arms."""
    eligible = {"known_family_new_member", "known_family_new_member_variant"}
    base_loci = defaultdict(set)
    for r in rows:
        if r.annotation:
            base_loci[base_annotation(r.annotation)].add(variant_by_row.get(r.idx, f"singleton:{r.idx}"))

    grouped = defaultdict(lambda: defaultdict(list))
    for r in rows:
        if r.status not in eligible:
            continue
        fam = family_of_name(r.annotation)
        if not fam:
            continue
        locus = variant_by_row.get(r.idx, f"singleton:{r.idx}")
        grouped[(r.source, fam, r.chr, r.strand)][locus].append(r)

    for (_source, _fam, _chr, strand), loci in grouped.items():
        if len(loci) < 2:
            continue
        items = []
        for locus, members in loci.items():
            bases = {base_annotation(r.annotation) for r in members}
            if len(bases) != 1:
                continue
            rep = sorted(members, key=lambda r: (-read_count(r.seq_id), -len(r.seq), r.file_line))[0]
            items.append({
                "locus": locus,
                "members": members,
                "base": next(iter(bases)),
                "chr": rep.chr,
                "strand": rep.strand,
                "h_start": min(r.h_start for r in members),
                "h_end": max(r.h_end for r in members),
                "m_start": min(r.m_start for r in members),
                "m_end": max(r.m_end for r in members),
                "rep": rep,
            })
        for i, left in enumerate(items):
            for right in items[i + 1:]:
                pair_loci = {left["locus"], right["locus"]}
                if base_loci[left["base"]] - pair_loci or base_loci[right["base"]] - pair_loci:
                    continue
                if not opposite_arm_loci_supported(left, right):
                    continue
                rep_item = sorted([left, right], key=lambda x: (-read_count(x["rep"].seq_id), -len(x["rep"].seq), x["rep"].file_line))[0]
                shared_base = mirna_locus(rep_item["base"])
                fivep_item = min([left, right], key=lambda x: x["m_start"]) if strand == "+" else max([left, right], key=lambda x: x["m_start"])
                for item in (left, right):
                    arm = "5p" if item is fivep_item else "3p"
                    is_anchor_item = item is rep_item
                    arm_base = shared_base if is_anchor_item else with_arm(shared_base, arm)
                    ordered = sorted(item["members"], key=lambda r: (-read_count(r.seq_id), -len(r.seq), r.file_line))
                    for n, r in enumerate(ordered):
                        if is_anchor_item:
                            r.status = "known_family_new_member" if n == 0 else "known_family_new_member_variant"
                        else:
                            r.status = "unreported_mature_arm" if n == 0 else "unreported_mature_arm_variant"
                        r.annotation = arm_base if n == 0 else variant_name(arm_base, n)
                        r.evidence += f"; merged with paired mature locus on the same predicted precursor as {arm}"


def merge_new_arm_with_known_precursor(rows: list[Row], variant_by_row: dict[int, str], cdhit: dict[str, str]):
    """Attach a new mature locus to a reported/known locus on the same precursor as an unreported arm."""
    known_status = {"reported", "known_member_variant"}
    new_status = {
        "known_family_new_member",
        "known_family_new_member_variant",
        "new_family_new_member",
        "new_family_new_member_variant",
    }
    grouped = defaultdict(lambda: defaultdict(list))
    for r in rows:
        if not family_of_name(r.annotation):
            continue
        locus = variant_by_row.get(r.idx, f"singleton:{r.idx}")
        grouped[(r.source, r.chr, r.strand)][locus].append(r)

    def cdhit_overlap(left, right) -> bool:
        return bool({cdhit.get(r.seq_id) for r in left["members"]} & {cdhit.get(r.seq_id) for r in right["members"]})

    def has_arm_evidence(item) -> bool:
        return any(re.search(r"-(?:5p|3p)$", r.matched_mature or r.annotation, flags=re.I) for r in item["members"])

    for (_source, _chr, strand), loci in grouped.items():
        items = []
        for locus, members in loci.items():
            rep = sorted(members, key=lambda r: (-read_count(r.seq_id), -len(r.seq), r.file_line))[0]
            statuses = {r.status for r in members}
            if statuses & known_status:
                kind = "known"
            elif statuses <= new_status:
                kind = "new"
            else:
                continue
            items.append({
                "locus": locus,
                "members": members,
                "kind": kind,
                "chr": rep.chr,
                "strand": rep.strand,
                "h_start": min(r.h_start for r in members),
                "h_end": max(r.h_end for r in members),
                "m_start": min(r.m_start for r in members),
                "m_end": max(r.m_end for r in members),
                "rep": rep,
            })
        for known in [x for x in items if x["kind"] == "known"]:
            known_base = base_annotation(known["rep"].annotation)
            if re.search(r"-(?:5p|3p)$", known_base, flags=re.I):
                continue
            for new in [x for x in items if x["kind"] == "new"]:
                same_family = family_of_name(known_base) == family_of_name(new["rep"].annotation)
                if not same_family and not (cdhit_overlap(known, new) and has_arm_evidence(new)):
                    continue
                if not opposite_arm_loci_supported(known, new):
                    continue
                fivep_item = min([known, new], key=lambda x: x["m_start"]) if strand == "+" else max([known, new], key=lambda x: x["m_start"])
                new_arm = "5p" if new is fivep_item else "3p"
                new_base = with_arm(known_base, new_arm)
                ordered = sorted(new["members"], key=lambda r: (-read_count(r.seq_id), -len(r.seq), r.file_line))
                for n, r in enumerate(ordered):
                    r.status = "unreported_mature_arm" if n == 0 else "unreported_mature_arm_variant"
                    r.annotation = new_base if n == 0 else variant_name(new_base, n)
                    r.evidence += f"; merged with {known_base} on the same predicted precursor as unreported {new_arm}"


def rescue_new_family_loci_to_free_precursor_anchors(
    rows: list[Row],
    variant_by_row: dict[int, str],
    precursors: list[Precursor],
    mature_seqs: dict[str, str],
    source: str,
    counters: dict[str, int],
):
    """Map remaining family-new loci to unused same-family database precursor anchors when possible."""
    candidate_status = {"known_family_new_member", "known_family_new_member_variant"}
    occupied_status = {
        "reported",
        "known_member_variant",
        "unreported_mature_arm",
        "unreported_mature_arm_variant",
    }
    db_bases = {mature_from_precursor(p.name) for p in precursors}
    occupied = set()
    for r in rows:
        if r.source != source or r.status not in occupied_status:
            continue
        b = re.sub(r"-(?:5p|3p)$", "", base_annotation(r.annotation), flags=re.I)
        if b in db_bases:
            occupied.add(b)

    # First group position variants into mature loci, then merge mature loci that share a predicted precursor.
    mature_loci = defaultdict(list)
    for r in rows:
        if r.source == source and r.status in candidate_status:
            mature_loci[variant_by_row.get(r.idx, f"singleton:{r.idx}")].append(r)

    items = []
    for locus, members in mature_loci.items():
        rep = sorted(members, key=lambda r: (-read_count(r.seq_id), -len(r.seq), r.file_line))[0]
        items.append({
            "locus": locus,
            "members": members,
            "rep": rep,
            "family": family_label(rep.annotation),
            "chr": rep.chr,
            "strand": rep.strand,
            "h_start": min(r.h_start for r in members),
            "h_end": max(r.h_end for r in members),
            "m_start": min(r.m_start for r in members),
            "m_end": max(r.m_end for r in members),
        })

    groups = []
    for item in sorted(items, key=lambda x: (x["chr"], x["strand"], x["h_start"], x["h_end"], x["locus"])):
        placed = False
        for group in groups:
            if group["family"] != item["family"] or group["chr"] != item["chr"] or group["strand"] != item["strand"]:
                continue
            if interval_overlap_ratio(group["h_start"], group["h_end"], item["h_start"], item["h_end"]) < 0.8:
                continue
            group["items"].append(item)
            group["h_start"] = min(group["h_start"], item["h_start"])
            group["h_end"] = max(group["h_end"], item["h_end"])
            placed = True
            break
        if not placed:
            groups.append({
                "family": item["family"],
                "chr": item["chr"],
                "strand": item["strand"],
                "h_start": item["h_start"],
                "h_end": item["h_end"],
                "items": [item],
            })

    def seq_identity(a: str, b: str) -> tuple[int, int]:
        matches = mature_region_overlap_score(a, b)
        return matches, -abs(len(a) - len(b))

    for group in sorted(groups, key=lambda g: (g["chr"], g["strand"], g["h_start"], g["h_end"])):
        cands = []
        for p in precursors:
            if not p.anchorable:
                continue
            base = mature_from_precursor(p.name)
            if base in occupied:
                continue
            if family_label(base) != group["family"]:
                continue
            if p.chr != group["chr"] or p.strand != group["strand"]:
                continue
            overlap = max(0, min(group["h_end"], p.end) - max(group["h_start"], p.start) + 1)
            if overlap >= MIN_PRECURSOR_OVERLAP_BP:
                cands.append((overlap, p))
        if not cands:
            continue
        _overlap, anchor = sorted(cands, key=lambda x: (-x[0], x[1].start, x[1].name))[0]
        anchor_base = mature_from_precursor(anchor.name)
        occupied.add(anchor_base)

        # Pick the mature locus most similar to the database mature sequence as the known-member side.
        known_item = None
        known_ann = anchor_base
        possible_matures = [anchor_base, f"{anchor_base}-5p", f"{anchor_base}-3p"]
        scored = []
        for item in group["items"]:
            for mature_id in possible_matures:
                seq = mature_seqs.get(mature_id)
                if not seq:
                    continue
                scored.append((seq_identity(item["rep"].seq, seq), item, mature_id))
        if scored:
            _score, known_item, known_ann = sorted(
                scored,
                key=lambda x: (x[0][0], x[0][1], read_count(x[1]["rep"].seq_id), len(x[1]["rep"].seq)),
                reverse=True,
            )[0]
        else:
            known_item = sorted(group["items"], key=lambda x: (-read_count(x["rep"].seq_id), -len(x["rep"].seq), x["rep"].file_line))[0]

        for item in group["items"]:
            ordered = sorted(item["members"], key=lambda r: (-read_count(r.seq_id), -len(r.seq), r.file_line))
            if item is known_item:
                known_seq = mature_seqs.get(known_ann, "")
                exact_ids = {r.idx for r in ordered if known_seq and r.seq == known_seq}
                variant_rank = 1
                dissimilar = []
                for n, r in enumerate(ordered, 1):
                    if r.idx in exact_ids:
                        r.status = "reported"
                        r.annotation = known_ann
                    elif mature_region_similar(r.seq, known_seq):
                        r.status = "known_member_variant"
                        r.annotation = variant_name(known_ann, variant_rank)
                        variant_rank += 1
                    else:
                        dissimilar.append(r)
                        continue
                    r.matched_mature = known_ann if known_ann in mature_seqs else r.matched_mature
                    r.matched_precursor = anchor.name
                    r.matched_chr = anchor.chr
                    r.matched_start = str(anchor.start)
                    r.matched_end = str(anchor.end)
                    r.matched_strand = anchor.strand
                    r.distance = str(precursor_overlap_bp(r, anchor))
                    r.evidence += f"; rescued to unused overlapping same-family precursor anchor {anchor.name}"
                if dissimilar:
                    fam = family_of_name(anchor_base)
                    prefix = "Gma" if source == "pmiren" else "gma"
                    new_base = next_family_name(prefix, fam, counters)
                    dissimilar.sort(key=lambda r: (-read_count(r.seq_id), -len(r.seq), r.file_line))
                    for n, r in enumerate(dissimilar):
                        r.status = "known_family_new_member" if n == 0 else "known_family_new_member_variant"
                        r.annotation = new_base if n == 0 else variant_name(new_base, n)
                        r.matched_mature = ""
                        r.matched_precursor = anchor.name
                        r.matched_chr = anchor.chr
                        r.matched_start = str(anchor.start)
                        r.matched_end = str(anchor.end)
                        r.matched_strand = anchor.strand
                        r.distance = str(precursor_overlap_bp(r, anchor))
                        r.evidence += (
                            f"; unused precursor rescue found overlapping anchor {anchor.name}, "
                            "but mature sequence is dissimilar to the mapped reference member; kept as same-family new locus"
                        )
            else:
                fivep_is_item = item["m_start"] < known_item["m_start"] if group["strand"] == "+" else item["m_start"] > known_item["m_start"]
                arm_base = with_arm(anchor_base, "5p" if fivep_is_item else "3p")
                for n, r in enumerate(ordered):
                    r.status = "unreported_mature_arm" if n == 0 else "unreported_mature_arm_variant"
                    r.annotation = arm_base if n == 0 else variant_name(arm_base, n)
                    r.matched_precursor = anchor.name
                    r.matched_chr = anchor.chr
                    r.matched_start = str(anchor.start)
                    r.matched_end = str(anchor.end)
                    r.matched_strand = anchor.strand
                    r.distance = str(precursor_overlap_bp(r, anchor))
                    r.evidence += f"; rescued as unreported arm of unused overlapping same-family precursor anchor {anchor.name}"


def renumber_family_new_members(rows: list[Row], mature_ids: list[str], precursor_names: list[str], prefix: str, source: str):
    """Compact remaining family-new member names after all rescue steps."""
    db_max = defaultdict(int)
    for name in list(mature_ids) + [mature_from_precursor(n) for n in precursor_names]:
        if prefix == "gma" and not name.startswith("gma-"):
            continue
        if prefix == "Gma" and not name.startswith("Gma-"):
            continue
        fam = family_of_name(name)
        if fam:
            db_max[fam] = max(db_max[fam], suffix_to_int(suffix_of_member(name)))

    statuses = {
        "known_family_new_member",
        "known_family_new_member_variant",
        "new_family_new_member",
        "new_family_new_member_variant",
    }
    groups = defaultdict(list)
    for r in rows:
        if r.source != source or r.status not in statuses:
            continue
        base = base_annotation(r.annotation)
        root = re.sub(r"-(?:5p|3p)$", "", base, flags=re.I)
        fam = family_of_name(root)
        if fam and root.startswith(f"{prefix}-"):
            groups[(fam, root)].append(r)

    def chr_key(chrom: str):
        return (0, int(chrom)) if chrom.isdigit() else (1, chrom)

    by_family = defaultdict(list)
    for (fam, root), members in groups.items():
        by_family[fam].append((root, members))

    root_map = {}
    for fam, items in by_family.items():
        items.sort(
            key=lambda item: (
                min(chr_key(r.chr) for r in item[1]),
                min(r.h_start for r in item[1]),
                min(r.h_end for r in item[1]),
                min(r.file_line for r in item[1]),
            )
        )
        next_idx = db_max[fam]
        for old_root, _members in items:
            next_idx += 1
            m = re.search(r"(miRNC|miRN|miR)(\d+)", fam, flags=re.I)
            root_map[old_root] = f"{prefix}-{m.group(1)}{m.group(2)}{int_to_suffix(next_idx)}"

    for r in rows:
        if r.source != source or r.status not in statuses:
            continue
        base = base_annotation(r.annotation)
        root = re.sub(r"-(?:5p|3p)$", "", base, flags=re.I)
        if root not in root_map:
            continue
        new_base = base.replace(root, root_map[root], 1)
        m = re.search(r"\.v(\d+)$", r.annotation)
        r.annotation = variant_name(new_base, int(m.group(1))) if m else new_base


def resolve_cross_source_name_collisions(rows: list[Row]):
    """Rename generated new members only when final output collides across sources."""
    generated = {
        "known_family_new_member",
        "known_family_new_member_variant",
        "new_family_new_member",
        "new_family_new_member_variant",
        "unreported_mature_arm",
        "unreported_mature_arm_variant",
    }

    def root_name(annotation: str) -> str:
        base = base_annotation(annotation)
        return re.sub(r"-(?:5p|3p)$", "", base, flags=re.I)

    def replace_root(annotation: str, old_root: str, new_root: str) -> str:
        base = base_annotation(annotation)
        new_base = base.replace(old_root, new_root, 1)
        m = re.search(r"\.v(\d+)$", annotation)
        return variant_name(new_base, int(m.group(1))) if m else new_base

    occupied_by_family = defaultdict(set)
    by_norm_root = defaultdict(list)
    for r in rows:
        if not r.annotation:
            continue
        root = root_name(r.annotation)
        fam = family_of_name(root)
        if not fam:
            continue
        occupied_by_family[fam].add(suffix_to_int(suffix_of_member(root)))
        by_norm_root[root.lower()].append(r)

    root_map = {}
    for _norm, members in by_norm_root.items():
        sources = {r.source for r in members}
        if "miRbase" not in sources or "pmiren" not in sources:
            continue
        reported_roots = {root_name(r.annotation) for r in members if r.status == "reported"}
        generated_roots = sorted({
            root_name(r.annotation)
            for r in members
            if r.status in generated and root_name(r.annotation) not in reported_roots
        })
        if not generated_roots:
            continue
        for old_root in generated_roots:
            fam = family_of_name(old_root)
            next_idx = max(occupied_by_family[fam] or {0}) + 1
            occupied_by_family[fam].add(next_idx)
            m = re.search(r"(miRNC|miRN|miR)(\d+)", fam, flags=re.I)
            prefix = old_root.split("-", 1)[0]
            root_map[old_root] = f"{prefix}-{m.group(1)}{m.group(2)}{int_to_suffix(next_idx)}"

    if not root_map:
        return

    for r in rows:
        if not r.annotation or r.status not in generated:
            continue
        root = root_name(r.annotation)
        if root in root_map:
            r.annotation = replace_root(r.annotation, root, root_map[root])


def find_database_precursor_for_annotation(row: Row, precursors: list[Precursor]) -> Precursor | None:
    keys = {row.matched_precursor.lower()} if row.matched_precursor else set()
    if row.annotation:
        keys.add(mirna_locus(row.annotation).lower())
    for p in precursors:
        if p.name.lower() in keys or mature_from_precursor(p.name).lower() in keys:
            return p
    return None


def demote_nonoverlapping_database_member_annotations(
    rows: list[Row],
    variant_by_row: dict[int, str],
    precursors_by_source: dict[str, list[Precursor]],
    counters_by_source: dict[str, dict[str, int]],
):
    """Withdraw member-level database annotations whose named precursor does not overlap."""
    anchored_statuses = {
        "reported",
        "known_member_variant",
    }
    groups = defaultdict(list)
    for r in rows:
        if r.source not in precursors_by_source or r.status not in anchored_statuses:
            continue
        p = find_database_precursor_for_annotation(r, precursors_by_source[r.source])
        if p and precursor_overlap_bp(r, p) < MIN_PRECURSOR_OVERLAP_BP:
            key = (
                r.source,
                family_of_name(r.annotation),
                variant_by_row.get(r.idx, f"singleton:{r.idx}"),
            )
            groups[key].append(r)

    for (source, fam, _locus), members in groups.items():
        if not fam:
            continue
        prefix = "Gma" if source == "pmiren" else "gma"
        base = next_family_name(prefix, fam, counters_by_source[source])
        members.sort(key=lambda r: (-read_count(r.seq_id), -len(r.seq), r.file_line))
        for n, r in enumerate(members):
            r.status = "known_family_new_member" if n == 0 else "known_family_new_member_variant"
            r.annotation = base if n == 0 else variant_name(base, n)
            r.matched_mature = ""
            r.matched_precursor = ""
            r.matched_chr = ""
            r.matched_start = ""
            r.matched_end = ""
            r.matched_strand = ""
            r.distance = ""
            r.evidence += "; withdrawn from member-level database annotation because the named database precursor does not overlap this ZH13 predicted precursor"


def fill_overlapping_database_anchor_fields(rows: list[Row], precursors_by_source: dict[str, list[Precursor]]):
    """Keep matched precursor fields for final database-backed members and arms."""
    anchored_statuses = {
        "reported",
        "known_member_variant",
        "unreported_mature_arm",
        "unreported_mature_arm_variant",
    }
    for r in rows:
        if r.source not in precursors_by_source or r.status not in anchored_statuses:
            continue
        p = find_database_precursor_for_annotation(r, precursors_by_source[r.source])
        if not p:
            continue
        overlap = precursor_overlap_bp(r, p)
        if overlap < MIN_PRECURSOR_OVERLAP_BP:
            continue
        r.matched_precursor = p.name
        r.matched_chr = p.chr
        r.matched_start = str(p.start)
        r.matched_end = str(p.end)
        r.matched_strand = p.strand
        r.distance = str(overlap)


def assert_database_precursor_overlaps(rows: list[Row], precursors_by_source: dict[str, list[Precursor]]):
    """Database-anchored annotations must keep a real ZH13 precursor overlap."""
    anchored_statuses = {
        "reported",
        "known_member_variant",
    }
    bad = []
    for r in rows:
        if r.source not in {"miRbase", "pmiren"} or r.status not in anchored_statuses:
            continue
        p = find_database_precursor_for_annotation(r, precursors_by_source.get(r.source, []))
        if not p:
            continue
        if precursor_overlap_bp(r, p) < MIN_PRECURSOR_OVERLAP_BP:
            bad.append(r)
    if bad:
        examples = "; ".join(
            f"{r.file_line}:{r.seq_id}->{r.annotation} overlap=0"
            for r in bad[:10]
        )
        raise RuntimeError(
            f"{len(bad)} database-anchored annotations have no ZH13 precursor overlap. "
            f"Examples: {examples}"
        )


def harmonize_unreported_arms_to_final_precursor_anchor(rows: list[Row], variant_by_row: dict[int, str]):
    """Make inferred arms follow the final annotation of their same-precursor anchor."""
    arm_status = {"unreported_mature_arm", "unreported_mature_arm_variant"}
    anchor_status = {
        "reported",
        "known_member_variant",
        "known_family_new_member",
        "known_family_new_member_variant",
        "new_family_new_member",
        "new_family_new_member_variant",
    }

    arm_loci = defaultdict(list)
    anchors = []
    for r in rows:
        if r.status in arm_status:
            arm_loci[variant_by_row.get(r.idx, f"singleton:{r.idx}")].append(r)
        elif r.status in anchor_status and r.annotation and not re.search(r"-(?:5p|3p)(?:\.v\d+)?$", r.annotation, flags=re.I):
            anchors.append(r)

    priority = {
        "reported": 0,
        "known_member_variant": 1,
        "known_family_new_member": 2,
        "known_family_new_member_variant": 3,
        "new_family_new_member": 4,
        "new_family_new_member_variant": 5,
    }

    for _locus, members in arm_loci.items():
        fam = family_of_name(members[0].annotation)
        if not fam:
            continue
        h_start = min(r.h_start for r in members)
        h_end = max(r.h_end for r in members)
        m_start = min(r.m_start for r in members)
        m_end = max(r.m_end for r in members)
        candidates = []
        for a in anchors:
            if a.source != members[0].source or family_of_name(a.annotation) != fam:
                continue
            member_bounds = {
                "chr": members[0].chr,
                "strand": members[0].strand,
                "h_start": h_start,
                "h_end": h_end,
                "m_start": m_start,
                "m_end": m_end,
            }
            if not opposite_arm_loci_supported(row_group_bounds([a]), member_bounds):
                continue
            candidates.append(a)
        if not candidates:
            continue
        anchor = sorted(
            candidates,
            key=lambda a: (priority.get(a.status, 99), -read_count(a.seq_id), -len(a.seq), a.file_line),
        )[0]
        anchor_base = base_annotation(anchor.annotation)
        fivep_is_arm = m_start < anchor.m_start if members[0].strand == "+" else m_start > anchor.m_start
        arm = "5p" if fivep_is_arm else "3p"
        arm_base = with_arm(anchor_base, arm)
        ordered = sorted(members, key=lambda r: (-read_count(r.seq_id), -len(r.seq), r.file_line))
        for n, r in enumerate(ordered):
            r.status = "unreported_mature_arm" if n == 0 else "unreported_mature_arm_variant"
            r.annotation = arm_base if n == 0 else variant_name(arm_base, n)
            r.evidence += f"; synchronized inferred arm to final same-precursor anchor {anchor_base}"


def demote_orphan_unreported_arms(
    rows: list[Row],
    variant_by_row: dict[int, str],
    counters: dict[str, int],
    database_bases: set[str],
):
    """Demote only orphan arms whose root is not an existing database member."""
    arm_status = {"unreported_mature_arm", "unreported_mature_arm_variant"}
    non_arm_status = {
        "reported",
        "known_member_variant",
        "known_family_new_member",
        "known_family_new_member_variant",
        "new_family_new_member",
        "new_family_new_member_variant",
        "new_family_member",
        "new_family_member_variant",
    }
    by_locus = defaultdict(list)
    for r in rows:
        if r.status in arm_status:
            by_locus[mirna_locus(r.annotation)].append(r)

    orphan_loci = {}
    orphan_ids = set()
    for locus, members in by_locus.items():
        if locus in database_bases:
            continue
        has_anchor = False
        for a in rows:
            if a.status not in non_arm_status or not a.annotation:
                continue
            if mirna_locus(a.annotation) != locus:
                continue
            for r in members:
                if a.chr != r.chr or a.strand != r.strand:
                    continue
                if interval_overlap_ratio(a.h_start, a.h_end, r.h_start, r.h_end) >= 0.8:
                    has_anchor = True
                    break
            if has_anchor:
                break
        if has_anchor:
            continue
        orphan_loci[locus] = members
        orphan_ids.update(r.idx for r in members)

    occupied_by_family = defaultdict(set)
    for r in rows:
        if r.idx in orphan_ids or not r.annotation:
            continue
        root = mirna_locus(r.annotation)
        fam = family_of_name(root)
        if fam:
            occupied_by_family[fam].add(suffix_to_int(suffix_of_member(root)))

    for locus, members in orphan_loci.items():
        fam = family_of_name(locus)
        if not fam:
            continue
        prefix = "Gma" if locus.startswith("Gma-") else "gma"
        next_idx = max(occupied_by_family[fam] or {0}) + 1
        occupied_by_family[fam].add(next_idx)
        m = re.search(r"(miRNC|miRN|miR)(\d+)", fam, flags=re.I)
        if not m:
            continue
        new_base = f"{prefix}-{m.group(1)}{m.group(2)}{int_to_suffix(next_idx)}"
        mature_loci = defaultdict(list)
        for r in members:
            mature_loci[variant_by_row.get(r.idx, f"singleton:{r.idx}")].append(r)
        loci = list(mature_loci.values())
        loci.sort(key=lambda rs: (min(r.m_start for r in rs), min(r.file_line for r in rs)))
        use_arms = len(loci) > 1
        main_locus = None
        if use_arms:
            main_locus = sorted(
                loci,
                key=lambda rs: (
                    -max(read_count(r.seq_id) for r in rs),
                    -max(len(r.seq) for r in rs),
                    min(r.file_line for r in rs),
                ),
            )[0]
        for loc_members in loci:
            loc_members.sort(key=lambda r: (-read_count(r.seq_id), -len(r.seq), r.file_line))
            if use_arms and loc_members is not main_locus:
                rep = loc_members[0]
                fivep = rep.m_start == min(min(r.m_start for r in x) for x in loci) if rep.strand == "+" else rep.m_start == max(max(r.m_start for r in x) for x in loci)
                base = with_arm(new_base, "5p" if fivep else "3p")
            else:
                base = new_base
            for n, r in enumerate(loc_members):
                if use_arms and loc_members is not main_locus:
                    r.status = "unreported_mature_arm" if n == 0 else "unreported_mature_arm_variant"
                else:
                    r.status = "known_family_new_member" if n == 0 else "known_family_new_member_variant"
                r.annotation = base if n == 0 else variant_name(base, n)
                r.evidence += "; demoted from orphan unreported arm because no same-precursor non-arm anchor was detected"


def attach_soymir_nohit_arms_to_annotated_precursors(rows: list[Row], variant_by_row: dict[int, str], cdhit: dict[str, str]):
    """Reclassify no-hit soymir loci as unreported arms when same-cluster evidence supports it."""
    anchor_status = {
        "reported",
        "known_member_variant",
        "known_family_new_member",
        "known_family_new_member_variant",
        "unreported_mature_arm",
        "unreported_mature_arm_variant",
    }
    soymir_status = {"new_family_member", "new_family_member_variant"}

    anchors = []
    for r in rows:
        if r.source == "soymir" or r.status not in anchor_status or not r.annotation:
            continue
        anchors.append(r)

    soymir_loci = defaultdict(list)
    for r in rows:
        if r.source == "soymir" and r.status in soymir_status:
            soymir_loci[variant_by_row.get(r.idx, f"singleton:{r.idx}")].append(r)

    priority = {
        "reported": 0,
        "known_member_variant": 1,
        "known_family_new_member": 2,
        "known_family_new_member_variant": 3,
        "unreported_mature_arm": 4,
        "unreported_mature_arm_variant": 5,
    }

    for _locus, members in soymir_loci.items():
        cluster_ids = {cdhit.get(r.seq_id) for r in members if cdhit.get(r.seq_id)}
        if not cluster_ids:
            continue
        h_start = min(r.h_start for r in members)
        h_end = max(r.h_end for r in members)
        m_start = min(r.m_start for r in members)
        m_end = max(r.m_end for r in members)
        candidates = []
        for a in anchors:
            if cdhit.get(a.seq_id) not in cluster_ids:
                continue
            member_bounds = {
                "chr": members[0].chr,
                "strand": members[0].strand,
                "h_start": h_start,
                "h_end": h_end,
                "m_start": m_start,
                "m_end": m_end,
            }
            if not opposite_arm_loci_supported(row_group_bounds([a]), member_bounds):
                continue
            candidates.append(a)
        if not candidates:
            continue
        anchor = sorted(
            candidates,
            key=lambda a: (priority.get(a.status, 99), -read_count(a.seq_id), -len(a.seq), a.file_line),
        )[0]
        anchor_base = base_annotation(anchor.annotation)
        fivep_is_soymir = m_start < anchor.m_start if members[0].strand == "+" else m_start > anchor.m_start
        arm_base = with_arm(anchor_base, "5p" if fivep_is_soymir else "3p")
        ordered = sorted(members, key=lambda r: (-read_count(r.seq_id), -len(r.seq), r.file_line))
        for n, r in enumerate(ordered):
            r.status = "unreported_mature_arm" if n == 0 else "unreported_mature_arm_variant"
            r.annotation = arm_base if n == 0 else variant_name(arm_base, n)
            r.source = anchor.source
            r.evidence += (
                f"; no-hit soymir candidate reclassified as unreported arm of {anchor_base} "
                f"using same CD-HIT cluster and same predicted precursor"
            )


def attach_pmiren_arms_to_mirbase_precursors(rows: list[Row], variant_by_row: dict[int, str], cdhit: dict[str, str]):
    """Use pmiREN-only mature evidence as an unreported arm of a matching miRBase precursor."""
    mirbase_status = {
        "reported",
        "known_member_variant",
        "known_family_new_member",
        "known_family_new_member_variant",
    }
    pmiren_status = {"reported", "known_member_variant"}

    anchors = [
        r for r in rows
        if r.source == "miRbase"
        and r.status in mirbase_status
        and r.annotation
        and cdhit.get(r.seq_id)
    ]

    pmiren_loci = defaultdict(list)
    for r in rows:
        if r.source == "pmiren" and r.status in pmiren_status:
            pmiren_loci[variant_by_row.get(r.idx, f"singleton:{r.idx}")].append(r)

    priority = {
        "reported": 0,
        "known_member_variant": 1,
        "known_family_new_member": 2,
        "known_family_new_member_variant": 3,
    }

    for _locus, members in pmiren_loci.items():
        cluster_ids = {cdhit.get(r.seq_id) for r in members if cdhit.get(r.seq_id)}
        if not cluster_ids:
            continue
        h_start = min(r.h_start for r in members)
        h_end = max(r.h_end for r in members)
        m_start = min(r.m_start for r in members)
        m_end = max(r.m_end for r in members)
        fam = family_of_name(members[0].annotation)
        candidates = []
        for a in anchors:
            if family_of_name(a.annotation) != fam:
                continue
            if cdhit.get(a.seq_id) not in cluster_ids:
                continue
            member_bounds = {
                "chr": members[0].chr,
                "strand": members[0].strand,
                "h_start": h_start,
                "h_end": h_end,
                "m_start": m_start,
                "m_end": m_end,
            }
            if not opposite_arm_loci_supported(row_group_bounds([a]), member_bounds):
                continue
            candidates.append(a)
        if not candidates:
            continue
        anchor = sorted(
            candidates,
            key=lambda a: (priority.get(a.status, 99), -read_count(a.seq_id), -len(a.seq), a.file_line),
        )[0]
        anchor_base = base_annotation(anchor.annotation)
        fivep_is_pmiren = m_start < anchor.m_start if members[0].strand == "+" else m_start > anchor.m_start
        arm_base = with_arm(anchor_base, "5p" if fivep_is_pmiren else "3p")
        ordered = sorted(members, key=lambda r: (-read_count(r.seq_id), -len(r.seq), r.file_line))
        for n, r in enumerate(ordered):
            r.status = "unreported_mature_arm" if n == 0 else "unreported_mature_arm_variant"
            r.annotation = arm_base if n == 0 else variant_name(arm_base, n)
            r.source = anchor.source
            r.evidence += (
                f"; pmiREN mature evidence reclassified as unreported arm of {anchor_base} "
                f"using same CD-HIT cluster and same predicted precursor"
            )


def attach_non_mirbase_arms_to_mirbase_precursors_by_overlap(rows: list[Row], variant_by_row: dict[int, str]):
    """Use precursor overlap to attach pmiREN/soymir loci as unreported arms of miRBase anchors.

    Opposite arms of the same precursor often fall into different CD-HIT clusters, so precursor
    overlap and non-overlapping mature intervals are the primary evidence here.
    """
    anchor_status = {
        "reported",
        "known_member_variant",
        "known_family_new_member",
        "known_family_new_member_variant",
    }
    query_status = {
        "reported",
        "known_member_variant",
        "known_family_new_member",
        "known_family_new_member_variant",
        "new_family_member",
        "new_family_member_variant",
    }

    anchors = [
        r for r in rows
        if r.source == "miRbase"
        and r.status in anchor_status
        and r.annotation
    ]

    query_loci = defaultdict(list)
    for r in rows:
        if r.source in {"pmiren", "soymir"} and r.status in query_status:
            query_loci[variant_by_row.get(r.idx, f"singleton:{r.idx}")].append(r)

    priority = {
        "reported": 0,
        "known_member_variant": 1,
        "known_family_new_member": 2,
        "known_family_new_member_variant": 3,
    }

    for _locus, members in query_loci.items():
        h_start = min(r.h_start for r in members)
        h_end = max(r.h_end for r in members)
        m_start = min(r.m_start for r in members)
        m_end = max(r.m_end for r in members)
        candidates = []
        for a in anchors:
            member_bounds = {
                "chr": members[0].chr,
                "strand": members[0].strand,
                "h_start": h_start,
                "h_end": h_end,
                "m_start": m_start,
                "m_end": m_end,
            }
            if not opposite_arm_loci_supported(row_group_bounds([a]), member_bounds):
                continue
            candidates.append(a)
        if not candidates:
            continue
        anchor = sorted(
            candidates,
            key=lambda a: (
                priority.get(a.status, 99),
                -interval_overlap_ratio(h_start, h_end, a.h_start, a.h_end),
                abs(((h_start + h_end) / 2) - a.precursor_mid),
                -read_count(a.seq_id),
                -len(a.seq),
                a.file_line,
            ),
        )[0]
        anchor_base = base_annotation(anchor.annotation)
        fivep_is_query = m_start < anchor.m_start if members[0].strand == "+" else m_start > anchor.m_start
        arm_base = with_arm(anchor_base, "5p" if fivep_is_query else "3p")
        ordered = sorted(members, key=lambda r: (-read_count(r.seq_id), -len(r.seq), r.file_line))
        for n, r in enumerate(ordered):
            original_source = r.source
            r.status = "unreported_mature_arm" if n == 0 else "unreported_mature_arm_variant"
            r.annotation = arm_base if n == 0 else variant_name(arm_base, n)
            r.source = "miRbase"
            r.evidence += (
                f"; {original_source} locus reclassified as unreported arm of miRBase anchor {anchor_base} "
                f"using high precursor overlap and non-overlapping mature intervals"
            )


def attach_unannotated_arms_to_mirbase_precursors_by_overlap(rows: list[Row], variant_by_row: dict[int, str]):
    """After miRBase annotation, attach no-hit loci as unreported arms of miRBase anchors."""
    anchor_status = {
        "reported",
        "known_member_variant",
        "known_family_new_member",
        "known_family_new_member_variant",
    }
    anchors = [
        r for r in rows
        if r.source == "miRbase"
        and r.status in anchor_status
        and r.annotation
    ]

    query_loci = defaultdict(list)
    for r in rows:
        if not r.status:
            query_loci[variant_by_row.get(r.idx, f"singleton:{r.idx}")].append(r)

    priority = {
        "reported": 0,
        "known_member_variant": 1,
        "known_family_new_member": 2,
        "known_family_new_member_variant": 3,
    }

    for _locus, members in query_loci.items():
        h_start = min(r.h_start for r in members)
        h_end = max(r.h_end for r in members)
        m_start = min(r.m_start for r in members)
        m_end = max(r.m_end for r in members)
        candidates = []
        for a in anchors:
            member_bounds = {
                "chr": members[0].chr,
                "strand": members[0].strand,
                "h_start": h_start,
                "h_end": h_end,
                "m_start": m_start,
                "m_end": m_end,
            }
            if not opposite_arm_loci_supported(row_group_bounds([a]), member_bounds):
                continue
            candidates.append(a)
        if not candidates:
            continue
        anchor = sorted(
            candidates,
            key=lambda a: (
                priority.get(a.status, 99),
                -interval_overlap_ratio(h_start, h_end, a.h_start, a.h_end),
                abs(((h_start + h_end) / 2) - a.precursor_mid),
                -read_count(a.seq_id),
                -len(a.seq),
                a.file_line,
            ),
        )[0]
        anchor_base = base_annotation(anchor.annotation)
        fivep_is_query = m_start < anchor.m_start if members[0].strand == "+" else m_start > anchor.m_start
        arm_base = with_arm(anchor_base, "5p" if fivep_is_query else "3p")
        ordered = sorted(members, key=lambda r: (-read_count(r.seq_id), -len(r.seq), r.file_line))
        for n, r in enumerate(ordered):
            r.status = "unreported_mature_arm" if n == 0 else "unreported_mature_arm_variant"
            r.annotation = arm_base if n == 0 else variant_name(arm_base, n)
            r.source = "miRbase"
            r.evidence = (
                f"no miRBase mature hit; reclassified before pmiREN annotation as unreported arm "
                f"of miRBase anchor {anchor_base} using high precursor overlap and non-overlapping mature intervals"
            )


def attach_unannotated_arms_to_pmiren_precursors_by_overlap(rows: list[Row], variant_by_row: dict[int, str]):
    """Before soymir naming, attach no-hit loci as unreported arms of pmiREN anchors."""
    anchor_status = {
        "reported",
        "known_member_variant",
        "known_family_new_member",
        "known_family_new_member_variant",
    }
    anchors = [
        r for r in rows
        if r.source == "pmiren"
        and r.status in anchor_status
        and r.annotation
    ]

    query_loci = defaultdict(list)
    for r in rows:
        if not r.status:
            query_loci[variant_by_row.get(r.idx, f"singleton:{r.idx}")].append(r)

    priority = {
        "reported": 0,
        "known_member_variant": 1,
        "known_family_new_member": 2,
        "known_family_new_member_variant": 3,
    }

    for _locus, members in query_loci.items():
        h_start = min(r.h_start for r in members)
        h_end = max(r.h_end for r in members)
        m_start = min(r.m_start for r in members)
        m_end = max(r.m_end for r in members)
        candidates = []
        for a in anchors:
            member_bounds = {
                "chr": members[0].chr,
                "strand": members[0].strand,
                "h_start": h_start,
                "h_end": h_end,
                "m_start": m_start,
                "m_end": m_end,
            }
            if not opposite_arm_loci_supported(row_group_bounds([a]), member_bounds):
                continue
            candidates.append(a)
        if not candidates:
            continue
        anchor = sorted(
            candidates,
            key=lambda a: (
                priority.get(a.status, 99),
                -interval_overlap_ratio(h_start, h_end, a.h_start, a.h_end),
                abs(((h_start + h_end) / 2) - a.precursor_mid),
                -read_count(a.seq_id),
                -len(a.seq),
                a.file_line,
            ),
        )[0]
        anchor_base = base_annotation(anchor.annotation)
        fivep_is_query = m_start < anchor.m_start if members[0].strand == "+" else m_start > anchor.m_start
        arm_base = with_arm(anchor_base, "5p" if fivep_is_query else "3p")
        ordered = sorted(members, key=lambda r: (-read_count(r.seq_id), -len(r.seq), r.file_line))
        for n, r in enumerate(ordered):
            r.status = "unreported_mature_arm" if n == 0 else "unreported_mature_arm_variant"
            r.annotation = arm_base if n == 0 else variant_name(arm_base, n)
            r.source = "pmiren"
            r.evidence = (
                f"no miRBase/pmiREN mature hit; reclassified before soymir naming as unreported arm "
                f"of pmiREN anchor {anchor_base} using high precursor overlap and non-overlapping mature intervals"
            )


def renumber_variants_within_position_loci(rows: list[Row], variant_by_row: dict[int, str]):
    groups = defaultdict(list)
    for r in rows:
        if re.search(r"\.v\d+$", r.annotation):
            groups[(base_annotation(r.annotation), variant_by_row.get(r.idx, f"singleton:{r.idx}"))].append(r)
    for (base, _), members in groups.items():
        members.sort(key=lambda r: (-read_count(r.seq_id), -len(r.seq), r.file_line))
        for n, r in enumerate(members, 1):
            r.annotation = variant_name(base, n)


def renumber_duplicate_unreported_arm_loci(rows: list[Row], variant_by_row: dict[int, str], cdhit: dict[str, str]):
    """Ensure one inferred arm locus has one base name and .v variants.

    Do not merge unrelated CD-HIT clusters merely because they sit on the same
    inferred precursor arm. Only records sharing a position-variant cluster or a
    CD-HIT cluster are kept as .v variants. Other duplicate arm names are split
    to independent known-family new-member names.
    """
    arm_status = {"unreported_mature_arm", "unreported_mature_arm_variant"}
    grouped = defaultdict(list)
    for r in rows:
        if r.status not in arm_status or not r.annotation:
            continue
        arm_base = base_annotation(r.annotation)
        if not re.search(r"-(?:5p|3p)$", arm_base, flags=re.I):
            continue
        anchor = r.matched_precursor or mirna_locus(arm_base)
        grouped[(r.source, anchor.lower(), arm_base.lower(), r.chr, r.strand)].append(r)

    occupied = {canonical_name_key(mirna_locus(r.annotation)) for r in rows if r.annotation}

    def linked(a: Row, b: Row) -> bool:
        va = variant_by_row.get(a.idx, "")
        vb = variant_by_row.get(b.idx, "")
        if va and va == vb:
            return True
        ca = cdhit.get(a.seq_id, "")
        cb = cdhit.get(b.seq_id, "")
        return bool(ca and ca == cb)

    for members in grouped.values():
        if len(members) < 2:
            continue
        clusters = []
        for r in sorted(members, key=lambda x: (x.m_start, x.m_end, x.file_line)):
            placed = False
            for cluster in clusters:
                if any(linked(r, x) for x in cluster):
                    cluster.append(r)
                    placed = True
                    break
            if not placed:
                clusters.append([r])

        clusters.sort(key=lambda cluster: (-max(read_count(r.seq_id) for r in cluster), min(r.file_line for r in cluster)))
        for cluster_idx, cluster in enumerate(clusters):
            cluster.sort(key=lambda x: (-read_count(x.seq_id), -len(x.seq), x.file_line))
            if cluster_idx == 0:
                base = base_annotation(cluster[0].annotation)
                for n, r in enumerate(cluster):
                    r.status = "unreported_mature_arm" if n == 0 else "unreported_mature_arm_variant"
                    r.annotation = base if n == 0 else variant_name(base, n)
                    if n > 0:
                        r.evidence += f"; renamed as variant of inferred arm {base} within the same variant/CD-HIT cluster"
            else:
                fam = family_of_name(cluster[0].annotation)
                prefix = "gma" if cluster[0].source == "miRbase" else "Gma"
                new_base = next_available_family_name(prefix, fam, occupied)
                for n, r in enumerate(cluster):
                    r.status = "known_family_new_member" if n == 0 else "known_family_new_member_variant"
                    r.annotation = new_base if n == 0 else variant_name(new_base, n)
                    r.evidence += (
                        "; split from duplicate inferred arm name because this record is not in the same "
                        "position-variant cluster or CD-HIT cluster"
                    )


def next_available_family_name(prefix: str, family: str, occupied: set[str]) -> str:
    m = re.search(r"(miRNC|miRN|miR)(\d+)", family, flags=re.I)
    n = 1
    while True:
        candidate = f"{prefix}-{m.group(1)}{m.group(2)}{int_to_suffix(n)}"
        key = canonical_name_key(candidate)
        if key not in occupied:
            occupied.add(key)
            return candidate
        n += 1


def canonical_name_key(name: str) -> str:
    return name.lower()


def promote_soymir_loci_by_cdhit_known_family(rows: list[Row], cdhit: dict[str, str], known_families_by_source: dict[str, set[str]]):
    """Use unambiguous CD-HIT known-family evidence to rescue no-hit soymir loci."""
    soymir_status = {
        "new_family_member",
        "new_family_member_variant",
        "unreported_mature_arm",
        "unreported_mature_arm_variant",
    }
    known_by_cluster = defaultdict(set)
    for r in rows:
        cid = cdhit.get(r.seq_id)
        if not cid or r.source not in {"miRbase", "pmiren"} or not r.annotation:
            continue
        fam = family_of_name(r.annotation)
        if not fam or fam.startswith("miRNC"):
            continue
        known_by_cluster[cid].add((fam, r.source))

    occupied = {canonical_name_key(mirna_locus(r.annotation)) for r in rows if r.annotation}
    soymir_loci = defaultdict(list)
    for r in rows:
        cid = cdhit.get(r.seq_id)
        if cid and r.source == "soymir" and r.status in soymir_status and r.annotation:
            soymir_loci[(cid, mirna_locus(r.annotation))].append(r)

    for (cid, _old_locus), members in sorted(soymir_loci.items(), key=lambda x: min(r.file_line for r in x[1])):
        known = known_by_cluster.get(cid, set())
        families = {fam for fam, _source in known}
        if len(families) != 1:
            continue
        fam = next(iter(families))
        source = "miRbase" if any(src == "miRbase" for _fam, src in known) else "pmiren"
        prefix = "gma" if source == "miRbase" else "Gma"
        new_base = next_available_family_name(prefix, fam, occupied)
        family_exists_in_soybean = fam in known_families_by_source.get(source, set())

        by_arm = defaultdict(list)
        for r in members:
            m = re.search(r"-(5p|3p)(?:\.v\d+)?$", r.annotation, flags=re.I)
            by_arm[m.group(1).lower() if m else ""].append(r)

        for arm, arm_members in by_arm.items():
            base = with_arm(new_base, arm) if arm else new_base
            ordered = sorted(arm_members, key=lambda r: (-read_count(r.seq_id), -len(r.seq), r.file_line))
            for n, r in enumerate(ordered):
                if r.status.startswith("unreported_mature_arm"):
                    r.status = "unreported_mature_arm" if n == 0 else "unreported_mature_arm_variant"
                else:
                    if family_exists_in_soybean:
                        r.status = "known_family_new_member" if n == 0 else "known_family_new_member_variant"
                    else:
                        r.status = "new_family_member" if n == 0 else "new_family_member_variant"
                r.annotation = base if n == 0 else variant_name(base, n)
                r.source = source
                r.evidence += f"; rescued from soymir by unambiguous CD-HIT known-family evidence for {source} {fam}"


def promote_soymir_loci_by_identical_mature_known_family(rows: list[Row], variant_by_row: dict[int, str], known_families_by_source: dict[str, set[str]]):
    """Use identical mature sequence as stronger family evidence than ambiguous CD-HIT clusters."""
    candidate_status = {"new_family_member", "new_family_member_variant"}
    anchor_status = {
        "reported",
        "known_member_variant",
        "known_family_new_member",
        "known_family_new_member_variant",
        "unreported_mature_arm",
        "unreported_mature_arm_variant",
    }
    by_seq = defaultdict(list)
    for r in rows:
        by_seq[r.seq].append(r)

    occupied = {canonical_name_key(mirna_locus(r.annotation)) for r in rows if r.annotation}
    source_priority = {"miRbase": 0, "pmiren": 1}

    for seq, members in by_seq.items():
        known = []
        for r in members:
            if r.source not in source_priority or r.status not in anchor_status or not r.annotation:
                continue
            fam = family_of_name(r.annotation)
            if fam and not fam.startswith("miRNC"):
                known.append((source_priority[r.source], r.source, fam))
        families = {fam for _prio, _source, fam in known}
        if len(families) != 1:
            continue
        _prio, source, fam = sorted(known)[0]
        family_exists_in_soybean = fam in known_families_by_source.get(source, set())
        if not family_exists_in_soybean:
            continue
        prefix = "gma" if source == "miRbase" else "Gma"

        loci = defaultdict(list)
        for r in members:
            if r.source == "soymir" and r.status in candidate_status:
                loci[variant_by_row.get(r.idx, f"singleton:{r.idx}")].append(r)
        for _locus, locus_members in sorted(loci.items(), key=lambda kv: min(r.file_line for r in kv[1])):
            base = next_available_family_name(prefix, fam, occupied)
            ordered = sorted(locus_members, key=lambda r: (-read_count(r.seq_id), -len(r.seq), r.file_line))
            for n, r in enumerate(ordered):
                r.status = "known_family_new_member" if n == 0 else "known_family_new_member_variant"
                r.annotation = base if n == 0 else variant_name(base, n)
                r.source = source
                r.evidence += f"; rescued from soymir by identical mature sequence already assigned to {source} {fam}"


def main():
    rows = load_rows(INPUT_DIR / "2814_precusor_miRNAs.txt")
    rows_by_key = {r.key: r.idx for r in rows}
    variant_by_row, variant_members = load_variant_clusters(INTERMEDIATE_DIR / "2814_precusor_miRNAs_mature_position_variant_clusters_overlap80.tsv", rows_by_key)
    cdhit = load_cdhit_clstr(INTERMEDIATE_DIR / "1588_mature_miRNAs_cdhit_est_c0.8.clstr")

    query_seqs = read_fasta(INPUT_DIR / "1588_mature_miRNAs.fasta")
    mirbase_seqs = read_fasta(INPUT_DIR / "miRbase-mature.fa")
    pmiren_seqs = read_fasta(INPUT_DIR / "pmiren_all_mature_clean.fa")
    mirbase_prec, mirbase_prec_by_name = load_mirbase_precursors(INPUT_DIR / "miRbase_gma_position.txt")
    pmiren_prec, pmiren_prec_by_name = load_pmiren_precursors(INPUT_DIR / "pmiren_gma_position.txt")
    mirbase_zh13_anchors = apply_zh13_hairpin_positions(
        mirbase_prec,
        mirbase_prec_by_name,
        load_zh13_hairpin_positions(INTERMEDIATE_DIR / "gma-hairpin_ZH13v2_genome_positions_best_hit_unique_interval.tsv"),
    )
    pmiren_zh13_anchors = apply_zh13_hairpin_positions(
        pmiren_prec,
        pmiren_prec_by_name,
        load_zh13_hairpin_positions(INTERMEDIATE_DIR / "pmiren_gmax_hairpin_ZH13v2_genome_positions_best_hit_unique_interval.tsv"),
    )

    mirbase_hits = load_blast(INTERMEDIATE_DIR / "1588_mature_miRNAs_vs_miRbase_mature_blastn_e1e-4.tsv", query_seqs, mirbase_seqs)
    pmiren_hits = load_blast(INTERMEDIATE_DIR / "1588_mature_miRNAs_vs_pmiren_mature_blastn_e1e-4.tsv", query_seqs, pmiren_seqs)
    conserved_queries = {
        q
        for hits_by_query in (mirbase_hits, pmiren_hits)
        for q, hits in hits_by_query.items()
        if any(not h.sseqid.lower().startswith("gma-") for h in hits)
    }
    conserved_species_counts = defaultdict(set)
    for hits_by_query in (mirbase_hits, pmiren_hits):
        for q, hits in hits_by_query.items():
            for h in hits:
                species = h.sseqid.split("-", 1)[0]
                if species.lower() != "gma":
                    conserved_species_counts[q].add(species)

    gma_counters = init_family_counters(list(mirbase_seqs), [p.name for p in mirbase_prec], "gma")
    pmiren_counters = init_family_counters(list(pmiren_seqs), [p.name for p in pmiren_prec], "Gma")

    for layer in ["full", "pident100_nonfull", "non100"]:
        changed = annotate_with_hits(rows, mirbase_hits, mirbase_prec, mirbase_prec_by_name, mirbase_seqs, "miRbase", "gma", gma_counters, layer)
        propagate_clusters(rows, variant_members, changed)

    # Current workflow: skip miRBase other-species mature family evidence.
    changed = []

    split_reported_precursors_across_loci(rows, variant_by_row, gma_counters)
    split_known_member_variants_across_loci(rows, variant_by_row, gma_counters)
    split_unreported_arms_across_loci(rows, variant_by_row, gma_counters)
    attach_unannotated_arms_to_mirbase_precursors_by_overlap(rows, variant_by_row)

    for layer in ["full", "pident100_nonfull", "non100"]:
        changed = annotate_with_hits(rows, pmiren_hits, pmiren_prec, pmiren_prec_by_name, pmiren_seqs, "pmiren", "Gma", pmiren_counters, layer)
        propagate_clusters(rows, variant_members, changed)

    # Current workflow: skip pmiREN other-species mature family evidence.
    changed = []

    split_reported_precursors_across_loci(rows, variant_by_row, gma_counters)
    split_reported_precursors_across_loci(rows, variant_by_row, pmiren_counters)
    split_known_member_variants_across_loci(rows, variant_by_row, gma_counters)
    split_known_member_variants_across_loci(rows, variant_by_row, pmiren_counters)
    split_unreported_arms_across_loci(rows, variant_by_row, gma_counters)
    split_unreported_arms_across_loci(rows, variant_by_row, pmiren_counters)
    attach_non_mirbase_arms_to_mirbase_precursors_by_overlap(rows, variant_by_row)
    attach_unannotated_arms_to_pmiren_precursors_by_overlap(rows, variant_by_row)
    assign_soymir(rows, variant_by_row, variant_members, cdhit)
    propagate_clusters(rows, variant_members, None, {"miRbase": mirbase_prec, "pmiren": pmiren_prec})
    split_reported_precursors_across_loci(rows, variant_by_row, gma_counters)
    split_reported_precursors_across_loci(rows, variant_by_row, pmiren_counters)
    split_known_member_variants_across_loci(rows, variant_by_row, gma_counters)
    split_known_member_variants_across_loci(rows, variant_by_row, pmiren_counters)
    split_unreported_arms_across_loci(rows, variant_by_row, gma_counters)
    split_unreported_arms_across_loci(rows, variant_by_row, pmiren_counters)
    rescue_tandem_ordered_known_members(rows, variant_by_row, mirbase_prec, mirbase_seqs, "miRbase", gma_counters)
    rescue_tandem_ordered_known_members(rows, variant_by_row, pmiren_prec, pmiren_seqs, "pmiren", pmiren_counters)
    merge_new_member_precursor_arms(rows, variant_by_row)
    merge_new_arm_with_known_precursor(rows, variant_by_row, cdhit)
    rescue_new_family_loci_to_free_precursor_anchors(rows, variant_by_row, mirbase_prec, mirbase_seqs, "miRbase", gma_counters)
    rescue_new_family_loci_to_free_precursor_anchors(rows, variant_by_row, pmiren_prec, pmiren_seqs, "pmiren", pmiren_counters)
    demote_nonoverlapping_database_member_annotations(
        rows,
        variant_by_row,
        {"miRbase": mirbase_prec, "pmiren": pmiren_prec},
        {"miRbase": gma_counters, "pmiren": pmiren_counters},
    )
    renumber_family_new_members(rows, list(mirbase_seqs), [p.name for p in mirbase_prec], "gma", "miRbase")
    renumber_family_new_members(rows, list(pmiren_seqs), [p.name for p in pmiren_prec], "Gma", "pmiren")
    attach_pmiren_arms_to_mirbase_precursors(rows, variant_by_row, cdhit)
    attach_soymir_nohit_arms_to_annotated_precursors(rows, variant_by_row, cdhit)
    harmonize_unreported_arms_to_final_precursor_anchor(rows, variant_by_row)
    database_bases = (
        {mirna_locus(x) for x in mirbase_seqs}
        | {mirna_locus(mature_from_precursor(p.name)) for p in mirbase_prec}
        | {mirna_locus(x) for x in pmiren_seqs}
        | {mirna_locus(mature_from_precursor(p.name)) for p in pmiren_prec}
    )
    demote_orphan_unreported_arms(rows, variant_by_row, gma_counters, database_bases)
    known_families_by_source = {
        "miRbase": {p.family for p in mirbase_prec} | {family_of_name(x) for x in mirbase_seqs if x.startswith("gma-")},
        "pmiren": {p.family for p in pmiren_prec} | {family_of_name(x) for x in pmiren_seqs if x.startswith("Gma-")},
    }
    promote_soymir_loci_by_cdhit_known_family(rows, cdhit, known_families_by_source)
    promote_soymir_loci_by_identical_mature_known_family(rows, variant_by_row, known_families_by_source)
    # If an anchor is promoted after initial variant propagation, all members of the same
    # mature-position variant cluster must inherit the promoted anchor's source/name.
    propagate_clusters(
        rows,
        variant_members,
        None,
        {"miRbase": mirbase_prec, "pmiren": pmiren_prec},
        {"new_family_member", "new_family_member_variant", "new_family_new_member", "new_family_new_member_variant"},
    )
    merge_new_member_precursor_arms(rows, variant_by_row)
    merge_new_arm_with_known_precursor(rows, variant_by_row, cdhit)
    renumber_family_new_members(rows, list(mirbase_seqs), [p.name for p in mirbase_prec], "gma", "miRbase")
    renumber_family_new_members(rows, list(pmiren_seqs), [p.name for p in pmiren_prec], "Gma", "pmiren")
    harmonize_unreported_arms_to_final_precursor_anchor(rows, variant_by_row)
    resolve_cross_source_name_collisions(rows)
    renumber_duplicate_unreported_arm_loci(rows, variant_by_row, cdhit)
    renumber_variants_within_position_loci(rows, variant_by_row)
    renumber_soymir_families_by_size(rows)
    demote_nonoverlapping_database_member_annotations(
        rows,
        variant_by_row,
        {"miRbase": mirbase_prec, "pmiren": pmiren_prec},
        {"miRbase": gma_counters, "pmiren": pmiren_counters},
    )
    fill_overlapping_database_anchor_fields(rows, {"miRbase": mirbase_prec, "pmiren": pmiren_prec})
    assert_database_precursor_overlaps(rows, {"miRbase": mirbase_prec, "pmiren": pmiren_prec})

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    full_out = RESULTS_DIR / "2814_precusor_miRNAs_annotated_precursor_overlap_unique_anchor_workflow.tsv"
    short_out = RESULTS_DIR / "2814_precusor_miRNAs_annotated_precursor_overlap_unique_anchor_workflow_short.tsv"
    summary_out = RESULTS_DIR / "2814_precusor_miRNAs_annotated_precursor_overlap_unique_anchor_workflow_summary.txt"

    full_header = [
        "File_Line", "Seq-ID", "Sequences_Mature", "Chr", "M_start", "M_end", "H_start", "H_end", "Strand",
        "Status", "Annotation", "miRNA_Locus", "Reported_Status", "Conservation", "Conserved_Species_Count",
        "Source", "Family", "Position_Variant_Cluster", "CDHIT_Cluster_ID",
        "Evidence", "Matched_mature", "Matched_precursor", "Matched_chr", "Matched_start", "Matched_end",
        "Matched_strand", "Matched_precursor_overlap_bp",
    ]
    with full_out.open("w") as f:
        f.write("\t".join(full_header) + "\n")
        for r in rows:
            status = output_status(r.status)
            annotation = normalize_mature_annotation_name(r.annotation)
            matched_mature = normalize_mature_annotation_name(r.matched_mature)
            evidence = normalize_mature_annotation_text(r.evidence)
            f.write("\t".join([
                str(r.file_line), r.seq_id, r.seq, r.chr, str(r.m_start), str(r.m_end), str(r.h_start), str(r.h_end), r.strand,
                status, annotation, mirna_locus(annotation), "reported" if status == "reported" else "unreported",
                "conserved" if r.seq_id in conserved_queries else "specific", str(len(conserved_species_counts.get(r.seq_id, set()))),
                r.source, output_family_label(annotation),
                variant_by_row.get(r.idx, ""), cdhit.get(r.seq_id, ""),
                evidence, matched_mature, r.matched_precursor, r.matched_chr, r.matched_start, r.matched_end,
                r.matched_strand, r.distance,
            ]) + "\n")

    short_rows = []
    for r in rows:
        status = output_status(r.status)
        annotation = normalize_mature_annotation_name(r.annotation)
        short_rows.append([
            r.seq_id, r.seq, r.chr, str(r.m_start), str(r.m_end), str(r.h_start), str(r.h_end), r.strand,
            status, annotation, mirna_locus(annotation), "reported" if status == "reported" else "unreported",
            "conserved" if r.seq_id in conserved_queries else "specific", str(len(conserved_species_counts.get(r.seq_id, set()))),
            r.source, output_family_label(annotation),
            variant_by_row.get(r.idx, ""), cdhit.get(r.seq_id, ""),
        ])
    short_rows.sort(key=lambda x: (x[15].lower(), x[9].lower(), int(x[2]) if x[2].isdigit() else 999, int(x[5]), int(x[6]), x[0]))
    with short_out.open("w") as f:
        f.write("\t".join([
            "Seq-ID", "Sequences_Mature", "Chr", "M_start", "M_end", "H_start", "H_end", "Strand",
            "Status", "Annotation", "miRNA_Locus", "Reported_Status", "Conservation", "Conserved_Species_Count",
            "Source", "Family", "Position_Variant_Cluster", "CDHIT_Cluster_ID",
        ]) + "\n")
        for r in short_rows:
            f.write("\t".join(r) + "\n")

    status_counts = Counter(output_status(r.status) for r in rows)
    source_counts = Counter(r.source for r in rows)
    reported_counts = Counter("reported" if output_status(r.status) == "reported" else "unreported" for r in rows)
    conservation_counts = Counter("conserved" if r.seq_id in conserved_queries else "specific" for r in rows)
    family_count = len({output_family_label(r.annotation) for r in rows})
    with summary_out.open("w") as f:
        f.write("Precursor-overlap unique-anchor workflow annotation summary\n")
        f.write(f"Input_records: {len(rows)}\n")
        f.write(f"Family_count_case_insensitive_mature_family_normalized: {family_count}\n")
        f.write(f"Position_variant_clusters_loaded: {len(variant_members)}\n")
        f.write(f"Rows_in_position_variant_clusters: {len(variant_by_row)}\n")
        f.write(f"CDHIT_query_ids_loaded: {len(cdhit)}\n")
        f.write(f"miRBase_ZH13_hairpin_anchors_loaded: {mirbase_zh13_anchors}\n")
        f.write(f"pmiREN_ZH13_hairpin_anchors_loaded: {pmiren_zh13_anchors}\n")
        f.write("Database_precursor_anchor_rule: hairpin BLAST E-value 1e-10; select one global best genome hit per database precursor; if multiple best hits tie, prefer the hit matching the original database chromosome and strand; require same chromosome, same strand, and precursor interval overlap on ZH13 for annotation\n")
        f.write(f"Full_output: {full_out.name}\n")
        f.write(f"Short_output: {short_out.name}\n")
        f.write("\nStatus_counts:\n")
        for k in sorted(status_counts):
            f.write(f"{k}: {status_counts[k]}\n")
        f.write("\nSource_counts:\n")
        for k in sorted(source_counts):
            f.write(f"{k}: {source_counts[k]}\n")
        f.write("\nReported_status_counts:\n")
        for k in sorted(reported_counts):
            f.write(f"{k}: {reported_counts[k]}\n")
        f.write("\nConservation_counts:\n")
        for k in sorted(conservation_counts):
            f.write(f"{k}: {conservation_counts[k]}\n")

    print(full_out)
    print(short_out)
    print(summary_out)
    print("Family_count_case_insensitive_mature_family_normalized", family_count)
    print("Status_counts")
    for k in sorted(status_counts):
        print(k, status_counts[k])
    print("Source_counts")
    for k in sorted(source_counts):
        print(k, source_counts[k])
    print("Reported_status_counts")
    for k in sorted(reported_counts):
        print(k, reported_counts[k])
    print("Conservation_counts")
    for k in sorted(conservation_counts):
        print(k, conservation_counts[k])


if __name__ == "__main__":
    main()
