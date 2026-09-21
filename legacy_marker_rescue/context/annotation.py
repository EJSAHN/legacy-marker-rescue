"""GFF3 intervals with exact queries and explicit annotation availability."""
from __future__ import annotations
from bisect import bisect_left, bisect_right
from collections import Counter, defaultdict
import gzip
from pathlib import Path
import re
from urllib.parse import unquote
from .io_utils import sha256, read_tsv, checked_relative

class GeneIndex:
    """Query all intersecting features, including long and nested intervals."""

    def __init__(self, intervals):
        self.intervals = sorted(set(((int(s), int(e), str(i)) for s, e, i in intervals)))
        if any((s < 0 or e <= s for s, e, _ in self.intervals)):
            raise ValueError('Invalid gene interval')
        self.starts = [r[0] for r in self.intervals]
        self.sorted_ends = sorted((r[1] for r in self.intervals))
        self.prefix_max_end = []
        maximum = -1
        for _, e, _ in self.intervals:
            maximum = max(maximum, e)
            self.prefix_max_end.append(maximum)

    def overlapping_ids(self, s: int, e: int):
        right = bisect_left(self.starts, e)
        left = bisect_right(self.prefix_max_end, s, 0, right)
        return {ident for gs, ge, ident in self.intervals[left:right] if ge > s and gs < e}

    def query(self, s: int, e: int, flank: int):
        if s < 0 or e <= s or flank < 0:
            raise ValueError('Invalid query coordinates or flank')
        ids = self.overlapping_ids(s, e)
        if ids:
            nearest = 0
        else:
            distances = []
            left = bisect_right(self.sorted_ends, s)
            if left:
                distances.append(s - self.sorted_ends[left - 1])
            right = bisect_left(self.starts, e)
            if right < len(self.starts):
                distances.append(self.starts[right] - e)
            nearest = min(distances) if distances else None
        return {'overlaps_annotated_gene': int(bool(ids)), 'nearest_gene_distance_bp': nearest, 'gene_count_within_flank': len(self.overlapping_ids(max(0, s - flank), e + flank))}

def choose_gff(project: Path, fasta: Path, acc: str, explicit_map: dict | None=None):
    """Never borrow a file from a different assembly or pick conflicting annotations."""
    if explicit_map and acc in explicit_map:
        chosen = checked_relative(project, explicit_map[acc])
        if not chosen.is_file():
            return (None, 'OVERRIDE_FILE_MISSING', [])
        return (chosen, 'EXPLICIT_MAPPING', [chosen])
    choices = sorted((p for p in fasta.parent.iterdir() if p.is_file() and any((p.name.lower().endswith(ext) for ext in ('.gff', '.gff3', '.gff.gz', '.gff3.gz')))))
    if not choices:
        for base in (project / 'annotations', project / 'data' / 'annotations'):
            if base.is_dir():
                for p in base.rglob('*'):
                    if p.is_file() and any((p.name.lower().endswith(ext) for ext in ('.gff', '.gff3', '.gff.gz', '.gff3.gz'))):
                        if acc in p.parts or re.match(re.escape(acc) + '(?:_|\\.)', p.name):
                            choices.append(p)
    choices = sorted(set(choices))
    if not choices:
        return (None, 'MISSING_GFF', [])
    if len(choices) > 1 and len({sha256(p) for p in choices}) != 1:
        return (None, 'AMBIGUOUS_GFF_FILES', choices)
    return (choices[0], 'ASSOCIATED_GFF', choices)

def parse_gff(path: Path, lengths: dict[str, int], assembly: str, feature_types=('gene', 'pseudogene')):
    features = defaultdict(list)
    regions = {}
    types = Counter()
    errors = []
    warnings = []
    duplicates = 0
    rows = 0
    circular = set()
    build_accessions = set()
    opener = gzip.open if path.name.lower().endswith('.gz') else open
    with opener(path, 'rt', encoding='utf-8-sig', errors='strict') as stream:
        for ln, raw in enumerate(stream, 1):
            line = raw.rstrip('\r\n')
            if not line:
                continue
            if line.startswith('##FASTA'):
                break
            if line.startswith('#!genome-build-accession'):
                build_accessions.update(re.findall('GC[AF]_\\d+\\.\\d+', line))
            if line.startswith('##sequence-region'):
                cols = line.split()
                if len(cols) != 4:
                    errors.append(f'line {ln}: malformed sequence-region')
                    continue
                seqid = unquote(cols[1])
                try:
                    start, end = (int(cols[2]) - 1, int(cols[3]))
                except ValueError:
                    errors.append(f'line {ln}: noninteger sequence-region')
                    continue
                if seqid not in lengths or start < 0 or end <= start or (end > lengths.get(seqid, 0)):
                    errors.append(f'line {ln}: sequence-region incompatible with FASTA: {seqid}')
                    continue
                if seqid in regions and regions[seqid] != (start, end):
                    errors.append(f'line {ln}: conflicting sequence-region: {seqid}')
                regions[seqid] = (start, end)
                continue
            if line.startswith('#'):
                continue
            cols = line.split('\t')
            if len(cols) != 9:
                errors.append(f'line {ln}: expected 9 GFF3 columns')
                continue
            seqid, _, ftype, start_s, end_s, _, _, _, attrs = cols
            seqid = unquote(seqid)
            types[ftype] += 1
            rows += 1
            attr = {}
            for part in attrs.split(';'):
                if '=' in part:
                    key, value = part.split('=', 1)
                    attr[key] = unquote(value)
            if ftype == 'region' and attr.get('Is_circular', '').lower() == 'true':
                circular.add(seqid)
            if ftype not in feature_types:
                continue
            try:
                s, e = (int(start_s) - 1, int(end_s))
            except ValueError:
                errors.append(f'line {ln}: noninteger gene coordinate')
                continue
            if seqid not in lengths:
                errors.append(f'line {ln}: gene seqid missing from FASTA: {seqid}')
                continue
            if s < 0 or e <= s or e > lengths[seqid]:
                errors.append(f'line {ln}: gene outside linear FASTA interval: {seqid}')
                continue
            ident = attr.get('ID') or f'{ftype}:line{ln}'
            record = (s, e, ident)
            features[seqid].append(record)
    if build_accessions and build_accessions != {assembly}:
        errors.append('GFF genome-build accession differs from selected FASTA assembly: ' + ','.join(sorted(build_accessions)))
    n_features = sum((len(v) for v in features.values()))
    if not n_features:
        warnings.append('No gene/pseudogene rows; transcripts and CDS are not substituted for genes.')
    usable = not errors and n_features > 0
    indexes, coverage, exports = ({}, {}, [])
    for seqid, length in lengths.items():
        entries = features.get(seqid, [])
        duplicates += len(entries) - len(set(entries))
        if not usable:
            state = 'invalid_gff' if errors else 'no_gene_feature_type'
        elif entries:
            state = 'gene_features_present'
        elif regions.get(seqid) == (0, length):
            state = 'declared_contig_no_gene_features'
        else:
            state = 'annotation_coverage_unknown'
        coverage[seqid] = state
        if usable and state in ('gene_features_present', 'declared_contig_no_gene_features'):
            indexes[seqid] = GeneIndex(entries)
            for s, e, ident in sorted(set(entries)):
                exports.append({'contig': seqid, 'start_0based': s, 'end_exclusive': e, 'feature_id': ident})
    return (indexes, coverage, exports, {'annotation_status': 'USABLE' if usable else 'INVALID_GFF' if errors else 'NO_GENE_FEATURES', 'gff_rows': rows, 'feature_type_counts': dict(types), 'gene_feature_rows': n_features, 'duplicate_gene_rows_collapsed': duplicates, 'gene_covered_contigs': len(indexes), 'gene_bearing_contigs': sum((bool(features.get(c)) for c in indexes)), 'circular_contigs_reported': sorted(circular), 'gff_build_accessions': sorted(build_accessions), 'errors': errors, 'warnings': warnings, 'scope': 'Overlap with listed gene/pseudogene features on explicitly covered contigs; missing annotation is not gene absence.'})
