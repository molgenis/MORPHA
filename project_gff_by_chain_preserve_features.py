#!/usr/bin/env python3
"""
Project GFF/GFF3 features through a UCSC chain without dropping features that span indel gaps.

Why this exists:
  CrossMap gff treats intervals that cross chain gaps as unmappable and may drop a whole
  exon/CDS. For consensus sequences, that is often wrong: an exon/CDS spanning a deletion
  should shrink, and one spanning an insertion should expand.

Intended use:
  Use this as a replacement for `CrossMap gff` after `bcftools consensus --chain`.
  It projects annotation coordinates from the reference sequence to the consensus sequence
  by applying all chain gap edits in order.

Assumptions and limits:
  - GFF coordinates are 1-based inclusive.
  - Chain coordinates are UCSC 0-based half-open.
  - Only plus-strand chains are supported. This is what bcftools consensus normally emits.
  - Insertions are assigned to a feature if the feature contains the left anchor base.
    Insertions exactly at annotation boundaries are biologically ambiguous; this is the
    useful default for coding insertions represented by left-anchored VCF records.
  - CDS phases are recalculated after projection unless --no-recalc-cds-phase is used.
"""

from __future__ import annotations

import argparse
import sys
from collections import defaultdict
from dataclasses import dataclass
from typing import Iterable, Optional


@dataclass(frozen=True)
class GapEdit:
    start: int      # 1-based coordinate of first affected reference base; for insertion, base after anchor
    ref_len: int    # target/reference gap length
    alt_len: int    # query/consensus gap length

    @property
    def delta(self) -> int:
        return self.alt_len - self.ref_len

    @property
    def end(self) -> int:
        return self.start + self.ref_len - 1

    @property
    def kind(self) -> str:
        if self.ref_len == 0 and self.alt_len > 0:
            return "insertion"
        if self.ref_len > 0 and self.alt_len == 0:
            return "deletion"
        if self.ref_len == self.alt_len:
            return "replacement_same_length"
        return "replacement"


@dataclass(frozen=True)
class ChainProjection:
    t_name: str
    q_name: str
    t_start0: int
    q_start0: int
    initial_offset: int
    edits: list[GapEdit]


def parse_chain(chain_path: str) -> ChainProjection:
    header: Optional[list[str]] = None
    body: list[list[int]] = []

    with open(chain_path, "rt") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith("chain"):
                if header is not None:
                    raise ValueError("Only one chain record is supported per chain file")
                header = line.split()
                continue
            if header is None:
                continue
            parts = [int(x) for x in line.split()]
            if len(parts) not in {1, 3}:
                raise ValueError(f"Malformed chain block line: {line!r}")
            body.append(parts)

    if header is None:
        raise ValueError(f"No chain header found in {chain_path}")
    if len(header) < 13:
        raise ValueError(f"Malformed chain header in {chain_path}: {' '.join(header)}")

    # chain score tName tSize tStrand tStart tEnd qName qSize qStrand qStart qEnd id
    t_name = header[2]
    t_strand = header[4]
    t_start0 = int(header[5])
    q_name = header[7]
    q_strand = header[9]
    q_start0 = int(header[10])

    if t_strand != "+" or q_strand != "+":
        raise ValueError(
            f"Only plus/plus chains are supported; got target={t_strand}, query={q_strand}"
        )

    t_pos0 = t_start0
    q_pos0 = q_start0
    edits: list[GapEdit] = []

    for parts in body:
        size = parts[0]
        t_pos0 += size
        q_pos0 += size

        if len(parts) == 3:
            dt, dq = parts[1], parts[2]
            if dt or dq:
                # Gap starts at the 0-based boundary t_pos0. In 1-based inclusive
                # coordinates, deleted/replaced reference bases start at t_pos0 + 1.
                edits.append(GapEdit(start=t_pos0 + 1, ref_len=dt, alt_len=dq))
            t_pos0 += dt
            q_pos0 += dq

    return ChainProjection(
        t_name=t_name,
        q_name=q_name,
        t_start0=t_start0,
        q_start0=q_start0,
        initial_offset=q_start0 - t_start0,
        edits=edits,
    )


def patch_interval_once(start: int, end: int, edit: GapEdit) -> Optional[tuple[int, int]]:
    """Project one 1-based inclusive interval through one current-coordinate edit."""
    if start > end:
        return None

    delta = edit.delta

    if edit.ref_len == 0 and edit.alt_len == 0:
        return start, end

    # Pure insertion. The insertion is after anchor = edit.start - 1.
    # Assign inserted sequence to any feature containing that anchor base.
    if edit.ref_len == 0:
        anchor = edit.start - 1

        if end < anchor:
            return start, end
        if start <= anchor <= end:
            return start, end + delta
        if start > anchor:
            return start + delta, end + delta
        return start, end

    edit_start = edit.start
    edit_end = edit.end

    if end < edit_start:
        return start, end

    if start > edit_end:
        return start + delta, end + delta

    # Build projected pieces and return their span. This preserves a single GFF
    # interval while shrinking/expanding across the edit.
    pieces: list[tuple[int, int]] = []

    # Feature portion before the edited reference interval.
    if start < edit_start:
        pieces.append((start, min(end, edit_start - 1)))

    # Replacement/query bases if this edit has query sequence.
    if edit.alt_len > 0:
        pieces.append((edit_start, edit_start + edit.alt_len - 1))

    # Feature portion after the edited reference interval.
    if end > edit_end:
        after_start_old = max(start, edit_end + 1)
        after_end_old = end
        pieces.append((after_start_old + delta, after_end_old + delta))

    if not pieces:
        return None

    new_start = min(p[0] for p in pieces)
    new_end = max(p[1] for p in pieces)

    if new_start < 1 or new_end < new_start:
        return None

    return new_start, new_end


def project_interval(start: int, end: int, projection: ChainProjection) -> Optional[tuple[int, int]]:
    """Project one original target interval through all chain gap edits."""
    current: Optional[tuple[int, int]] = (
        start + projection.initial_offset,
        end + projection.initial_offset,
    )

    cumulative_delta = projection.initial_offset
    for original_edit in projection.edits:
        if current is None:
            return None

        current_edit = GapEdit(
            start=original_edit.start + cumulative_delta,
            ref_len=original_edit.ref_len,
            alt_len=original_edit.alt_len,
        )
        current = patch_interval_once(current[0], current[1], current_edit)
        cumulative_delta += original_edit.delta

    return current


def parse_attrs(attrs: str) -> dict[str, str]:
    result: dict[str, str] = {}
    for raw in attrs.split(";"):
        item = raw.strip()
        if not item or "=" not in item:
            continue
        key, value = item.split("=", 1)
        result[key] = value
    return result


def transcript_id_for_cds(fields: list[str]) -> str:
    attrs = parse_attrs(fields[8])
    return attrs.get("Parent") or attrs.get("transcript_id") or fields[8]


def recalc_cds_phases(lines: list[str]) -> list[str]:
    parsed: list[tuple[str, Optional[list[str]]]] = []
    cds_by_tx: dict[str, list[int]] = defaultdict(list)

    for idx, line in enumerate(lines):
        if line.startswith("#"):
            parsed.append((line, None))
            continue

        fields = line.rstrip("\n").split("\t")
        if len(fields) != 9:
            parsed.append((line, None))
            continue

        parsed.append((line, fields))
        if fields[2] == "CDS":
            cds_by_tx[transcript_id_for_cds(fields)].append(idx)

    for _tx, idxs in cds_by_tx.items():
        if not idxs:
            continue
        first_fields = parsed[idxs[0]][1]
        assert first_fields is not None
        strand = first_fields[6]

        if strand == "-":
            idxs.sort(
                key=lambda i: (int(parsed[i][1][4]), int(parsed[i][1][3])),  # type: ignore[index]
                reverse=True,
            )
        else:
            idxs.sort(
                key=lambda i: (int(parsed[i][1][3]), int(parsed[i][1][4]))  # type: ignore[index]
            )

        cds_bases_seen = 0
        for i in idxs:
            fields = parsed[i][1]
            assert fields is not None
            length = int(fields[4]) - int(fields[3]) + 1
            fields[7] = str((3 - (cds_bases_seen % 3)) % 3)
            cds_bases_seen += length

    output: list[str] = []
    for original, fields in parsed:
        if fields is None:
            output.append(original)
        else:
            output.append("\t".join(fields) + "\n")
    return output


def project_gff(
    gff_in: str,
    chain_path: str,
    gff_out: str,
    recalc_phase: bool = True,
    set_seqname_to_query: bool = False,
) -> tuple[int, int, ChainProjection]:
    projection = parse_chain(chain_path)
    kept = 0
    dropped = 0
    out_lines: list[str] = []

    with open(gff_in, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                out_lines.append(line)
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                out_lines.append(line)
                continue

            if fields[0] != projection.t_name:
                out_lines.append(line)
                kept += 1
                continue

            try:
                start = int(fields[3])
                end = int(fields[4])
            except ValueError:
                out_lines.append(line)
                kept += 1
                continue

            projected = project_interval(start, end, projection)
            if projected is None:
                dropped += 1
                continue

            new_start, new_end = projected
            fields[3] = str(new_start)
            fields[4] = str(new_end)
            if set_seqname_to_query:
                fields[0] = projection.q_name
            out_lines.append("\t".join(fields) + "\n")
            kept += 1

    if recalc_phase:
        out_lines = recalc_cds_phases(out_lines)

    with open(gff_out, "wt") as out:
        out.writelines(out_lines)

    return kept, dropped, projection


def main(argv: Optional[Iterable[str]] = None) -> int:
    parser = argparse.ArgumentParser(
        description="Project GFF through a consensus chain while preserving features that span indels."
    )
    parser.add_argument("chain")
    parser.add_argument("gff_in")
    parser.add_argument("gff_out")
    parser.add_argument(
        "--no-recalc-cds-phase",
        action="store_true",
        help="Do not recalculate CDS phases after coordinate projection.",
    )
    parser.add_argument(
        "--set-seqname-to-query",
        action="store_true",
        help="Set column 1 of projected rows to the query name from the chain header.",
    )
    parser.add_argument("--quiet", action="store_true", help="Suppress summary on stderr.")

    args = parser.parse_args(argv)

    kept, dropped, projection = project_gff(
        args.gff_in,
        args.chain,
        args.gff_out,
        recalc_phase=not args.no_recalc_cds_phase,
        set_seqname_to_query=args.set_seqname_to_query,
    )

    if not args.quiet:
        n_ins = sum(1 for e in projection.edits if e.ref_len == 0 and e.alt_len > 0)
        n_del = sum(1 for e in projection.edits if e.ref_len > 0 and e.alt_len == 0)
        n_repl = sum(1 for e in projection.edits if e.ref_len > 0 and e.alt_len > 0)
        sys.stderr.write(
            f"Projected GFF through chain: target={projection.t_name}, query={projection.q_name}, "
            f"edits={len(projection.edits)} "
            f"(insertions={n_ins}, deletions={n_del}, replacements={n_repl}), "
            f"kept={kept}, dropped={dropped}\n"
        )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
