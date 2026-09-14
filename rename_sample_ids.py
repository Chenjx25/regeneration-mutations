#!/usr/bin/env python3
"""Replace exact sample IDs in a tab-delimited VCF-like table.

The mapping file must contain two tab-delimited columns named old_id and new_id.
Multiple IDs in one field may be separated by semicolons.
"""

import argparse
import csv
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--mapping", required=True, type=Path)
    parser.add_argument("--column", default="ID", help="Header name or zero-based index")
    return parser.parse_args()


def read_mapping(path):
    with path.open(encoding="utf-8", newline="") as handle:
        rows = csv.DictReader(handle, delimiter="\t")
        if rows.fieldnames is None or not {"old_id", "new_id"}.issubset(rows.fieldnames):
            raise SystemExit("Mapping file must contain old_id and new_id columns")
        mapping = {row["old_id"]: row["new_id"] for row in rows if row["old_id"]}
    if not mapping:
        raise SystemExit("Mapping file is empty")
    return mapping


def replace_field(value, mapping):
    return ";".join(mapping.get(token, token) for token in value.split(";"))


def main():
    args = parse_args()
    mapping = read_mapping(args.mapping)
    column_index = int(args.column) if args.column.isdigit() else None
    header_seen = False

    with args.input.open(encoding="utf-8", newline="") as source, args.output.open(
        "w", encoding="utf-8", newline=""
    ) as destination:
        for raw_line in source:
            line = raw_line.rstrip("\r\n")
            if line.startswith("##"):
                destination.write(line + "\n")
                continue

            fields = line.split("\t")
            if not header_seen and (line.startswith("#") or column_index is None):
                names = [name.lstrip("#") for name in fields]
                if column_index is None:
                    if args.column not in names:
                        raise SystemExit(f"Column not found: {args.column}")
                    column_index = names.index(args.column)
                header_seen = True
                destination.write(line + "\n")
                continue

            if column_index is None or column_index >= len(fields):
                raise SystemExit("Could not determine the requested column")
            fields[column_index] = replace_field(fields[column_index], mapping)
            destination.write("\t".join(fields) + "\n")


if __name__ == "__main__":
    main()

