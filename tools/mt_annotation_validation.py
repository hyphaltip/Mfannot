#!/usr/bin/env python3
"""Validation helpers for mitochondria annotation submissions."""

from __future__ import annotations

import argparse
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


class ValidationError(RuntimeError):
    """Raised when validation preconditions or validation execution fails."""


@dataclass
class Table2AsnValidationResult:
    command: list[str]
    returncode: int
    stdout: str
    stderr: str
    validation_errors: list[str]

    @property
    def ok(self) -> bool:
        return self.returncode == 0 and not self.validation_errors


def _validate_input_file(path: str, label: str, suffixes: Iterable[str]) -> Path:
    file_path = Path(path)
    if not file_path.exists():
        raise ValidationError(f"{label} file does not exist: {file_path}")
    if not file_path.is_file():
        raise ValidationError(f"{label} path is not a file: {file_path}")
    if file_path.stat().st_size == 0:
        raise ValidationError(f"{label} file is empty: {file_path}")

    allowed = {suffix.lower() for suffix in suffixes}
    if allowed and file_path.suffix.lower() not in allowed:
        raise ValidationError(
            f"{label} file extension must be one of {sorted(allowed)}: {file_path}"
        )
    return file_path.resolve()


def build_table2asn_command(
    fasta_path: str,
    feature_table_path: str,
    submit_block_path: str,
    table2asn_bin: str = "table2asn",
) -> list[str]:
    fasta = _validate_input_file(fasta_path, "FASTA", {".fa", ".fna", ".fasta"})
    feature_table = _validate_input_file(feature_table_path, "Feature table", {".tbl"})
    submit_block = _validate_input_file(submit_block_path, "Submit-block", {".sbt", ".asn"})

    return [
        table2asn_bin,
        "-i",
        str(fasta),
        "-f",
        str(feature_table),
        "-t",
        str(submit_block),
        "-V",
        "vb",
    ]


def extract_validation_errors(stdout: str, stderr: str) -> list[str]:
    error_lines: list[str] = []
    for line in (stdout + "\n" + stderr).splitlines():
        upper = line.upper()
        if "ERROR" in upper or "FATAL" in upper:
            error_lines.append(line.strip())
    return error_lines


def run_table2asn_validation(
    fasta_path: str,
    feature_table_path: str,
    submit_block_path: str,
    table2asn_bin: str = "table2asn",
    timeout_seconds: int = 300,
) -> Table2AsnValidationResult:
    if shutil.which(table2asn_bin) is None:
        raise ValidationError(
            f"table2asn executable not found: {table2asn_bin}"
        )

    command = build_table2asn_command(
        fasta_path=fasta_path,
        feature_table_path=feature_table_path,
        submit_block_path=submit_block_path,
        table2asn_bin=table2asn_bin,
    )

    process = subprocess.run(
        command,
        check=False,
        capture_output=True,
        text=True,
        timeout=timeout_seconds,
    )

    errors = extract_validation_errors(process.stdout, process.stderr)
    result = Table2AsnValidationResult(
        command=command,
        returncode=process.returncode,
        stdout=process.stdout,
        stderr=process.stderr,
        validation_errors=errors,
    )

    if not result.ok:
        problem_lines = errors[:5]
        raise ValidationError(
            "table2asn validation failed. "
            f"return_code={result.returncode}; "
            f"errors={problem_lines if problem_lines else 'none captured'}"
        )

    return result


def _build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Run table2asn validation on MT annotation outputs."
    )
    parser.add_argument("--fasta", required=True, help="FASTA input file (.fa/.fna/.fasta).")
    parser.add_argument(
        "--tbl", required=True, help="Feature table file generated for submission (.tbl)."
    )
    parser.add_argument("--sbt", required=True, help="NCBI submit-block template (.sbt).")
    parser.add_argument(
        "--table2asn-bin",
        default="table2asn",
        help="Name or absolute path to the table2asn executable.",
    )
    return parser


def main() -> int:
    args = _build_arg_parser().parse_args()
    try:
        run_table2asn_validation(
            fasta_path=args.fasta,
            feature_table_path=args.tbl,
            submit_block_path=args.sbt,
            table2asn_bin=args.table2asn_bin,
        )
    except ValidationError as exc:
        print(f"Validation failed: {exc}")
        return 2

    print("Validation passed: table2asn reported no errors.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
