#!/usr/bin/env python3
"""Python scaffold pipeline for MFannot."""

from __future__ import annotations

import argparse
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

from tools.mt_annotation_validation import ValidationError, run_table2asn_validation


class PipelineError(RuntimeError):
    """Raised when scaffold pipeline validation or execution fails."""


@dataclass
class PipelineStep:
    name: str
    command: list[str]
    reason: str


@dataclass
class PipelineConfig:
    input_masterfile: Path
    output_dir: Path
    genome_fasta: Path | None
    execute: bool
    min_ani_for_liftoff: float
    enable_liftoff_first_pass: bool
    reference_fasta: Path | None
    liftoff_gff: Path | None
    genetic_code: int
    flip_bin: str
    blast_bin: str
    exonerate_bin: str
    hmmer_bin: str
    erpin_bin: str
    liftoff_bin: str
    fastani_bin: str
    table2asn_bin: str
    submit_tbl: Path | None
    submit_sbt: Path | None


def _validate_file(path: Path, label: str, suffixes: Iterable[str] | None = None) -> Path:
    if not path.exists():
        raise PipelineError(f"{label} does not exist: {path}")
    if not path.is_file():
        raise PipelineError(f"{label} is not a file: {path}")
    if path.stat().st_size == 0:
        raise PipelineError(f"{label} is empty: {path}")
    if suffixes:
        allowed = {s.lower() for s in suffixes}
        if path.suffix.lower() not in allowed:
            raise PipelineError(f"{label} must have one of {sorted(allowed)} extensions: {path}")
    return path.resolve()


def _validate_bin(bin_name: str) -> None:
    if shutil.which(bin_name) is None:
        raise PipelineError(f"Required executable was not found on PATH: {bin_name}")


def validate_pipeline_config(config: PipelineConfig) -> PipelineConfig:
    if not (0.0 <= config.min_ani_for_liftoff <= 1.0):
        raise PipelineError("--min-ani-for-liftoff must be between 0.0 and 1.0")

    _validate_file(config.input_masterfile, "Input masterfile")

    if config.genome_fasta is not None:
        config.genome_fasta = _validate_file(
            config.genome_fasta,
            "Genome FASTA",
            suffixes={".fa", ".fna", ".fasta"},
        )

    if config.enable_liftoff_first_pass:
        if config.reference_fasta is None or config.genome_fasta is None:
            raise PipelineError(
                "LiftOff first-pass requires both --reference-fasta and --genome-fasta"
            )
        config.reference_fasta = _validate_file(
            config.reference_fasta,
            "LiftOff reference FASTA",
            suffixes={".fa", ".fna", ".fasta"},
        )
        if config.liftoff_gff is None:
            raise PipelineError("LiftOff first-pass requires --liftoff-gff")
        config.liftoff_gff = _validate_file(
            config.liftoff_gff,
            "LiftOff reference GFF",
            suffixes={".gff", ".gff3"},
        )

    if config.submit_tbl is not None:
        config.submit_tbl = _validate_file(config.submit_tbl, "Submission TBL", suffixes={".tbl"})
    if config.submit_sbt is not None:
        config.submit_sbt = _validate_file(config.submit_sbt, "Submission SBT", suffixes={".sbt", ".asn"})

    config.output_dir.mkdir(parents=True, exist_ok=True)

    if config.execute:
        core_bins = (
            config.flip_bin,
            config.blast_bin,
            config.exonerate_bin,
            config.hmmer_bin,
            config.erpin_bin,
        )
        for required_bin in core_bins:
            _validate_bin(required_bin)
        if config.enable_liftoff_first_pass:
            _validate_bin(config.liftoff_bin)
            _validate_bin(config.fastani_bin)
        if config.submit_tbl and config.submit_sbt:
            _validate_bin(config.table2asn_bin)

    return config


def parse_fastani_stdout(stdout: str) -> float:
    """Parse ANI value from fastANI stdout and normalize to 0.0-1.0 scale."""
    for line in stdout.splitlines():
        fields = line.strip().split("\t")
        if len(fields) < 3:
            continue
        try:
            value = float(fields[2])
        except ValueError:
            continue
        if value > 1.0:
            value = value / 100.0
        if 0.0 <= value <= 1.0:
            return value
    raise PipelineError("Could not parse ANI value from fastANI output")


def compute_ani(query_fasta: Path, reference_fasta: Path, fastani_bin: str) -> float:
    proc = subprocess.run(
        [fastani_bin, "--query", str(query_fasta), "--ref", str(reference_fasta), "--stdout"],
        check=False,
        capture_output=True,
        text=True,
    )
    if proc.returncode != 0:
        raise PipelineError(f"fastANI failed: {proc.stderr.strip() or proc.stdout.strip()}")
    return parse_fastani_stdout(proc.stdout)


def should_run_liftoff(ani_value: float, threshold: float) -> bool:
    return ani_value >= threshold


def build_pipeline_steps(config: PipelineConfig, ani_value: float | None = None) -> list[PipelineStep]:
    steps: list[PipelineStep] = []

    if config.enable_liftoff_first_pass:
        if ani_value is None:
            raise PipelineError("ANI value required when LiftOff first-pass is enabled")
        if should_run_liftoff(ani_value, config.min_ani_for_liftoff):
            cmd = [
                config.liftoff_bin,
                str(config.genome_fasta),
                str(config.reference_fasta),
                "-g",
                str(config.liftoff_gff),
                "-o",
                str(config.output_dir / "liftoff_first_pass.gff3"),
            ]
            steps.append(
                PipelineStep(
                    name="liftoff_first_pass",
                    command=cmd,
                    reason=f"ANI={ani_value:.4f} meets threshold {config.min_ani_for_liftoff:.4f}",
                )
            )
        else:
            steps.append(
                PipelineStep(
                    name="liftoff_skipped",
                    command=[],
                    reason=f"ANI={ani_value:.4f} below threshold {config.min_ani_for_liftoff:.4f}",
                )
            )

    flip_output = config.output_dir / "flip_output.pep"
    blast_output = config.output_dir / "blast_output.tsv"
    hmmer_output = config.output_dir / "hmmer.tbl"

    steps.extend(
        [
            PipelineStep(
                "flip_orf_detection",
                [config.flip_bin, str(config.input_masterfile), str(flip_output)],
                "Generate ORFs",
            ),
            PipelineStep(
                "blast_first_pass",
                [config.blast_bin, "-query", str(flip_output), "-out", str(blast_output)],
                "Initial protein homology scan",
            ),
            PipelineStep(
                "exonerate_refinement",
                [config.exonerate_bin, "--model", "protein2genome", "--showtargetgff", "yes", "--showalignment", "no"],
                "Intron/exon structure refinement",
            ),
            PipelineStep(
                "hmmer_rna_scan",
                [config.hmmer_bin, "--tblout", str(hmmer_output), str(config.input_masterfile)],
                "Detect RNA and conserved motifs",
            ),
            PipelineStep(
                "erpin_splice_scan",
                [config.erpin_bin, str(config.input_masterfile)],
                "Splice motif structure checks",
            ),
        ]
    )

    return steps


def execute_pipeline_steps(steps: list[PipelineStep], execute: bool) -> None:
    for step in steps:
        if not step.command:
            print(f"[skip] {step.name}: {step.reason}")
            continue

        print(f"[plan] {step.name}: {' '.join(step.command)}")
        if not execute:
            continue

        proc = subprocess.run(step.command, check=False, text=True, capture_output=True)
        if proc.returncode != 0:
            raise PipelineError(
                f"Step '{step.name}' failed with code {proc.returncode}: "
                f"{proc.stderr.strip() or proc.stdout.strip()}"
            )


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Python scaffold for the MFannot pipeline")
    parser.add_argument("masterfile", help="Input masterfile path")
    parser.add_argument("--output-dir", default="mfannot_py_out", help="Output directory")
    parser.add_argument("--genome-fasta", help="Genome FASTA path used for ANI/LiftOff checks")
    parser.add_argument("--execute", action="store_true", help="Execute commands instead of dry-run planning")
    parser.add_argument("--genetic-code", type=int, default=1, help="Genetic code id")

    parser.add_argument("--enable-liftoff-first-pass", action="store_true", help="Enable LiftOff pre-annotation path")
    parser.add_argument("--reference-fasta", help="Reference FASTA for ANI and LiftOff")
    parser.add_argument("--liftoff-gff", help="Reference GFF/GFF3 for LiftOff")
    parser.add_argument("--min-ani-for-liftoff", type=float, default=0.95, help="ANI threshold to run LiftOff (0.0-1.0)")

    parser.add_argument("--flip-bin", default="flip")
    parser.add_argument("--blast-bin", default="blastp")
    parser.add_argument("--exonerate-bin", default="exonerate")
    parser.add_argument("--hmmer-bin", default="hmmsearch")
    parser.add_argument("--erpin-bin", default="erpin")
    parser.add_argument("--liftoff-bin", default="liftoff")
    parser.add_argument("--fastani-bin", default="fastANI")

    parser.add_argument("--tbl", help="Submission feature table for table2asn validation")
    parser.add_argument("--sbt", help="Submission template for table2asn validation")
    parser.add_argument("--table2asn-bin", default="table2asn")

    return parser


def config_from_args(argv: list[str] | None = None) -> PipelineConfig:
    args = _build_parser().parse_args(argv)
    return PipelineConfig(
        input_masterfile=Path(args.masterfile),
        output_dir=Path(args.output_dir),
        genome_fasta=Path(args.genome_fasta) if args.genome_fasta else None,
        execute=bool(args.execute),
        min_ani_for_liftoff=args.min_ani_for_liftoff,
        enable_liftoff_first_pass=bool(args.enable_liftoff_first_pass),
        reference_fasta=Path(args.reference_fasta) if args.reference_fasta else None,
        liftoff_gff=Path(args.liftoff_gff) if args.liftoff_gff else None,
        genetic_code=args.genetic_code,
        flip_bin=args.flip_bin,
        blast_bin=args.blast_bin,
        exonerate_bin=args.exonerate_bin,
        hmmer_bin=args.hmmer_bin,
        erpin_bin=args.erpin_bin,
        liftoff_bin=args.liftoff_bin,
        fastani_bin=args.fastani_bin,
        table2asn_bin=args.table2asn_bin,
        submit_tbl=Path(args.tbl) if args.tbl else None,
        submit_sbt=Path(args.sbt) if args.sbt else None,
    )


def run_scaffold(argv: list[str] | None = None) -> int:
    config = validate_pipeline_config(config_from_args(argv))

    ani_value = None
    if config.enable_liftoff_first_pass:
        if config.genome_fasta is None or config.reference_fasta is None:
            raise PipelineError("LiftOff first-pass requires --genome-fasta and --reference-fasta")
        if config.execute:
            ani_value = compute_ani(config.genome_fasta, config.reference_fasta, config.fastani_bin)
            print(f"[info] ANI between query and reference: {ani_value:.4f}")
        else:
            ani_value = config.min_ani_for_liftoff
            print(
                "[plan] Dry-run mode uses ANI threshold as placeholder for planning. "
                "Use --execute for real ANI gating."
            )

    steps = build_pipeline_steps(config, ani_value=ani_value)
    execute_pipeline_steps(steps, execute=config.execute)

    if config.submit_tbl and config.submit_sbt:
        if config.genome_fasta is None:
            raise PipelineError("--genome-fasta is required when --tbl/--sbt validation is requested")

        try:
            run_table2asn_validation(
                fasta_path=str(config.genome_fasta),
                feature_table_path=str(config.submit_tbl),
                submit_block_path=str(config.submit_sbt),
                table2asn_bin=config.table2asn_bin,
            )
            print("[ok] table2asn validation passed")
        except ValidationError as exc:
            raise PipelineError(str(exc)) from exc

    return 0
