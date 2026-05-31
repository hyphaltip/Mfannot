#!/usr/bin/env python3
"""Python scaffold entrypoint for MFannot."""

from __future__ import annotations

from tools.mfannot_pipeline import PipelineError, run_scaffold


def main() -> int:
    try:
        return run_scaffold()
    except PipelineError as exc:
        print(f"mfannot.py failed: {exc}")
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
