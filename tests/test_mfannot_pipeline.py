import tempfile
import unittest
from pathlib import Path

from tools.mfannot_pipeline import (
    PipelineConfig,
    PipelineError,
    build_pipeline_steps,
    parse_fastani_stdout,
    should_run_liftoff,
    validate_pipeline_config,
)


class MfannotPipelineTests(unittest.TestCase):
    def setUp(self) -> None:
        self.tmp = tempfile.TemporaryDirectory()
        self.base = Path(self.tmp.name)

        self.master = self.base / "input.masterfile"
        self.master.write_text(">ctg\nATGC\n", encoding="utf-8")

        self.query = self.base / "query.fna"
        self.ref = self.base / "ref.fna"
        self.gff = self.base / "ref.gff3"

        self.query.write_text(">q\nATGC\n", encoding="utf-8")
        self.ref.write_text(">r\nATGC\n", encoding="utf-8")
        self.gff.write_text("##gff-version 3\n", encoding="utf-8")

    def tearDown(self) -> None:
        self.tmp.cleanup()

    def _config(self, **kwargs) -> PipelineConfig:
        values = dict(
            input_masterfile=self.master,
            output_dir=self.base / "out",
            genome_fasta=self.query,
            execute=False,
            min_ani_for_liftoff=0.95,
            enable_liftoff_first_pass=True,
            reference_fasta=self.ref,
            liftoff_gff=self.gff,
            genetic_code=1,
            flip_bin="flip",
            blast_bin="blastp",
            exonerate_bin="exonerate",
            hmmer_bin="hmmsearch",
            erpin_bin="erpin",
            liftoff_bin="liftoff",
            fastani_bin="fastANI",
            table2asn_bin="table2asn",
            submit_tbl=None,
            submit_sbt=None,
        )
        values.update(kwargs)
        return PipelineConfig(**values)

    def test_validate_config_creates_output_dir(self) -> None:
        cfg = validate_pipeline_config(self._config())
        self.assertTrue(cfg.output_dir.exists())

    def test_invalid_ani_threshold_raises(self) -> None:
        with self.assertRaises(PipelineError):
            validate_pipeline_config(self._config(min_ani_for_liftoff=1.5))

    def test_parse_fastani_stdout_percent(self) -> None:
        parsed = parse_fastani_stdout("q\tr\t98.6\t120\t120\n")
        self.assertAlmostEqual(parsed, 0.986)

    def test_parse_fastani_stdout_fraction(self) -> None:
        parsed = parse_fastani_stdout("q\tr\t0.96\t120\t120\n")
        self.assertAlmostEqual(parsed, 0.96)

    def test_should_run_liftoff(self) -> None:
        self.assertTrue(should_run_liftoff(0.96, 0.95))
        self.assertFalse(should_run_liftoff(0.90, 0.95))

    def test_build_steps_includes_liftoff_when_ani_high(self) -> None:
        steps = build_pipeline_steps(self._config(), ani_value=0.97)
        self.assertEqual(steps[0].name, "liftoff_first_pass")

    def test_build_steps_marks_liftoff_skipped_when_ani_low(self) -> None:
        steps = build_pipeline_steps(self._config(), ani_value=0.80)
        self.assertEqual(steps[0].name, "liftoff_skipped")


if __name__ == "__main__":
    unittest.main()
