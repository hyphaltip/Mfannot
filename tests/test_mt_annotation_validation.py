import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from tools.mt_annotation_validation import (
    ValidationError,
    build_table2asn_command,
    extract_validation_errors,
    run_table2asn_validation,
)


class BuildCommandTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temp_dir = tempfile.TemporaryDirectory()
        self.base = Path(self.temp_dir.name)

        self.fasta = self.base / "input.fna"
        self.tbl = self.base / "features.tbl"
        self.sbt = self.base / "template.sbt"

        self.fasta.write_text(">contig1\nATGCGT\n", encoding="utf-8")
        self.tbl.write_text(">Feature contig1\n1\t6\tgene\n", encoding="utf-8")
        self.sbt.write_text("Submit-block ::= {}\n", encoding="utf-8")

    def tearDown(self) -> None:
        self.temp_dir.cleanup()

    def test_build_table2asn_command(self) -> None:
        command = build_table2asn_command(
            fasta_path=str(self.fasta),
            feature_table_path=str(self.tbl),
            submit_block_path=str(self.sbt),
        )

        self.assertEqual(command[0], "table2asn")
        self.assertIn("-i", command)
        self.assertIn("-f", command)
        self.assertIn("-t", command)
        self.assertIn("-V", command)

    def test_rejects_missing_file(self) -> None:
        with self.assertRaises(ValidationError):
            build_table2asn_command(
                fasta_path=str(self.base / "missing.fna"),
                feature_table_path=str(self.tbl),
                submit_block_path=str(self.sbt),
            )

    def test_rejects_empty_file(self) -> None:
        self.fasta.write_text("", encoding="utf-8")
        with self.assertRaises(ValidationError):
            build_table2asn_command(
                fasta_path=str(self.fasta),
                feature_table_path=str(self.tbl),
                submit_block_path=str(self.sbt),
            )

    def test_rejects_wrong_extension(self) -> None:
        wrong = self.base / "input.txt"
        wrong.write_text(">x\nATG\n", encoding="utf-8")
        with self.assertRaises(ValidationError):
            build_table2asn_command(
                fasta_path=str(wrong),
                feature_table_path=str(self.tbl),
                submit_block_path=str(self.sbt),
            )


class RunValidationTests(unittest.TestCase):
    @patch("tools.mt_annotation_validation.shutil.which", return_value=None)
    def test_missing_table2asn_binary(self, _which) -> None:
        with self.assertRaises(ValidationError):
            run_table2asn_validation("a.fna", "b.tbl", "c.sbt")

    def test_extract_validation_errors(self) -> None:
        lines = extract_validation_errors(
            stdout="INFO all good\nERROR overlap issue\n",
            stderr="warning only\nFATAL parse failure\n",
        )
        self.assertEqual(lines, ["ERROR overlap issue", "FATAL parse failure"])

    @patch("tools.mt_annotation_validation.subprocess.run")
    @patch("tools.mt_annotation_validation.shutil.which", return_value="/usr/bin/table2asn")
    def test_raises_on_table2asn_error_output(self, _which, run_mock) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            base = Path(tmp)
            fasta = base / "input.fna"
            tbl = base / "features.tbl"
            sbt = base / "template.sbt"
            fasta.write_text(">c\nATG\n", encoding="utf-8")
            tbl.write_text(">Feature c\n1\t3\tgene\n", encoding="utf-8")
            sbt.write_text("Submit-block ::= {}\n", encoding="utf-8")

            run_mock.return_value.returncode = 0
            run_mock.return_value.stdout = "ERROR invalid product"
            run_mock.return_value.stderr = ""

            with self.assertRaises(ValidationError):
                run_table2asn_validation(str(fasta), str(tbl), str(sbt))


if __name__ == "__main__":
    unittest.main()
