#!/usr/bin/env python

import argparse
import csv
import logging
import sys
from collections import Counter
from pathlib import Path

logger = logging.getLogger()


class RowChecker:
    """
    Define a service that can validate and transform each given row.

    Attributes:
        modified (list): A list of dicts, where each dict corresponds to a previously
            validated and transformed row. The order of rows is maintained.

    """

    VALID_FORMATS = (
        ".bam",
        ".bai",
    )
    VALID_STATUSES = (
        "tumour",
        "tumor",
        "control",
        "normal",
    )

    def __init__(
        self,
        patient_col="patient",
        sample_col="sample",
        id_col="id",
        status_col="status",
        bam_col="bam",
        bai_col="bai",
        **kwargs,
    ):
        """
        Initialize the row checker with the expected column names.

        Args:
            patient_col (str): The name of the column that contains the patient name
                (default "patient").
            sample_col (str): The name of the column that contains the sample name
                (default "sample").
            id_col(str): The id of the column that contains the id
                (default "id")
            status_col (str): The name of the column that contains the status of the file
                (default "status")
            bam_col (str): The name of the column that contains the bam file
                BAM file path (default "bam").
            bai_col (str): The name of the column that contains the bai file
                BAI file path (default "bai").

        """
        super().__init__(**kwargs)
        self._patient_col = patient_col
        self._sample_col = sample_col
        self._id_col = id_col
        self._status_col = status_col
        self._bam_col = bam_col
        self._bai_col = bai_col
        self._seen = set()
        self.modified = []

    def validate_and_transform(self, row):
        """
        Perform all validations on the given row and insert the read pairing status.

        Args:
            row (dict): A mapping from column headers (keys) to elements of that row
                (values).

        """
        self._validate_patient(row)
        self._validate_sample(row)
        self._validate_status(row)
        self._validate_bam(row)
        self._validate_bai(row)
        self._seen.add(
            (row[self._patient_col], row[self._sample_col], row[self._bam_col]))
        self.modified.append(row)

    def _validate_patient(self, row):
        """Assert that the sample name exists and convert spaces to underscores."""
        assert len(row[self._patient_col]) > 0, "Patient input is required."
        # Sanitize samples slightly.
        row[self._patient_col] = row[self._patient_col].replace(" ", "_")

    def _validate_sample(self, row):
        """Assert that the sample name exists and convert spaces to underscores."""
        assert len(row[self._sample_col]) > 0 or len(
            row[self._id_col]) > 0, "Sample or ID input is required."
        # Sanitize samples slightly.
        row[self._sample_col] = row[self._sample_col].replace(" ", "_")

    def _validate_status(self, row):
        """Assert that the sample name exists and convert spaces to underscores."""
        assert len(row[self._status_col]) > 0, "Status input is required."
        self._validate_status_format(row)

    def _validate_bam(self, row):
        """Assert that the BAM entry is non-empty and has the right format."""
        assert len(row[self._bam_col]) > 0, "The BAM file is required."
        if len(row[self._bam_col]) > 0:
            self._validate_bam_format(row[self._bam_col])

    def _validate_bai(self, row):
        """Assert that the BAI entry has the right format if it exists."""
        assert len(row[self._bai_col]) > 0, "The BAM file is required."
        if len(row[self._bai_col]) > 0:
            self._validate_bam_format(row[self._bai_col])

    def _validate_bam_format(self, filename):
        """Assert that a given filename has one of the expected BAM extensions."""
        assert any(filename.endswith(extension) for extension in self.VALID_FORMATS), (
            f"The BAM file has an unrecognized extension: {filename}\n"
            f"It should be one of: {', '.join(self.VALID_FORMATS)}"
        )

    def _validate_status_format(self, row):
        """Assert that a status is provided."""
        assert any(row[self._status_col] == statuses for statuses in self.VALID_STATUSES), (
            f"Unrecognized status: {row[self._status_col]}\n"
            f"It should be one of: {', '.join(self.VALID_STATUSES)}"
        )
        if any(row[self._status_col] == statuses for statuses in ["control", "normal"]):
            row[self._status_col] = "normal"
        elif any(row[self._status_col] == statuses for statuses in ["tumour", "tumor"]):
            row[self._status_col] = "tumour"

    def validate_unique_samples(self):
        """
        Assert that the combination of sample name and BAM filename is unique.

        """
        assert len(self._seen) == len(
            self.modified), "The pair of patient and sample name must be unique."


def check_samplesheet(file_in, file_out):
    """
    Check that the tabular samplesheet has the structure expected by nf-core pipelines.

    Validate the general shape of the table, expected columns, and each row. Also add
    an additional column which records whether one or two FASTQ reads were found.

    Args:
        file_in (pathlib.Path): The given tabular samplesheet. The format can be either
            CSV, TSV, or any other format automatically recognized by ``csv.Sniffer``.
        file_out (pathlib.Path): Where the validated and transformed samplesheet should
            be created; always in CSV format.

    Example:
        This function checks that the samplesheet follows the following structure,
        see also the `viral recon samplesheet`_::

            patient,sample,status,bam,bai
            PATIENT_1,SAMPLE_1,tumour,patient_1_sample_1.bam,patient_1_sample_1.bam.bai
            PATIENT_1,SAMPLE_2,tumour,patient_1_sample_2.bam,patient_1_sample_2.bam.bai
            PATIENT_1,SAMPLE_3,control,patient_1_sample_3.bam,patient_1_sample_3.bam.bai
            PATIENT_2,SAMPLE_1,tumour,patient_2_sample_1.bam,patient_2_sample_1.bam.bai
            PATIENT_2,SAMPLE_2,control,patient_2_sample_2.bam,patient_2_sample_2.bam.bai


    .. _viral recon samplesheet:
        https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/samplesheet/samplesheet_test_illumina_amplicon.csv

    """
    required_columns = {"patient", "sample", "status", "bam", "bai"}
    # See https://docs.python.org/3.9/library/csv.html#id3 to read up on `newline=""`.
    with file_in.open(newline="") as in_handle:
        reader = csv.DictReader(in_handle)
        # Validate the existence of the expected header columns.
        if not required_columns.issubset(reader.fieldnames):
            logger.critical(
                f"The sample sheet **must** contain the column headers: {', '.join(required_columns)}.")
            sys.exit(1)
        # Validate each row.
        checker = RowChecker()
        for i, row in enumerate(reader):
            try:
                checker.validate_and_transform(row)
            except AssertionError as error:
                logger.critical(f"{str(error)} On line {i + 2}.")
                sys.exit(1)
        checker.validate_unique_samples()
    header = list(reader.fieldnames)
    # See https://docs.python.org/3.9/library/csv.html#id3 to read up on `newline=""`.
    with file_out.open(mode="w", newline="") as out_handle:
        writer = csv.DictWriter(out_handle, header, delimiter=",")
        writer.writeheader()
        for row in checker.modified:
            writer.writerow(row)


def parse_args(argv=None):
    """Define and immediately parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Validate and transform a tabular samplesheet.",
        epilog="Example: python check_samplesheet.py samplesheet.csv samplesheet.valid.csv",
    )
    parser.add_argument(
        "file_in",
        metavar="FILE_IN",
        type=Path,
        help="Tabular input samplesheet in CSV or TSV format.",
    )
    parser.add_argument(
        "file_out",
        metavar="FILE_OUT",
        type=Path,
        help="Transformed output samplesheet in CSV format.",
    )
    parser.add_argument(
        "-l",
        "--log-level",
        help="The desired log level (default WARNING).",
        choices=("CRITICAL", "ERROR", "WARNING", "INFO", "DEBUG"),
        default="WARNING",
    )
    return parser.parse_args(argv)


def main(argv=None):
    """Coordinate argument parsing and program execution."""
    args = parse_args(argv)
    logging.basicConfig(level=args.log_level,
                        format="[%(levelname)s] %(message)s")
    if not args.file_in.is_file():
        logger.error(f"The given input file {args.file_in} was not found!")
        sys.exit(2)
    args.file_out.parent.mkdir(parents=True, exist_ok=True)
    check_samplesheet(args.file_in, args.file_out)


if __name__ == "__main__":
    sys.exit(main())
