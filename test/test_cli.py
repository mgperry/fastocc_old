import os

import pytest
from click.testing import CliRunner
from fastocc.fastocc import fastocc


@pytest.fixture(autouse=True, scope="class")
def cleanup_files():
    files = ["test.h5wig", "test_nfrs.bed", "test_occ.bw"]
    yield
    for file in files:
        if os.path.exists(file):
            os.remove(file)


class TestCli:
    def test_fastocc_call(self):
        if os.path.exists("test.h5wig"):
            os.remove("test.h5wig")

        runner = CliRunner()
        result = runner.invoke(fastocc, ["call", "--bed", "test.bed", "--bam", "test.bam", "--output", "test.h5wig"])

        assert os.path.exists("test.h5wig")
        assert result.exit_code == 0

    def test_fastocc_nfr(self):
        if not os.path.exists("test.h5wig"):
            pytest.skip("file 'test.h5wig' required to test `fastocc nfr`")

        runner = CliRunner()
        result = runner.invoke(fastocc, ["nfr", "--fragments", "test.h5wig", "--output", "test_nfrs.bed"])

        assert os.path.exists("test_nfrs.bed")
        assert result.exit_code == 0

    def test_fastocc_export(self):
        if not os.path.exists("test.h5wig"):
            pytest.skip("file 'test.h5wig' required to test `fastocc export`")

        runner = CliRunner()
        result = runner.invoke(fastocc, ["export", "--fragments", "test.h5wig", "--output", "test_occ.bw"])

        assert os.path.exists("test_occ.bw")
        assert result.exit_code == 0


def test_fastocc_nps():
    if os.path.exists("test.h5wig"):
        os.remove("test.h5wig")

    runner = CliRunner()
    result = runner.invoke(fastocc, ["call", "--bed", "test.bed", "--bam", "test.bam", "--output", "test.h5wig", "--nps"])

    assert os.path.exists("test.h5wig")
    assert result.exit_code == 0


def test_fastocc_opts():
    if os.path.exists("test.h5wig"):
        os.remove("test.h5wig")

    runner = CliRunner()
    result = runner.invoke(fastocc, ["call", "--bed", "test.bed", "--bam", "test.bam", "--output", "test.h5wig", "--options", "test_opts.json"])

    assert os.path.exists("test.h5wig")
    assert result.exit_code == 0
