from click.testing import CliRunner
from fastocc.fastocc import fastocc

def test_fastocc_cli():
    runner = CliRunner()
    result = runner.invoke(fastocc, ["occ", "--bed", "test.bed", "--bam", "test.bam"])
    assert result.exit_code == 0

