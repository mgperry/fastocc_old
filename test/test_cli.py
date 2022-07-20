from click.testing import CliRunner
from fastocc.fastocc import cli

def test_fastocc_cli():
    runner = CliRunner()
    result = runner.invoke(cli, ["--bed", "test.bed", "--bam", "test.bam"])
    assert result.exit_code == 0

