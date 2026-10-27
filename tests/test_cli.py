import subprocess
import sys
import unittest
from types import SimpleNamespace
from unittest.mock import patch

from rhapsodist import cli


class CliTests(unittest.TestCase):
    def snakemake_command(self, *arguments):
        with patch.object(sys, 'argv', ['rhapsodist', *arguments]), \
                patch.object(subprocess, 'run', return_value=SimpleNamespace(returncode=0)) as run:
            with self.assertRaises(SystemExit) as exit_info:
                cli.main()
        self.assertEqual(exit_info.exception.code, 0)
        return run.call_args.args[0]

    def test_config_overrides_are_forwarded_to_snakemake(self):
        cmd = self.snakemake_command('--configfile', 'a.yaml', '--config', 'mode=usa')
        self.assertEqual(cmd[cmd.index('--configfile') + 1], 'a.yaml')
        self.assertEqual(cmd[-2:], ['--config', 'mode=usa'])

    def test_dry_run_and_cores(self):
        cmd = self.snakemake_command('--configfile', 'a.yaml', '--cores', '3', '-n')
        self.assertEqual(cmd[cmd.index('--cores') + 1], '3')
        self.assertIn('--dry-run', cmd)
