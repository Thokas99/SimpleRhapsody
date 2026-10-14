"""MkDocs hook: replaces the placeholder in cli.md with the output of `rhapsodist --help`."""

import subprocess
import sys

PLACEHOLDER = '<!-- rhapsodist-help -->'


def on_page_markdown(markdown, page, **kwargs):
    if PLACEHOLDER not in markdown:
        return markdown
    help_text = subprocess.run([sys.executable, '-m', 'rhapsodist.cli', '--help'],
                               capture_output=True, text=True, check=True).stdout
    return markdown.replace(PLACEHOLDER, '```\n' + help_text + '```')
