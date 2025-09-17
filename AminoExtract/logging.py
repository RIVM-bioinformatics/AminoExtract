"""This module contains the logging configuration and the custom ArgumentParser classes."""

import logging
import re as _re
import shutil
import textwrap
from argparse import SUPPRESS, ArgumentParser, RawTextHelpFormatter
from typing import IO, Optional

import rich
from rich.highlighter import NullHighlighter
from rich.logging import RichHandler

# Central logging object using Rich's logging library
FORMAT = "%(message)s"
logging.basicConfig(
    level="NOTSET",
    format=FORMAT,
    datefmt="[%x %X]",
    handlers=[
        RichHandler(
            show_path=False,
            omit_repeated_times=False,
            rich_tracebacks=True,
            markup=True,
            highlighter=NullHighlighter(),
        )
    ],
)
log = logging.getLogger("rich")


# FlexiFormatter is a custom ArgumentParser formatter that:
# 1. Fixes spacing in help text
# 2. Properly handles bullet points and lists
# 3. Improves multi-line help text display with default values
#
# This implementation is based on the 'argparse-formatter' module by Dave Steele:
# Original source: https://github.com/davesteele/argparse_formatter
# Specific file: argparse_formatter/flexi_formatter.py
# Commit: a15d89a99e20b0cad4c389a2aa490c551cef4f9c
#
# The code has been refactored to work better with Python 3 and our specific needs.
class FlexiFormatter(RawTextHelpFormatter):
    """
    Help message formatter which respects paragraphs and bulleted lists.
    Only the name of this class is considered a public API.
    """

    def _split_lines(self, text, width):
        return self._para_reformat(text, width)

    def _fill_text(self, text, width, indent):
        lines = self._para_reformat(text, width)
        return "\n".join(lines)

    def _indents(self, line):
        """Return line indent level and "sub_indent" for bullet list text."""

        indent = len(_re.match(r"( *)", line).group(1))
        if list_match := _re.match(r"( *)(([*\-+>]+|\w+\)|\w+\.) +)", line):
            sub_indent = indent + len(list_match.group(2))
        else:
            sub_indent = indent

        return (indent, sub_indent)

    def _split_paragraphs(self, text):
        """Split text in to paragraphs of like-indented lines."""

        text = textwrap.dedent(text).strip()
        text = _re.sub("\n\n\n+", "\n\n", text)

        last_sub_indent = None
        paragraphs = []
        for line in text.splitlines():
            (indent, sub_indent) = self._indents(line)
            is_text = _re.search(r"[^\s]", line) is not None

            if is_text and indent == sub_indent == last_sub_indent:
                paragraphs[-1] += f" {line}"
            else:
                paragraphs.append(line)

            if is_text:
                last_sub_indent = sub_indent
        return paragraphs

    def _para_reformat(self, text, width):
        """Reformat text, by paragraph."""

        paragraphs = []
        for paragraph in self._split_paragraphs(text):

            (indent, sub_indent) = self._indents(paragraph)

            paragraph = self._whitespace_matcher.sub(" ", paragraph).strip()
            new_paragraphs = textwrap.wrap(
                text=paragraph,
                width=width,
                initial_indent=" " * indent,
                subsequent_indent=" " * sub_indent,
            )

            # Blank lines get eaten by textwrap, put it back with [' ']
            paragraphs.extend(new_paragraphs or [" "])

        return paragraphs


class QuickArgFormatter(FlexiFormatter):
    """
    A FlexiFormatter subclass that improves help text formatting.

    This formatter:
    1. Adjusts the maximum width of help text based on terminal width
    2. Adds default values to help text when not explicitly mentioned
    3. Maintains proper alignment and readability of command-line help output
    """

    def __init__(self, prog):
        term_width = shutil.get_terminal_size().columns
        max_help_position = min(max(24, term_width // 2), 80)
        super().__init__(prog, max_help_position=max_help_position)

    def _get_help_string(self, action):
        """
        If the default value is not SUPPRESS and the word "default" is not in the help text,
        then add the default value to the help text
        Parameters
        ----------
        action
            The action object that is being processed.
        Returns
        -------
            The help text for the action.
        """
        help_text = action.help
        if action.default != SUPPRESS and "default" not in help_text.lower() and action.default is not None:
            help_text += f"\n ([underline]default: {str(action.default)}[/underline])"
        return help_text


class RichParser(ArgumentParser):
    """
    A subclass of `argparse.ArgumentParser` that overrides the `_print_message` method to use
    `rich.print` instead of `print`
    """

    # We ignore the override of the typing of the _file variable because we cannot
    # add the correct typing (SupportsWrite) because it is internally defined.
    # We dont use it anyway.
    def _print_message(self, message: str, _file: Optional[IO[str]] = None) -> None:  # type: ignore[override]
        return rich.print(message)
