"""Small module for file utility functions"""

from pathlib import Path

# magic has library typing issues and pylint import issues
import magic  # type: ignore # pylint: disable=import-error


class FileUtils:
    """Grouping of file utility functions"""

    @staticmethod
    def is_gzipped(file_path: Path) -> bool:
        """Check if the file is gzipped"""
        return magic.from_file(file_path, mime=True) == "application/gzip"

    @staticmethod
    def is_readable(file_path: Path) -> bool:
        """Check if the file is readable"""
        mime_type = magic.from_file(file_path, mime=True)
        return mime_type in ("text/plain", "application/gzip")
