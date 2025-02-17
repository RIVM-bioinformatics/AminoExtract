from pathlib import Path

# magic has library typing issues
import magic  # type: ignore


class FileUtils:
    """File utility functions"""

    @staticmethod
    def is_gzipped(file_path: Path) -> bool:
        return magic.from_file(file_path, mime=True) == "application/gzip"

    @staticmethod
    def is_readable(file_path: Path) -> bool:
        mime_type = magic.from_file(file_path, mime=True)
        return mime_type in ("text/plain", "application/gzip")
