"""
COBRAk is a general COBRA/COBRA-k suite written as Python package. For more about it, visit its repository:
https://github.com/klamt-lab/COBRAk

The __init__ file of COBRAk initializes the rich text output & tracebacks as well as its logger. Furthermore, graceful
shutdown of user-induced shutdowns is enabled.
"""

import logging
import signal
import sys
from types import FrameType
from typing import Any, Optional

from rich.console import Console
from rich.logging import RichHandler
from rich.traceback import install


# PUBLIC FUNCTIONS SECTION #
def exit_signal_handler(sig: int, frame: FrameType | None) -> None:
    """Handles the exit signal by printing a shutdown message and exiting the program.

    Args:
        sig (int): The signal number.
        frame (Optional[FrameType]): The current frame.
    """
    print(
        "COBRAk received user signal to terminate (this message may appear multiple times in parallelized contexts). Shutting down..."
    )
    sys.exit(0)


def setup_rich_tracebacks(show_locals: bool) -> None:
    """Sets up rich tracebacks with the given options.

    Args:
        show_locals (bool): Whether to show local variables in the traceback.
    """
    install(show_locals=show_locals)


def set_logging_level(level: int) -> None:
    """Sets the logging level.

    E.g. INFO, ERROR, WARNING and CRITICAL from Python's logging module.

    Args:
        level (int): The logging level.
    """
    logger.setLevel(level)


def set_logging_handler(
    show_path: bool = False,
    show_time: bool = False,
    show_level: bool = True,
    keywords: list[str] = [
        "info",
        "warning",
        "error",
        "critical",
    ],
    **args: Any,
) -> RichHandler:
    """
    Sets up the logging handler with the given options.

    Args:
        show_path (bool, optional): Whether to show the path. Defaults to False.
        show_time (bool, optional): Whether to show the time. Defaults to False.
        show_level (bool, optional): Whether to show the level. Defaults to True.
        keywords (Dict[str, str], optional): The keywords to highlight. Defaults to ["info", "warning", "error", "critical"]
        **args (Any, optional): Additional Rich handler arguments.
    """
    return RichHandler(
        show_path=show_path,
        show_time=show_time,
        show_level=show_level,
        keywords=keywords,
        **args,
    )


# ACTIVE SCRIPT SECTION #
signal.signal(signal.SIGINT, exit_signal_handler)
signal.signal(signal.SIGTERM, exit_signal_handler)

setup_rich_tracebacks(False)

# Create a logger
logger = logging.getLogger(__name__)

# Set the logging level
set_logging_level(logging.INFO)

# Add the handler to the logger
logger.addHandler(set_logging_handler())

console = Console()
