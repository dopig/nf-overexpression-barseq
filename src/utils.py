import time
import logging
from pathlib import Path

DEFAULT_LOG_PATH = Path("logs/project.log")


def setup_logging(message: str = None, file_path: Path = None) -> None:
    """
    Initializes or reinitializes the logging session.

    `message` should be a brief description of the logging session.
    If no path is provided, it uses DEFAULT_LOG_PATH.

    :param message: A brief description of the logging session.
    :param file_path: The path to the log file. If None, uses DEFAULT_LOG_PATH.
    """
    if file_path is None:
        file_path = DEFAULT_LOG_PATH

    if message is None:
        message = "Resuming logging"

    file_path.parent.mkdir(parents=True, exist_ok=True)

    logging.basicConfig(
        filename=file_path,
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S"
    )

    message = message.strip()

    if len(message) < 43:
        padding = 47 - 2 - len(message)
        message = "*" * (padding // 2) + " " + message + " " + "*" * (padding - padding // 2)

    logging.info("*" * 47)
    logging.info(message)

def loginfo(message: str, print_out: bool = None, debut: str = None) -> None:
    """Logs an info message and optionally prints it to the console."""
    logging.info(message)
    if print_out:
        print(debut if debut is not None else '', message)

def logerror(message: str, print_out: bool = None, debut: str = None) -> None:
    """Logs an error message, optionally prints it, and raises an exception."""
    logging.error(message)
    if print_out:
        print(debut if debut is not None else '', message)
    raise RuntimeError(message)
