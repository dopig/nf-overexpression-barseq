import logging
import sys
from pathlib import Path
from pprint import pformat

DEFAULT_LOG_PATH = Path("logs/project.log")


def setup_logging(message: str = None, file_path: Path = None) -> None:
    """
    Initializes or reinitializes the logging session.
    """
    if file_path is None:
        file_path = DEFAULT_LOG_PATH

    if message is None:
        message = "Resuming logging"

    file_path.parent.mkdir(parents=True, exist_ok=True)

    # Create a file handler for all levels
    file_handler = logging.FileHandler(file_path)
    file_handler.setLevel(logging.DEBUG)  # Log all levels to file
    file_formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s", datefmt="%Y-%m-%d %H:%M:%S")
    file_handler.setFormatter(file_formatter)

    # Create a stream handler for INFO and above
    stream_handler = logging.StreamHandler(sys.stdout)
    stream_handler.setLevel(logging.INFO)  # Log INFO and above to console
    stream_formatter = logging.Formatter("%(asctime)s - %(message)s", datefmt="%H:%M:%S")
    stream_handler.setFormatter(stream_formatter)

    # Get the root logger
    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)  # Allow all levels to be processed
    logger.handlers.clear() # removes old handlers if setup_logging is called multiple times.

    # Add the handlers to the logger
    logger.addHandler(file_handler)
    logger.addHandler(stream_handler)

    message = message.strip()

    if len(message) < 43:
        padding = 47 - 2 - len(message)
        message = "*" * (padding // 2) + " " + message + " " + "*" * (padding - padding // 2)

    logging.debug("*" * 47)
    logging.debug(message)


def format_output(result_string: str, variety=None) -> str:
    """
    Formats strings for the log.
    Typical result_string is result.stdout from a module like ccs.
    """
    if variety == "wide":
        opening = "           :::::::: - "
    else:
        opening = "  :: "
    return "\n".join([opening + line for line in result_string.splitlines() if line.strip() != ""])

def log_args(arg_dict: dict) -> None:
    nice_args = pformat(arg_dict, indent=0)[1:-1]
    nicer_args = format_output(nice_args, variety="wide")
    logging.debug(f"Running with arguments:\n{nicer_args}")
