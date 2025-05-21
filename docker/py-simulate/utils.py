import logging
import sys
from pathlib import Path
from pprint import pformat
from typing import TypedDict

class Plasmid(TypedDict):
    id: int
    barcode_sequence: str
    direction: str
    start: int
    end: int

DEFAULT_LOG_PATH = Path("logs/project.log")
HEADER_WIDTH = 47

def setup_logging(message: str | None = None, file_path: Path | None = None) -> None:
    """
    Initializes or reinitializes the application's logging system.

    This function sets up a dual logging approach:
    1.  **Detailed file logging:** All messages (DEBUG level and above) are saved to a file,
        providing a complete record for debugging and analysis.
    2.  **Concise console output:** Only important messages (INFO level and above) are printed
        to the console, offering real-time progress updates without clutter.

    It ensures the log directory exists and prevents duplicate log entries if called multiple times.

    Args:
        message: An optional string to display as an initial banner in the log.
                 Defaults to "Initiating logging" if the log file doesn't exist,
                 or "Resuming logging" if it does.
        file_path: An optional path for the log file. Defaults to "logs/project.log".
    """
    if file_path is None:
        file_path = DEFAULT_LOG_PATH

    if message is None:
        if file_path.exists():
            message = "Resuming logging"
        else:
            message = "Initiating logging"

    file_path.parent.mkdir(parents=True, exist_ok=True)

    # Create a file handler for all levels
    file_handler = logging.FileHandler(file_path)
    file_handler.setLevel(logging.DEBUG)  # Log all levels to file
    file_formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s", datefmt="%Y-%m-%d %H:%M:%S")
    file_handler.setFormatter(file_formatter)

    # Create a stream handler for INFO and above, so they are logged to console
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

    # Add a header to the log file

    message = message.strip()
    if len(message) < HEADER_WIDTH - 4:
        padding = HEADER_WIDTH - 2 - len(message)
        message = "*" * (padding // 2) + " " + message + " " + "*" * (padding - padding // 2)
    logging.debug("*" * HEADER_WIDTH)
    logging.debug(message)

def format_output(result_string: str, variety: str | None = None) -> str:
    """
    Formats multi-line strings for improved readability within log entries.

    This function prepends a consistent indentation or prefix to each non-empty
    line of the input string, making blocks of output (like subprocess results)
    easier to distinguish and read in the log file.

    Args:
        result_string: The raw multi-line string to format.
        variety: An optional string ("wide") to choose a different prefix style.

    Returns:
        The formatted string with prefixes applied to each line.
    """
    if variety == "wide":
        opening = "           :::::::: - "
    else:
        opening = "  :: "
    return "\n".join([opening + line for line in result_string.splitlines() if line.strip() != ""])

def log_args(arg_dict: dict) -> None:
    """
    Logs the command-line arguments.

    The argument dictionary is first pretty-printed for readability, then
    formatted with a line prefix using `format_output`, and finally logged
    at the DEBUG level. This ensures a clear, detailed record of inputs
    is present in the log file.

    Args:
        arg_dict: A dictionary containing the command-line arguments.
    """
    nice_args = pformat(arg_dict, indent=0).strip("{}")
    nicer_args = format_output(nice_args, variety="wide")
    logging.debug(f"Running with arguments:\n{nicer_args}")

def begin_log(message: str | None, log_path: Path, arg_dict: dict) -> None:
    """
    Sets up logging and saves initial information to the log file.
    """
    setup_logging("Simulating Barseq Reads", file_path=log_path)
    logging.info(f"Saving all output, including a more detailed version of this log, to {log_path}")
    raw_command_line = " ".join(sys.argv)
    logging.debug(f"Command run: {raw_command_line}")
    log_args(arg_dict)
