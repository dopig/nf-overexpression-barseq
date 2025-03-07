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


def logprint(message: str, print_out: bool = None, debut: str = None, message_type: str = None) -> None:
    """
    Logs a message and optionally prints it to the console with an optional prefix.

    :param message: The message to log and potentially print.
    :param print_out: A boolean indicating whether to print the message.
    :param debut: An optional prefix to print before the message if print_out is True.
    :param message_type: Default is logging.info, but if value is 'error', uses logging.error.
    """

    if message_type == 'error':
        logging.error(message)
    else:
        logging.info(message)

    if print_out:
        # Print the message to the console with an optional prefix
        print(debut if debut is not None else '', message)
