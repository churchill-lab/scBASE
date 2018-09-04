# -*- coding: utf-8 -*-

import logging


logging.basicConfig(format='[scBASE::%(funcName)s][%(asctime)s] %(message)s', datefmt='%m/%d/%Y %I:%M:%S%p')


def get_logger():
    """
    Get the :class:`logging.Logger`.

    :return: :class:`logging.Logger`
    """
    return logging.getLogger(__name__)


def configure_logging(level):
    """
    Configure the :class:`Logger`.

    - 0 = WARNING
    - 1 = INFO
    - 2 = DEBUG

    :func:`get_logger`.

    :param int level: logging level
    :return: None
    """
    if level == 0:
        get_logger().setLevel(logging.WARN)
    elif level == 1:
        get_logger().setLevel(logging.INFO)
    elif level > 1:
        get_logger().setLevel(logging.DEBUG)


def format_time(start, end):
    """
    Format length of time between ``start`` and ``end``.

    :param start: the start time
    :param end: the end time
    :return: a formatted string of hours, minutes, and seconds
    """
    hours, rem = divmod(end - start, 3600)
    minutes, seconds = divmod(rem, 60)
    return "{:0>2}:{:0>2}:{:05.2f}".format(int(hours), int(minutes), seconds)
