# Copyright 2020 Jacob D. Durrant
#
# Licensed under the Apache License, Version 2.0 (the "License"); you may not
# use this file except in compliance with the License. You may obtain a copy
# of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS, WITHOUT
# WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the
# License for the specific language governing permissions and limitations
# under the License.

from textwrap import fill


def msg(txt, params):
    """Logs a message.

    :param txt: The text to log.
    :type txt: str
    :param params: A dictionary of the user parameters and filters.
    :type params: dict
    """

    msg = "\nMESSAGE: " + txt
    if not params["internal_test"]:
        print(msg)
    params["debug_msg"] = params["debug_msg"] + msg


def warn(txt, params):
    """Logs a warning.

    :param txt: The text to log.
    :type txt: str
    :param params: A dictionary of the user parameters and filters.
    :type params: dict
    """

    msg = "\nWARNING: " + txt
    if not params["internal_test"]:
        print(msg)
    params["debug_warn"] = params["debug_warn"] + msg


def error(txt, params):
    """Logs an error.

    :param txt: The text to log.
    :type txt: str
    :param params: A dictionary of the user parameters and filters.
    :type params: dict
    """

    err_msg = "\n" + fill("ERROR: " + txt)
    if not params["internal_test"]:
        print(err_msg)
    params["debug_error"] = params["debug_error"] + err_msg

    if not params["internal_test"]:
        raise Exception(err_msg)
