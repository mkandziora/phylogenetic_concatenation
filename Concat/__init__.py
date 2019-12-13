"""
Concat: automatically concatenate single gene alignments.
Copyright (C) 2019  Martha Kandziora
martha.kandziora@yahoo.com

Package to automatically concatenate single gene alignments and to calculate concatenated phylogeny.

All classes and methods are distributed under the following license.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
"""

import sys
import os
import contextlib


_DEBUG_MK = 0


def debug(msg):
    """short debugging command
    """
    if _DEBUG_MK == 1:
        print(msg)


@contextlib.contextmanager
def cd(path):
    cwd = os.getcwd()
    os.chdir(path)
    try:
        yield
    except:
        print('Exception caught: ', sys.exc_info()[0])
    finally:
        os.chdir(cwd)

