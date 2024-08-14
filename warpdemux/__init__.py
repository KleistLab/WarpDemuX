"""
WarpDemuX

Copyright (c) 2023 by Wiep K. van der Toorn
Contact: w.vandertoorn@fu-berlin.de

"""

import sys
import os

here = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(here, "adapted"))

from warpdemux._version import __version__

__all__ = ["__version__"]
