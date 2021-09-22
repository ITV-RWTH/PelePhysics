#!/usr/bin/env python
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#                               Michael A.G. Aivazis
#                        California Institute of Technology
#                        (C) 1998-2003  All Rights Reserved
#
# <LicenseText>
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#

from __future__ import absolute_import

from .SI import micro, milli, nano, pico, second

#
# Definitions of common time units
# Data taken from Appendix F of Halliday, Resnick, Walker, "Fundamentals of Physics",
#     fourth edition, John Willey and Sons, 1993

picosecond = pico * second
nanosecond = nano * second
microsecond = micro * second
millisecond = milli * second

# aliases

s = second
ps = picosecond
ns = nanosecond
us = microsecond
ms = millisecond

# other common units

minute = 60 * second
hour = 60 * minute
day = 24 * hour
year = 365.25 * day


# version
__id__ = "$Id$"

#
# End of file