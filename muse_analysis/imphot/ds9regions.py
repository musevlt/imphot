#!/usr/bin/env python

from __future__ import print_function, division

# This module lets one read circle, ellipse and box regions from a ds9
# region file. Beware that it assumes that the coordinate system of
# the central coordinates of each region is fk5.

import re

__all__ = ['Ds9Regions', 'Ds9Region', 'Ds9Circle', 'Ds9Ellipse', 'Ds9Box']

class Ds9Regions(object):
    """Parse circle, ellipse and box regions from a ds9 region file.

    Only circle, ellipse, and box regions are recognized. Other types
    of ds9 region and configuration parameters are simply ignored. The
    center coordinates of each region are assumed to be in equatorial
    coordinates and sizes are assumed to be angular sizes, not pixel
    sizes.

    The Ds9Regions object can be treated like an indexable tuple. For
    example::

       regions = ds9regions.Ds9Regions("myfile.reg")

       # Indexing and len() are supported.

       print("First region = ", regions[0])
       print("Last region = ", regions[-1])

       # Repeatable iteration like a tuple is supported:

       for region in regions:
          print(region)

    Parameters
    ----------
    regions : filename or iterable
        Either a string that contains the name of a ds9 region file,
        or an iterable that returns successive lines of a ds9 file.

    """

    def __init__(self, regions):

        # Start with an empty list of regions.

        self._regions = []

        # Have we been given the name of a ds9 region file, or an
        # iterable.
        if isinstance(regions, str):
            lines = open(regions)
        else:
            lines = regions

        # Parse the lines of the region file to build up the list of
        # regions.

        region = self._parse_regions(lines)

        # Close the region file if we opened one.

        if isinstance(regions, str):
            lines.close()

    def __str__(self):
        s = ""
        for region in self._regions:
            s += "%s\n" % region.__str__()
        return s

    def __len__(self):
        return len(self._regions)

    def __getitem__(self, index):
        return self._regions[index]

    def _parse_regions(self, lines):

        # Set the default coordinate system to match ds9's default.

        system = "physical"

        # Parse one line at a time.

        for line_number, line in enumerate(lines):
            try:

                # No region has been parsed from this line yet.

                region = None

                # Each line can start with a sequence of configuration
                # specifiers that are separated from each other by
                # semicolons.

                last = ""
                for field in re.split(' *; *', line):

                    # Trim white-space and newlines.

                    field = field.strip()

                    # Check for a coordinate-system specifier.

                    if field.lower() in ["physical","image","fk4","b1950","fk5",
                                         "j2000", "galactic","ecliptic","icrs",
                                         "linear", "amplifier","detector"]:
                        system = field.lower()

                    # Keep a record of the last non-empty field that we haven't
                    # already processed.

                    elif len(field) != 0:
                        last = field

                # Discard all but the last unprocessed field.

                line = last

                # If the line is empty, then there is no region here.

                if len(line) < 1:
                    continue

                # Regions that start with a minus sign are to be excluded.
                # If a minus sign is present, record this and remove it.

                if line[0] == '-':
                    exclude = True
                    line = line[1:]
                else:
                    exclude = False

                # All lines that describe regions specify the shape followed by
                # an open parenthesis

                if re.match('[a-z]+\(', line) is None:
                    continue

                # Split the line into its fields.

                fields = re.split(' *[(),] *', line)

                # All fields have a shape name, followed by the x and y axis
                # coordinates of the center of the region.

                shape = fields[0]
                x = _parse_x(fields[1])
                y = _parse_y(fields[2])

                # Parse the distinct parameters of different shapes.

                if shape == 'circle':
                    region = Ds9Circle(exclude, system, x, y,
                                       radius=_parse_size(fields[3]))
                elif shape == 'ellipse':
                    region = Ds9Ellipse(exclude, system, x, y,
                                        width=_parse_size(fields[3]),
                                        height=_parse_size(fields[4]),
                                        pa=float(fields[5]))
                elif shape == 'box':
                    region = Ds9Box(exclude, system, x, y,
                                    width=_parse_size(fields[3]),
                                    height=_parse_size(fields[4]),
                                    pa=float(fields[5]))

                if region is not None:
                    self._regions.append(region)
            except ValueError as e:
                raise ValueError("Error on line %d (%s)" % (line_number+1, e.message))

class Ds9Region(object):
    """A generic Ds9Region.

    Parameters
    ----------
    shape   : str
       The name of the shape ("circle", "ellipse", or "box").
    exclude : bool
       True if the region is to be excluded.
    system  : str
       The name of the coordinate system.
    x       : float
       The x coordinate of the center of the shape (degrees).
    y       : float
       The y coordinate of the center of the shape (degrees).

    Attributes
    ----------
    shape   : str
       The name of the shape ("circle", "ellipse", or "box").
    exclude : bool
       True if the region is to be excluded.
    system  : str
       The name of the coordinate system.
    x       : float
       The x coordinate of the center of the shape (degrees).
    y       : float
       The y coordinate of the center of the shape (degrees).

    """

    def __init__(self, shape, exclude, system, x, y):
        self.shape = shape
        self.exclude = exclude
        self.system = system
        self.x = x
        self.y = y

    def __str__(self):
        return "%s %s %s x: %.6f  y: %.6f" % ("exclude" if self.exclude else "include", self.system, self.shape, self.x, self.y)


class Ds9Circle(Ds9Region):
    """A circular ds9 region.

    Parameters
    ----------
    exclude : bool
       True if the region is to be excluded.
    system  : str
       The name of the coordinate system.
    x      : float
       The x-axis coordinate of the center of the shape (degrees).
    y     : float
       The declination of the center of the shape (degrees).
    radius  : float
       The radius of the circle (degrees).

    Attributes
    ----------
    shape   : str
       The name of the shape ("circle", "ellipse", or "box").
    exclude : bool
       True if the region is to be excluded.
    system  : str
       The name of the coordinate system.
    x       : float
       The x coordinate of the center of the shape (degrees).
    y       : float
       The y coordinate of the center of the shape (degrees).
    radius  : float
       The radius of the circle (degrees).

    """

    def __init__(self, exclude, system, x, y, radius):
        Ds9Region.__init__(self, "circle", exclude, system, x, y)
        self.radius = radius
    def __str__(self):
        return Ds9Region.__str__(self) + ("  r: %g" % self.radius)

class Ds9Ellipse(Ds9Region):
    """An elliptical ds9 region.

    Parameters
    ----------
    exclude : bool
       True if the region is to be excluded.
    system  : str
       The name of the coordinate system.
    x      : float
       The x-axis coordinate of the center of the shape (degrees).
    y     : float
       The declination of the center of the shape (degrees).
    width   : float
       The width (degrees) of the ellipse when the position angle is 0.
    height  : float
       The height (degrees) of the ellipse when the position angle is 0.
    pa      : float
       The position angle of the rectangle, east of celestial north for
       sky coordinate systems.

    Attributes
    ----------
    shape   : str
       The name of the shape ("circle", "ellipse", or "box").
    exclude : bool
       True if the region is to be excluded.
    system  : str
       The name of the coordinate system.
    x       : float
       The x coordinate of the center of the shape (degrees).
    y       : float
       The y coordinate of the center of the shape (degrees).
    width   : float
       The width (degrees) of the ellipse when the position angle is 0.
    height  : float
       The height (degrees) of the ellipse when the position angle is 0.
    pa      : float
       The position angle of the rectangle, east of celestial north for
       sky coordinate systems.

    """

    def __init__(self, exclude, system, x, y, width, height, pa):
        Ds9Region.__init__(self, "ellipse", exclude, system, x, y)
        self.width = width
        self.height = height
        self.pa = pa
    def __str__(self):
        return Ds9Region.__str__(self) + ("  w: %g  h: %g  PA: %g" % (self.width, self.height, self.pa))

class Ds9Box(Ds9Region):
    """A rectangular ds9 region.

    Parameters
    ----------
    exclude : bool
       True if the region is to be excluded.
    system  : str
       The name of the coordinate system.
    x      : float
       The x-axis coordinate of the center of the shape (degrees).
    y     : float
       The declination of the center of the shape (degrees).
    width   : float
       The width (degrees) of the rectangle when the position angle is 0.
    height  : float
       The height (degrees) of the rectangle when the position angle is 0.
    pa      : float
       The position angle of the rectangle, east of celestial north for
       sky coordinate systems.

    Attributes
    ----------
    shape   : str
       The name of the shape ("circle", "ellipse", or "box").
    exclude : bool
       True if the region is to be excluded.
    system  : str
       The name of the coordinate system.
    x       : float
       The x coordinate of the center of the shape (degrees).
    y       : float
       The y coordinate of the center of the shape (degrees).
    width   : float
       The width (degrees) of the rectangle when the position angle is 0.
    height  : float
       The height (degrees) of the rectangle when the position angle is 0.
    pa      : float
       The position angle of the rectangle, east of celestial north for
       sky coordinate systems.

    """
    def __init__(self, exclude, system, x, y, width, height, pa):
        Ds9Region.__init__(self, "box", exclude, system, x, y)
        self.width = width
        self.height = height
        self.pa = pa
    def __str__(self):
        return Ds9Region.__str__(self) + ("  w: %g  h: %g  PA: %g" % (self.width, self.height, self.pa))

def _parse_x(s):
    """Parse an x-axis coordinate from a string.

    This can be a simple number, or it can be a sexagesimal number of
    hours, with up to three fields (hours, minutes, seconds) that are
    separated by either colons or spaces. In the later case the
    sexagesimal number is converted from sexagesimal hours to decimal
    degrees.

    """
    try:
        return float(s)
    except:
        return _parse_sexagesimal(s) * 15.0

def _parse_y(s):
    """Parse a y-axis coordinate from a string.

    This can be a simple number, or it can be a sexagesimal number of
    degrees, with up to three fields (degrees, minutes, seconds) that
    are separated by either colons or spaces. In the later case the
    sexagesimal number is converted to decimal degrees.

    """
    try:
        return float(s)
    except:
        return _parse_sexagesimal(s)

def _parse_sexagesimal(s):
    """Parse a sexagesimal number from a string. This can have up
    to numeric fields separated by either colons or spaces.
    """

    # Get the sign of the number, and remove any sign character
    # from the string.

    if s[0] == '-' or s[0] == '+':
        sign = -1.0 if s[0] == '-' else 1.0
        s = s[1:]
    else:
        sign = 1.0

    # Split the sexagesimal number into fields.

    numbers = re.split(' *[ :] *', s)

    # Complain if there are less than 1 or more than 3 numeric fields.

    if len(numbers) < 1 or len(numbers) > 3:
        raise ValueError("parse_sexagesimal expects from 1 to 3 numeric fields")

    # Accumulate the angle in the units of the largest term.

    angle = float(numbers[0])
    scale = 1.0 / 60.0
    for n in numbers[1:]:
        angle += float(n) * scale
        scale /= 60.0

    # Scale the result by the sign.

    return sign * angle

def _parse_size(s):

    """Parse a size from a string. This can be a simple number, or
    an angle in arc-minutes followed by "'", or an angle in arc-seconds
    followed by '"'. Arcseconds and arcminutes are converted to degrees.
    """

    if s[-1] == "'":   # Arcmin
        return float(s[:-1]) / 60.0
    elif s[-1] == '"': # Arcsec
        return float(s[:-1]) / 3600.0
    else:
        return float(s)

# When this module is invoked as a standalone program, parse the
# contents of any region file that is specified on the command line
# and print out what is parsed from that file to stdout.

if __name__ == '__main__':
    import sys
    if len(sys.argv) > 1:
        regs = Ds9Regions(sys.argv[1])
        print(regs)
    else:
        print("Usage: %s <filename>" % sys.argv[0])
        exit(1)
