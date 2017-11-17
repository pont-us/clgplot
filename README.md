# clgplot: a utility to plot IRM acquisition data

clgplot reads two types of file: raw isothermal remanent magnetization
(IRM) acquisition data (as a two-column file containing applied field
and measured magnetization), and output files from the IrmUnmix program
[1] which contain a representation of IRM acquisition data as a sum of
cumulative log-Gaussian (CLG) curves. clgplot produces a plot of the
data and/or curves using pyplot. The plot can be viewed on-screen and
saved to a file.

clgplot also calculates and prints H'cr, the applied field which
produces a remanent magnetization half as intense as the saturation
remanent magnetization. For this calculation, the maximum magnetization
value in the input file is assumed to be the saturation remanent
magnetization.

If run with no arguments, clgplot starts in GUI mode. Various options
for command-line usage are also available; supplying the `-h` or
`--help` argument will detail the available command-line arguments.

clgplot is copyright 2012 by Pontus Lurcock (pont -at- talvi -dot- net)
and released under the MIT license (see source code for full license
text).

[1] See http://www.geo.uu.nl/~forth/Software/irmunmix/ and
http://dx.doi.org/10.1046/j.0956-540x.2001.01558.x
