
The GRASS psd that is part of the ESE distribution is a simple "bare-bones" psd.
It is primarily meant for teaching purposes. The relative simplicity of this psd compared to 
other product psd's helps one to easily understand the structure of a psd and eventually aids
researchers to write their own psd's or modify product psd's suitably.

The cost of simplification of the grass psd implies that it can run only with a limited
set of protocols. One such protocol (grassProtocols) is part of the ESE distribution.
To use grassProtocols create a directory grass in the location /usr/g/bin in the scanner
and copy the grass executables (grass and grass.psd.o) to /usr/g/bin/grass directory.

Some known limitations of the grass psd are:
1) Fractional NEX is not supported.
2) Data Acquisition errors are sometimes seen when the number of slices is greater than 1.
