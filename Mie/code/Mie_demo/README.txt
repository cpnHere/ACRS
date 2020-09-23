!!README!!

Running Wiscombe Mie Code

Inside this demo directory you will find several .dat input files for the mie_single_size.exe code.  Everyone in the group should already have this code inserted into their path in their '~/.bashrc' file, but in case you need to do that for some reason or you want to poke around in the source code you can find it in the linked directory Mie_code (which points to '/Users/zhibo/Research/Library/Mie_code/').

The .dat files are unformatted text files that the mie_single_size.exe code reads before executing and creating the netCDF output file.  With the mie_single_size.exe code inserted into your path you can run it from any directory, just so long as you also have the following files (with the same filenames) in that directory.  If you don't have all of these files and have them properly named you can run it but you will get an error.

Each file serves the following purposes:

filename.dat -- States the filename of the output netCDF file.  i.e. output.ncdf
ang.dat --  Scattering Angles, first line denotes number of following lines to read, 
	    all following lines are angles in degrees; hit return between each element
size.dat -- Particle diameters, first line denotes number of following lines to read,
	    all following line are diameters in micrometers; hit return between each element.
wl_ref.dat -- Wavelength and Refractive Index, first line denotes number of following lines to read,
	      all following lines are wavelength in micrometers, blank space, real part of refractive index,
	      blank space, complex part of refractive index.

Once you have all of these files created.  Simply type 'mie_single_size.exe' and it will run on the files in your directory and leave the output in that same directory.
