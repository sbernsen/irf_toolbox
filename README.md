## Impulse Response Function Toolbox

Created by Steven Bernsen  
University of Maine  

---
This toolbox has functionality for computing the Green's function or impulse
reponse function via cross correlation and deconvolution. Virtual gathers are
also supported and the program to compute them will be added soon. The toolbox
is currently under construction and many of the MATLAB functions are being
translated to Julia and will become obsolete.

Examples and edits will continue to be added along with more documentation as the project progresses. Check the

#### Usage ####

R and julia programming languages along with Bash need to be installed. To install the required packages, run from a terminal:  

$ bash install_packages

After installing required packages, put the files 'EGF' and 'Preprocessing' into the julia path. This can be done by putting the modules into the directory given by running the command  

\> Pkg.Dir.path()  

or creating a ~/.juliarc.jl file and appending:  @everywhere push!(LOAD_PATH,"/directory/path/to/the/modules/").


The tab delimited text file, 'station_pairs.txt' is required for station-station cross-correlation. The columns are:

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[filepath_tx]
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[filepath_rx]
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[location_number]
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[location_x_meters]  

where [location_x_meters] is relative to the source for an active source. For an
arbitrary source location use 'virtual_gathers'. Since the distance from source
may have duplicates, use 'location_number' to refer to the
relative order of the stations.



Feel free to use any code but please acknowledge me in your work.


---  
## Updates
__README Last Updated__ 01/27/2018  

01/25/2018 - Original files pushed  

01/27/2018 -  Modified frequency-domain cross correlation functions, added time-domain cross correlation, added new functions for preprocessing and 
filtering
