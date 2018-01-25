## Impulse Response Function Toolbox

Created by Steven Bernsen  
University of Maine  

---
This toolbox has functionality for computing the Green's function or impulse reponse function via cross correlation and deconvolution. Virtual gathers are also supported and the program to compute them will be added soon. The toolbox is currently under construction and many of the MATLAB functions are being translated to Julia and will become obsolete.



Examples and edits will continue to be added along with more documentation as the project progresses. Check the


The text file, 'station_pairs.txt' is required for station-station cross-
correlation. The columns are:

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[filepath_tx]
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[filepath_rx]
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[location_number]
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[location_x_meters]  

where [location_x_meters] is relative to the source for an active source. For an
arbitrary source location use 'virtual_gathers'.  

Feel free to use any code but please acknowledge me in your work.


---  
## Updates
__README Last Updated__ 01/25/2018  
- Original files pushed
