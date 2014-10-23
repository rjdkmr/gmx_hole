##gmx_hole
***

###About

To calculate radius of protein channel using hole program for GROMACS MD trajectory. gmx\_hole uses [hole program]                                                       (http://www.csb.yale.edu/userguides/graphics/hole/doc/hole_d00.html). <strong>Please cite the original publication of hole: </strong>                         
O.S. Smart, J.M. Goodfellow and B.A. Wallace (1993)                     
The Pore Dimensions of Gramicidin A                                     
_Biophysical Journal_ 65:2455-2460.

***

###Requirements
* To compile and install, GROMACS libraries <code> libgmx, libmd, libgmxana </code> are required.
* To use <code>gmx_hole</code>, <code>hole</code> program should be already installed.

***

###Download
<pre><code>git clone https://github.com/rjdkmr/gmx_hole
</code></pre>
***

###Installation
<pre><code>cd gmx_hole
mkdir build
cd build
cmake ..  -DGMX_PATH=/opt/gromacs -DCMAKE_INSTALL_PREFIX=/opt/gmxhole
make
make install
</code></pre>

Directory <code>/opt/gromacs</code> should contains <code>include</code> and <code> lib </code> directories. If these directories are in seprate locations, use followings:
<pre><code>cmake ..  -DGMX_LIB=/path/to/lib -DGMX_INCLUDE=/path/to/include -DCMAKE_INSTALL_PREFIX=/opt/gmxhole
</code></pre>

If fftw library <code> libfftw3f.so or libfftw3f.a </code> are not present in standard locations:
<pre><code>-DFFTW_LIB=/path/to/fftw3/lib</code></pre>
***

###Usage

To calculate channel radius using hole program. <code>hole</code> program should be already installed and present in <code>$PATH</code> environment variable.


<pre><code>gmx_hole -h
</code></pre>
***
