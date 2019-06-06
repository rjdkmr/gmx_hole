.. |HOLE| raw:: html

   <a href="http://www.holeprogram.org/" target="_blank">HOLE</a>



This tool is moved to `gmx_clusterByFeatures <https://gmx-clusterbyfeatures.readthedocs.io/>`_ as sub-command `hole <https://gmx-clusterbyfeatures.readthedocs.io/en/latest/commands/hole.html>`_.
===================================================================================================================================================================================================

`gmx_clusterByFeatures <https://gmx-clusterbyfeatures.readthedocs.io/>`_ is easy to install by one command: ``sudo pip3 install gmx-clusterByFeatures``. `hole <https://gmx-clusterbyfeatures.readthedocs.io/en/latest/commands/hole.html>`_ command is same as ``gmx_hole``. **HELP** and **USAGE** for ``gmx_clusterByFeatures hole`` is provided in this `link <https://gmx-clusterbyfeatures.readthedocs.io/en/latest/commands/hole.html>`_. Moreover, `gmx_clusterByFeatures <https://gmx-clusterbyfeatures.readthedocs.io/>`_ include `holeplot <https://gmx-clusterbyfeatures.readthedocs.io/en/latest/commands/holeplot.html>`_, which can be used to plot the data obtained from `hole <https://gmx-clusterbyfeatures.readthedocs.io/en/latest/commands/hole.html>`_.

******


gmx_hole
========

It can be used to calculate radius of protein channel/cavity for GROMACS MD
trajectory. ``gmx_hole`` uses |HOLE| program to calculate radius of cavity/channel
and dumps the output to a text file as a function of tiime. It also extract
channel's outlining residues and dumps to same output file. This output file
can be further read to perform final statistcal operations and plotting.

Please cite the original publication of hole:
  O.S. Smart, J.M. Goodfellow and B.A. Wallace (1993). The Pore Dimensions of Gramicidin A. Biophysical Journal 65:2455-2460.


******


Requirements
------------

* To compile and install, **GROMACS** libraries are required.
  Presently, Gromacs **4.5.x**, **4.6.x**, **5.0.x**, **5.1.x**, **2016.x** and **2018.x**
  versions are supported.

* To use ``gmx_hole``, |HOLE| program should be already installed.


******

Download
--------

::

    git clone https://github.com/rjdkmr/gmx_hole



******

Installation
------------

.. code-block:: bash

  cd gmx_hole
  mkdir build
  cd build
  export CMAKE_PREFIX_PATH=/path/to/installed/gromacs/directory
  cmake .. -DCMAKE_INSTALL_PREFIX=/opt/gmx_hole
  make
  make install


If fftw library ``libfftw3f.so`` or ``libfftw3f.a`` are not present in standard locations:
  ::

      -DFFTW_LIB=/path/to/fftw3/lib


******

Usage
-----

To calculate channel radius using hole program. |HOLE| program should be
already installed and present in ``$PATH`` environment variable.

.. code-block:: bash

  gmx_hole -h


******
