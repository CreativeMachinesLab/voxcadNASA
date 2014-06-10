voxelyze-14-07-09
==================

Created: Jun 9, 2014



Description:

This directory contains a copy of JonH's voxelyze library,
version 0.9.92.  It has some edits to get it to compile on
a Linux system and a Makefile has been added to compile it
into a library which can be linked to.

The library can be created by going into the 'Voxelyze'
library and running: $ make



Installation:

Compiling the library is done with:
  $ make



This library and needed include files can be installed either into
your user account ($(HOME)) or into the system account.  To install
into your user account this assumes that the following directories
exist:
o $(HOME)/include
o $(HOME)/lib



These installs will put copies of files in the following directories:

installusr:
  headers: ~/include/voxelyze.X.Y
  library: ~/lib/voxelyze.X.Y.a

installglobal:
  headers: /usr/local/include/voxelyze.X.Y
  library: /usr/local/lib/voxelyze.X.Y.a


In addition, links will be created to these as simply:
   voxelyze.X or voxelyze.X.a
