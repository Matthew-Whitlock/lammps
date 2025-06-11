.. index:: kill

kill command
============

Syntax
""""""

.. code-block:: LAMMPS

   kill rank kill_rank

* kill_rank = rank to kill

Examples
""""""""

.. code-block:: LAMMPS

   kill rank 1

Description
"""""""""""

Used for testing fault handling. Generates an artificial failure on kill_rank
by raising the SIGTERM signal. Calls std::abort() if SIGTERM is ignored.

Related commands
""""""""""""""""

:doc:`fix fail`
