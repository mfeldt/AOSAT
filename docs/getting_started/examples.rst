====================
Running the Examples
====================

Example files
=============

When you have pulled AOSAT from GitHub, an '''examples''' directory has been
checked out along with the source.

Following the installation instructions actually installs the examples in the
destination directory.

To find the examples, you can do the following in python (running of course, in
your activated AOSAT environment):

.. code-block::

  >>> import os
  >>> os.path.join(os.path.split(aosat.__file__)[0],'examples')


A full "tearsheet"
==================

A common use case is to produce a so-called "tearsheet", a double-sided page
summarizing the most interesting properties of a given AO simulation.

This can be achieved easily for the provided close-loop example like so:

.. code-block::

  >>> import aosat
  >>> example_path = os.path.join(os.path.split(aosat.__file__)[0],'examples')
  >>> example_file = os.path.join(example_path,'example_analyze_closed_loop.setup')
  >>> aosat.analyze.tearsheet(example_file)

This will produce two files in the current working directory:

.. code-block::

  ts_test.txt
  ts_test.pdf

Guess what the extensions mean and look at them with the appropriate tool for each
to see what they are about.

Running a specific analyzer
===========================

Sometimes (well, actually many times) one is not interested in the full tearsheet
thingy, and may just want to perform a particular analysis.  Let's say you only
want to know about how your simulation behaved with respect to pupil
fragmentation. In this case, for the same example data set, you would do:

.. code-block::

  >>> import aosat
  >>> from aosat.analyze import frg_analyzer, run, setup
  >>> import os
  >>> example_path = os.path.join(os.path.split(aosat.__file__)[0],'examples')
  >>> example_file = os.path.join(example_path,'example_analyze_closed_loop.setup')
  >>> aosat_cfg.CFG_SETTINGS = aosat_cfg.configure(example_file)
  >>> sd=setup()
  >>> analyzers=[frg_analyzer(sd)]
  >>> run(analyzers)

This will run only the fragmentation analyzer on th example set.
To plot the result, you could do

.. code-block::

  >>> fig = analyzers[0].make_plot()
  >>> fig.savefig('fragmentation.pdf')

Or, you could look at individual results by e.g doing

.. code-block::

  >>> import matplotlib.pyplot as plt
  >>> plt.plot(analyzers[0].pistont[:,0])

To get the behaviour of the piston of the 0th pupil fragment. What kind of data
is made available by individual analyzers can be looked up in the respective
documentation, in this case of the :doc:`Pupil Fragmentation Analyzer<../analyzers/frg_analyzer>`.
