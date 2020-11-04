====================
Community Guidelines
====================

Getting Help
============

If you do not understand how something is supposed to work, how to set up a particular analysis, or how a particular result comes about, please try to find guidance first in the relevant sections of this documentation. If all else fails, contact the authors by email. (see https://github.com/mfeldt/AOSAT/blob/master/AUTHORS.rst)


Issues
======

If AOSAT does not perform as expected, delivers clearly wrong results, or plainly crashes, you can open an issue over on github: https://github.com/mfeldt/AOSAT


Contributing
============

AOSAT is meant to be easily extensible by e.g. integrating in simulation environments, or by contributing new analyzers.
If you intend to do so without having an idea where to start, you can

* look at an existing piece of code like e.g. an analyzer to derive your own from there
* read the relevant sections of this documentation, i.e. :doc:`Frameserver <frameserver>` and :doc:`Analyzers <analyzers>`

Adding your contribution
------------------------

If you make new analyzer that you deem useful to the community at large (we're e.g. still struggling to make a nice temporal power spectrum...), it's worth considering to add it to AOSAT.

To so, consider the following:

AOSAT follows a development process close to the gitflow model.  To contribute, clone the dev branch like so:

.. code::

  git clone -b dev https://github.com/mfeldt/AOSAT.git

Then, create a new branch and give it a useful name identifying your contributed feature:

.. code::

  git checkout -b feature_nameOfFeature

Then add, commit, and push changes.  If you think your feature is ready, send the authors an email and we'll integrate it into the dev branch, and into the main branch at the next release.
