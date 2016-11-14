===========================
Documentation
===========================

This subdirectory contains the `reStructuredText`_ documentation that can be built with `Sphinx`_.

To build the documentation, you will need both `Sphinx`_ and the `sphinx-argparse`_ to be installed. If you don't have these, install them with ``pip``, such as by ``pip install sphinx-argparse``.

Then simply make the documentation with::

    make html

and the HTML documentation will be installed in ``./_build/html/``.

If you are getting an error, it is probably because you forget to install `sphinx-argparse`_.

.. _`reStructuredText`: http://docutils.sourceforge.net/docs/user/rst/quickref.html
.. _`Sphinx`: http://sphinx-doc.org/
.. _`sphinx-argparse`: http://sphinx-argparse.readthedocs.org
