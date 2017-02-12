.. _installation:

================================
Installation
================================

`dms_tools`_ is freely available. Here are instructions on how to install the package.

.. contents::

Quick installation
---------------------
If your system already has the appropriate required software, you can install `dms_tools`_ with the simple command::

    pip install dms_tools --user

If this command fails, then read the instructions below.

``dms_tools`` requires ``Python``
-----------------------------------

`dms_tools`_ is written in `Python`_, and is compatible with `Python`_ version 2.7. So before installing `dms_tools`_, making sure you have `Python`_ version 2.7 installed.

Where to install
---------------------------------------
You need to figure out where you want to install `dms_tools`_.
Global installation using ``sudo`` `is not recommended for Python packages in general <http://stackoverflow.com/questions/21055859/what-are-the-risks-of-running-sudo-pip/21056000#21056000>`_.

The simplest solution is to install locally via the ``--user`` option to ``pip``, which by default on Linux will install into the ``~/.local/`` directory.

In order for locally installed programs to be accessible, you need to add ``~/.local/bin/`` to the ``PATH`` variable, and ``~/.local/lib/`` to the ``PYTHONPATH`` variable. If you are using the `bash shell`_, you would do this by adding the following lines to your ``~/.bashrc`` file::

    PATH=$HOME/.local/bin/:$PATH
    export PYTHONPATH=$HOME/.local/lib/python2.7:$PATH

You then want to make sure that your ``~/.bash_profile`` file simple sources your ``~/.bashrc`` file as `described here <http://www.joshstaiger.org/archives/2005/07/bash_profile_vs.html>`_ by making ``~/.bash_profile`` consist of the following contents::

    if [ -f ~/.bashrc ]; then
        source ~/.bashrc
    fi

On Mac OS X, the default directory for ``--user`` may be ``$HOME/Library/Python/x.y/`` rather than ``~/.local/`` where ``x.y`` indicates the version number (e.g., ``2.7``.

Installing with ``pip`` and ``--user``
----------------------------------------------------
Once you have taken care of the steps above, you can then simply install with ``pip``.

First, make sure ``pip`` is installed
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

First, check if you already have `pip`_ installed. You can do this by typing at the command line::

    pip -h

If you get the help message for `pip`_, then `pip`_ is already installed and you can move to the next step. If you instead get an error message such as ``-bash: pip: command not found`` then you need to install `pip`_.

If you have ``easy_install`` installed, then you can simply install `pip`_ with::

    easy_install pip --user

If this fails (e.g., you don't have ``easy_install`` available either), then install `pip`_ by following the `instructions here <https://pip.pypa.io/en/latest/installing.html>`_.

Next, use ``pip`` to install ``dms_tools``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Once `pip`_ is installed, you can do a local installation with::

    pip install dms_tools --user

Using a virtual environment
-----------------------------
The other good option is to use ``pip`` to install into a virtual environment `as described here <http://docs.python-guide.org/en/latest/dev/virtualenvs/>`_

Upgrading your ``dms_tools`` version with ``pip``
--------------------------------------------------
If you have previously installed `dms_tools`_ but are not sure that you have the latest version, you can upgrade using `pip`_. To do this for a local installation, use::

    pip install dms_tools --user --upgrade


Install ``dms_tools`` from source code
-----------------------------------------------------------------------
You can also install from the `dms_tools source code`_ on GitHub. If you do this, you may have to manually make sure that you have also installed the `Other required software`_.

To install from source, first download the `dms_tools source code`_ from GitHub.
After unpacking the main ``dms_tools`` directory that contains this source, install locally with::

    cd dms_tools
    python setup.py install --user


Other required software
------------------------------------------
`dms_tools requires Python`_. In addition it requires the following external `Python`_ packages. If you are installing with `pip`_, these external packages will automatically be installed. If you are installing from source, you may need to install these packages yourself.

The required packages are listed in the ``setup.py`` file.


.. _license:

License
-----------
`dms_tools source code`_ is available on GitHub under an open-source `GPLv3`_ license. Part of the code utilized by `dms_tools`_ is based on `weblogo`_, which is licensed under the GPL-compatible BSD 3-clause license.


.. include:: weblinks.txt
