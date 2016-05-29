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

or::

    sudo pip install dms_tools

depending on whether you are `Installing locally versus globally`_.

If this command fails, then read the instructions below.

``dms_tools`` requires ``Python``
-----------------------------------

`dms_tools`_ is written in `Python`_, and is known to be compatabile with `Python`_ version 2.7 (compatibility with other `Python`_ versions has not been checked). So before installing `dms_tools`_, making sure you have `Python`_ version 2.7 installed.

Installing locally versus globally
---------------------------------------
You need to figure out whether you want to install `dms_tools`_ globally or locally. The instructions in the subsequent sections will show the options for both local and global installations -- use whichever is correct for your situation.

Install globally
~~~~~~~~~~~~~~~~~~~~

If you are the administrator (super-user) for your computer, then you can install programs globally. In that case, you will typically need to run the install commands prefaced with ``sudo``, and will be asked for your password. If you are able to install globally, then all users will be able to use `dms_tools`_.

Install locally
~~~~~~~~~~~~~~~~~~~~
If you are using a shared computing cluster, then you are probably not the administrator for the cluster. In that case, you will need to install `dms_tools`_ locally so that it can be used by you even though it will not be available to other users. The `Python convention for local installation`_ on Linux (and Mac OS X) is to place installations in the directory ``~/.local/``, typically by adding the option ``--user`` to installation commands.

In order for locally installed programs to be accessible, you need to add ``~/.local/bin/`` to the ``PATH`` variable, and ``~/.local/lib/`` to the ``PYTHONPATH`` variable. If you are using the `bash shell`_, you would do this by adding the following lines to your ``~/.bashrc`` file::

    PATH=$HOME/.local/bin/:$PATH
    export PYTHONPATH=$HOME/.local/lib/python2.7:$PATH

You then want to make sure that your ``~/.bash_profile`` file simple sources your ``~/.bashrc`` file as `described here <http://www.joshstaiger.org/archives/2005/07/bash_profile_vs.html>`_ by making ``~/.bash_profile`` consist of the following contents::

    if [ -f ~/.bashrc ]; then
        source ~/.bashrc
    fi

Recommended: installing ``dms_tools`` with ``pip``
----------------------------------------------------
The easiest way to install `dms_tools`_ is with `pip`_, which will automatically take care of installing the `Other software required by dms_tools`_.

First, make sure ``pip`` is installed
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

First, check if you already have `pip`_ installed. You can do this by typing at the command line::

    pip -h

If you get the help message for `pip`_, then `pip`_ is already installed and you can move to the next step. If you instead get an error message such as ``-bash: pip: command not found`` then you need to install `pip`_.

If you have ``easy_install`` installed, then you can simply install `pip`_ either globally with::

    sudo easy_install pip

or locally with::

    easy_install pip --user

If those commands also fail (i.e. you don't have ``easy_install`` available either), then install `pip`_ by following the `instructions here <https://pip.pypa.io/en/latest/installing.html>`_.

Next, use ``pip`` to install ``dms_tools``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Once `pip`_ is installed, you can do a local installation with::

    pip install dms_tools --user

or a global installation with::

    sudo pip install dms_tools

Upgrading your ``dms_tools`` version with ``pip``
--------------------------------------------------
If you have previously installed `dms_tools`_ but are not sure that you have the latest version, you can upgrade using `pip`_. To do this for a local installation, use::

    pip install dms_tools --user --upgrade

or globally::

    sudo pip install dms_tools --upgrade

Install ``dms_tools`` from source code
-----------------------------------------------------------------------------------
You can also install from the `dms_tools source code`_ on GitHub. If you do this, you will have to manually make sure that you have also installed the `Other software required by dms_tools`_.

To install from source, first download the `dms_tools source code`_ from GitHub.
After unpacking the main ``dms_tools`` directory that contains this source, install locally with::

    cd dms_tools
    python setup.py install --user

or globally with::

    cd dms_tools
    sudo python setup.py install 


Other software required by ``dms_tools``
------------------------------------------
`dms_tools requires Python`_. In addition it requires the following external `Python`_ packages. If you are installing with `pip`_, these external packages will automatically be installed. If you are installing from source, you will need to install these packages yourself.

The required packages are listed in the ``setup.py`` file.


.. _license:

License
-----------
`dms_tools source code`_ is available on GitHub under an open-source `GPLv3`_ license. Part of the code utilized by `dms_tools`_ is based on `weblogo`_, which is licensed under the GPL-compatible BSD 3-clause license.


.. include:: weblinks.txt
