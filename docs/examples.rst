=======================================
Examples
=======================================

.. contents::

Overview
------------
Several examples analyses are provided in the ``./examples/`` subdirectory of the `dms_tools source code`_. You can either examine these `examples via the GitHub web interface <https://github.com/jbloomlab/dms_tools/tree/master/examples>`_ or download the `dms_tools source code`_ and then look at the examples.

These examples all utilize previously published deep mutational scanning data. One of the analyses also looks at data simulated based on parameters from an actual study to confirm that `dms_tools`_ is able to accurately infer the parameters used to simulate the data. 

Melnikov et al analysis of Tn5
---------------------------------
`Melnikov et al (2014)`_ performed deep mutational scanning of all amino-acid mutations to a Tn5 transposon. In the ``./examples/Melnikov_et_al_Tn5/`` subdirectory (https://github.com/jbloomlab/dms_tools/tree/master/examples/Melnikov_et_al_Tn5), `dms_tools`_ is used to analyze and visualize the site-specific amino-acid preferences. 

Based on data from the experiment of `Melnikov et al (2014)`_, additional deep mutational scanning data is then simulated using known site-specific amino-acid preferences and differential preferences. `dms_tools`_ is then used to analyze this simulated data to confirm that (at sufficient read depth), the programs can accurately infer the known parameters used for the simulations. These simulations are in further subdirectories at https://github.com/jbloomlab/dms_tools/tree/master/examples/Melnikov_et_al_Tn5.

Thyagarajan and Bloom analysis of influenza hemagglutinin
-----------------------------------------------------------
`Thyagarajan and Bloom (2014)`_ performed deep mutational scanning on codon mutants of an H1 influenz hemagglutinin. In the  ``./examples/ThyagarajanBloom2014_WSN_HA/`` subdirectory (https://github.com/jbloomlab/dms_tools/tree/master/examples/ThyagarajanBloom2014_WSN_HA), `dms_tools`_ is used to analyze and visualize the site-specific amino-acid preferences.

Wu et al analysis of influenza NS1 in two different conditions
------------------------------------------------------------------
`Wu et al (2014)`_ created a library of nucleotide mutants of the NS1 gene of influenza hemagglutinin, and then subjected viruses carrying this library to selection in the presence and absence of type I interferon. In the ``./examples/Wu_et_al_NS1/`` subdirectory (https://github.com/jbloomlab/dms_tools/tree/master/examples/Wu_et_al_NS1), `dms_tools`_ is used to analyze and visualize the differential preferences for mutations in the two different conditions.

.. include:: weblinks.txt
