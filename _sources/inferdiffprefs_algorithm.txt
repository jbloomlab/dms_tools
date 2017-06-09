.. _inferdiffprefs_algorithm:
  
===================================================
Algorithm to infer differential preferences
===================================================
This is a description of the algorithm that :ref:`dms_inferdiffprefs` uses to infer the differential site-specific preferences between two different selection conditions.

.. contents::
   :depth: 3

Should you use differential selection instead?
-------------------------------------------------
For many applications, we now prefer the much simpler analysis performed by :ref:`dms_diffselection`. In particular, if you are specifically looking at changes away from a defined wildtype (e.g., antibody selection), :ref:`dms_diffselection` may be better.

The definition of differential preferences.
--------------------------------------------
The script :ref:`dms_inferprefs` infers the site-specific amino-acid preferences :math:`\pi_{r,x}` for character :math:`x` at site :math:`r` using the :ref:`inferprefs_algorithm`. These site-specific preferences are inferred by comparing an unselected mutant library to variants that have passed through some selection for some function.

One could then imagine further selecting these functional variants in two different selection conditions, :math:`s1` and :math:`s2`. We might want to identify mutations that are differentially selected in the two conditions. One way to do this would simply be to infer the preferences :math:`\pi_{r,x}^{s1}` and :math:`\pi_{r,x}^{s2}` for each selection condition relative to the initial mutant library, and then compare these preferences. The disadvantage of this approach is that the comparison between the two conditions is indirect, and in general we expect the preferences for the two conditions to be mostly similar at most sites. So a more direct way is to directly compare the results for the two conditions, and estimate a parameter :math:`\Delta\pi_{r,x}` that represents the difference in preference between selection condition :math:`s1` and :math:`s2`. The advantage here is that we can directly place a prior over :math:`\Delta\pi_{r,x}` that favors the *a priori* expectation that there is no difference between the selection conditions, and then assess the posterior probability that there is a difference between the two conditions (as well as the magnitude of this difference).

We begin with some starting library, which will typically be a library of random mutants that has already gone through some initial selection for basic functionality and is now about to be further subjected to selections :math:`s1` and :math:`s2`. Let :math:`f^{\textrm{start}}_{r,x}` be the frequency of character :math:`x` at site :math:`r` in the starting library, let :math:`f^{s1}_{r,x}` be the frequency after the starting library has passed through selection :math:`s1`, and let :math:`f^{s2}_{r,x}` be the frequency after the starting library has passed through selection :math:`s2`. We take the convention of letting :math:`s1` be the "control" or "reference" selection, and looking for changes in site-specific preferences in :math:`s2` relative to :math:`s1`. Specifically, we define the differential preferences by the following equations:

.. math::
   :label: controlpi

   f_{r,x}^{s1} = \frac{f_{r,x}^{\textrm{start}} \times \pi_{r,x}^{s1}}{\sum\limits_y f_{r,y}^{\textrm{start}} \times \pi_{r,y}^{s1}}

.. math::
   :label: deltapi

   f_{r,x}^{s2} = \frac{f_{r,x}^{\textrm{start}} \times \left(\pi_{r,x}^{s1} + \Delta\pi_{r,x}\right)}{\sum\limits_y f_{r,y}^{\textrm{start}} \times \left(\pi_{r,y}^{s1} + \Delta\pi_{r,y}\right)}

where we have the constraints that

.. math::
   :label: diffprefsum

   0 = \sum\limits_x \Delta\pi_{r,x}

and

.. math::
   :label: diffprefpluspref

   0 \le \pi_{r,x}^{s1} + \Delta\pi_{r,x} \le 1.

If there is no difference between the site-specific preferences in selections :math:`s1` and :math:`s2`, then :math:`\Delta\pi_{r,x} = 0`. If :math:`x` at :math:`r` is more preferred in selection :math:`s2` than :math:`s1`, then :math:`\Delta\pi_{r,x} > 0`; conversely if :math:`x` at :math:`r` is more preferred in selection :math:`s1` than :math:`s2`, then :math:`\Delta\pi_{r,x} < 0`. 

The experimental data from deep sequencing.
-------------------------------------------
In the starting library, each site :math:`r` is potentially mutated from the wildtype identity :math:`\operatorname{wt}\left(r\right)` to any of the other possible characters. We deep sequence this starting library; let :math:`N_r^{\textrm{start}}` be the total number of reads covering :math:`r`, and let :math:`n_{r,x}^{\textrm{start}}` be the number that report :math:`x`.
We then subject the starting library to selections :math:`s1` and :math:`s2`, and sequence the resulting selected variants. Let :math:`N_r^{s1}` and :math:`N_r^{s2}` be the total number of reads covering :math:`r` in :math:`s1` and :math:`s2`, respectively; and let :math:`n_{r,x}^{s1}` and :math:`n_{r,x}^{s2}` be the number that report :math:`x`.

We allow for the possibility that the deep sequencing is not completely
accurate, and so some of the reads that report a mutant character might reflect
a sequencing error rather than a true mutation. Here we assume that the error rate
is the same for all three libraries (the starting one, the one selected in :math:`s1`, and the one selected in :math:`s2`). We quantify the error rate by control sequencing of an unmutated gene that has passed through similar experimental steps to the three libraries. 
Let :math:`N_r^{\textrm{err}}` be the total number
of sequencing reads covering :math:`r` in this error
control, and let :math:`n_{r,x}^{\textrm{err}}` be the number that report :math:`x`.
We define :math:`\xi_{r,x}` to the be the frequency at which :math:`r` is identified as :math:`x` in this error control.
As in :ref:`inferprefs_algorithm`, we assume that mutations and errors are both rare, and so will neglect as negligible the possibility that the same site experiences both a mutation and an error.

For notational compactness, we define vectors that contain the counts of
all characters at each site :math:`r` for each sample:
:math:`\mathbf{n_r^{\textbf{start}}}= \left(\cdots, \mbox{$n_{r,x}^{\textrm{start}}$}, \cdots\right)`,
:math:`\mathbf{n_r^{s1}}= \left(\cdots, \mbox{$n_{r,x}^{s1}$}, \cdots\right)`,
and
:math:`\mathbf{n_r^{\textbf{err}}}= \left(\cdots, \mbox{$n_{r,x}^{\textrm{err}}$}, \cdots\right)`.

Actual frequencies and likelihoods of observing experimental data.
------------------------------------------------------------------
The notation below follows that in :ref:`inferprefs_algorithm`, so look there for definitions of variable that are not introduced here.

We define vectors of the character frequencies, preferences, and differential preferences
at each site :math:`r`:
:math:`\boldsymbol{\mathbf{f_r^{\textbf{start}}}}= \left(\cdots, f_{r,x}^{\textrm{start}}, \cdots\right)`,
:math:`\boldsymbol{\mathbf{f_r^{s1}}}= \left(\cdots, f_{r,x}^{\textrm{s1}}, \cdots\right)`,
:math:`\boldsymbol{\mathbf{f_r^{s2}}}= \left(\cdots, f_{r,x}^{\textrm{s2}}, \cdots\right)`,
:math:`\boldsymbol{\mathbf{\xi_r}}= \left(\cdots, \xi_{r,x}, \cdots\right)`,
:math:`\boldsymbol{\mathbf{\pi^{s1}_r}}= \left(\cdots, \pi^{s1}_{r,x}, \cdots\right)`,
and
:math:`\boldsymbol{\mathbf{\Delta\pi_r}}= \left(\cdots, \Delta\pi_{r,x}, \cdots\right)`,

Equations :eq:`controlpi` and :eq:`deltapi` imply

.. math::
   :label: controlpir

   \boldsymbol{\mathbf{f_r^{s1}}}= \frac{\boldsymbol{\mathbf{f_r^{\textbf{start}}}}\circ \boldsymbol{\mathbf{\pi_r^{s1}}}}{\boldsymbol{\mathbf{f_r^{\textbf{start}}}}\cdot \boldsymbol{\mathbf{\pi_r^{s1}}}}

.. math::
   :label: deltapir

   \boldsymbol{\mathbf{f_r^{s2}}} = \frac{\boldsymbol{\mathbf{f_r^{\textbf{start}}}}\circ \left(\boldsymbol{\mathbf{\pi_r^{s1}}} + \boldsymbol{\mathbf{\Delta\pi_r}}\right)}{\boldsymbol{\mathbf{f_r^{\textbf{start}}}}\cdot \left(\boldsymbol{\mathbf{\pi_r^{s1}}} + \boldsymbol{\mathbf{\Delta\pi_r}}\right)}.

As described in :ref:`inferprefs_algorithm`, the likelihood of
observing some specific set of deep-sequencing counts will be multinomially distributed
around the actual frequencies. So

.. math::
   :label: pr_nrstart

   \Pr\left(\mathbf{n_r^{\textbf{start}}}\mid N_r^{\textrm{start}}, \boldsymbol{\mathbf{f_r^{\textbf{start}}}}, \boldsymbol{\mathbf{\xi_r}}\right) = \operatorname{Multinomial}\left(\mathbf{n_r^{\textbf{start}}}; N_r^{\textrm{start}}, \boldsymbol{\mathbf{f_r^{\textbf{start}}}} + \boldsymbol{\mathbf{\xi_r}}- \boldsymbol{\mathbf{\delta_r}}\right)

.. math::
   :label: pr_nrerr

   \Pr\left(\mathbf{n_r^{\textbf{err}}}\mid N_r^{\textrm{err}}, \boldsymbol{\mathbf{\xi_r}}\right) = \operatorname{Multinomial}\left(\mathbf{n_r^{\textbf{err}}}; N_r^{\textrm{err}}, \boldsymbol{\mathbf{\xi_r}} \right)


.. math::
   :label: pr_nrs1

   \begin{eqnarray}
   \Pr\left(\mathbf{n_r^{s1}}\mid N_r^{s1}, \boldsymbol{\mathbf{f_r^{\textbf{start}}}}, \boldsymbol{\mathbf{\pi_r^{s1}}}, \boldsymbol{\mathbf{\xi_r}}\right) &=& \operatorname{Multinomial}\left(\mathbf{n_r^{s1}}; N_r^{s1}, \boldsymbol{\mathbf{f_r^{s1}}} + \boldsymbol{\mathbf{\xi_r}}- \boldsymbol{\mathbf{\delta_r}}\right) \\
   &=& \operatorname{Multinomial}\left(\mathbf{n_r^{s1}}; N_r^{s1}, \frac{\boldsymbol{\mathbf{f_r^{\textbf{start}}}}\circ \boldsymbol{\mathbf{\pi_r^{s1}}}}{\boldsymbol{\mathbf{f_r^{\textbf{start}}}}\cdot \boldsymbol{\mathbf{\pi_r^{s1}}}} + \boldsymbol{\mathbf{\xi_r}}- \boldsymbol{\mathbf{\delta_r}}\right) \\
   \end{eqnarray}

.. math::
   :label: pr_nrs2

   \begin{eqnarray}
   \Pr\left(\mathbf{n_r^{s2}}\mid N_r^{s2}, \boldsymbol{\mathbf{f_r^{\textbf{start}}}}, \boldsymbol{\mathbf{\pi_r^{s1}}}, \boldsymbol{\mathbf{\Delta\pi_r}}, \boldsymbol{\mathbf{\xi_r}}\right) &=& \operatorname{Multinomial}\left(\mathbf{n_r^{s2}}; N_r^{s2}, \boldsymbol{\mathbf{f_r^{s2}}} + \boldsymbol{\mathbf{\xi_r}}- \boldsymbol{\mathbf{\delta_r}}\right) \\
   &=& \operatorname{Multinomial}\left(\mathbf{n_r^{s2}}; N_r^{s2}, \frac{\boldsymbol{\mathbf{f_r^{\textbf{start}}}}\circ \left(\boldsymbol{\mathbf{\Delta\pi_r}} + \boldsymbol{\mathbf{\pi_r^{s1}}}\right)}{\boldsymbol{\mathbf{f_r^{\textbf{start}}}}\cdot \left(\boldsymbol{\mathbf{\Delta\pi_r}} + \boldsymbol{\mathbf{\pi_r^{s1}}}\right)} + \boldsymbol{\mathbf{\xi_r}}- \boldsymbol{\mathbf{\delta_r}}\right) \\
   \end{eqnarray}


Priors over the unknown parameters.
-----------------------------------

We specify `Dirichlet`_ priors over the four parameter vectors for each
site :math:`r`:

.. math::
   :label: pr_pir

   \Pr\left(\boldsymbol{\mathbf{\pi_r^{s1}}}\right) = \operatorname{Dirichlet}\left(\boldsymbol{\mathbf{\pi_r^{s1}}}; \alpha_{\pi^{s1}} \times \mathbf{1}\right)

.. math::
   :label: pr_xir

   \Pr\left(\boldsymbol{\mathbf{\xi_r}}\right) = \operatorname{Dirichlet}\left( \boldsymbol{\mathbf{\xi_r}}; \alpha_{\xi} \times \mathcal{N}_x \times \boldsymbol{\mathbf{a_{r,\xi}}}\right)

.. math::
   :label: pr_frstart

   \Pr\left(\boldsymbol{\mathbf{f_r^{\textbf{start}}}}\right) = \operatorname{Dirichlet}\left( \boldsymbol{\mathbf{f_r^{\textbf{start}}}}; \alpha_{\textrm{start}} \times \mathcal{N}_x \times \boldsymbol{\mathbf{a_{r,\textbf{start}}}}\right) 

.. math::
   :label: pr_deltapi

   \Pr\left(\boldsymbol{\mathbf{\Delta\pi_r}} \mid \boldsymbol{\mathbf{\pi_r^{s1}}}\right) = \operatorname{Dirichlet}\left(\boldsymbol{\mathbf{\Delta\pi_r}}; \alpha_{\Delta\pi} \times \mathcal{N}_x \times \boldsymbol{\mathbf{\pi_r^{s1}}}\right) - \boldsymbol{\mathbf{\pi_r^{s1}}}

where the :math:`\alpha` parameters 
are scalar concentration parameters with values
:math:`> 0` specified by the user, and the :math:`\mathbf{a_r}` vectors
have entries :math:`> 0` that sum to one. You might consider choosing a value slightly greater than one for :math:`\alpha_{\Delta\pi}` to enforce the prior expectation that :math:`\Delta\pi_{r,x}` is close to zero.

We specify the prior vectors :math:`\mathbf{a_{r,\textbf{start}}}` and :math:`\mathbf{a_{r,\xi}}` in terms of the average per-site mutation or error rates over the entire library. 

Our prior assumption is that the rate of sequencing errors depends on the number of nucleotides being changed -- for nucleotide characters all mutations have only one nucleotide changed, but for codon characters there can be one, two, or three nucleotides changed. 
Specifically, the average per-site error rate for mutations with :math:`m` nucleotide changes :math:`\overline{\xi_m}` is

.. math::
   :label: avgxim

   \overline{\xi_m} = \frac{1}{L}\sum\limits_r \frac{1}{N_r^{\textrm{err}}}\sum\limits_{x} n_{r,x}^{\textrm{err}}\times \delta_{m,D_{x,\operatorname{wt}\left(r\right)}} 

where :math:`r` ranges over all :math:`L` sites with deep mutational scanning data, :math:`x` ranges over all characters (nucleotides or codons), and :math:`D_{x,\operatorname{wt}\left(r\right)}` is the number of nucleotide differences between :math:`x` and the wildtype character :math:`\operatorname{wt}\left(r\right)`. 
We define the prior vector as

.. math::
   :label: arxi

   \boldsymbol{\mathbf{a_{r,\xi}}} = \left(\cdots, \sum\limits_m \frac{\overline{\xi}_{m}}{\mathcal{C}_m} \times \delta_{m,D_{x,\operatorname{wt}\left(r\right)}},\cdots\right),

where :math:`\mathcal{C}_m` is the number of mutant characters with :math:`m` changes relative to the wildtype.

Our prior assumption is that all mutations are at equal frequency in the starting library (although this assumption is unlikely to be true if the starting library has already been subjected to some functional selection, we lack a rationale for another prior in the absence of the deep mutational scanning). In estimating the average frequency of these mutations, we need to account for the fact that the observed counts of non-wildtype characters in the starting library will be inflated by the error rate represented by :math:`\overline{\xi_m}`. We therefore estimate the average mutation frequency

.. math::
   :label: avgfstart

   \overline{f^{\textrm{start}}} = \left(\frac{1}{L}\sum\limits_r \frac{1}{N_r^{\textrm{start}}}\sum\limits_{x\ne \operatorname{wt}\left(r\right)} n_{r,x}^{\textrm{start}}\right) - \sum\limits_{m \ge 1} \overline{\xi_m}.

Note that a value of :math:`\overline{f^{\textrm{start}}} \le 0` would suggest that mutations are no more prevalent in the starting library then the error control. Such a situation would violate the assumptions of the experiment, and so the algorithm will halt if this is the case. The prior vector for the mutagenesis rate is 

.. math::
   :label: arstart

   \boldsymbol{\mathbf{a_{r,\textbf{start}}}} = \left(\cdots, \frac{\overline{f^{\textrm{start}}}}{\mathcal{N}_x - 1} + \delta_{x,\operatorname{wt}\left(r\right)} \times \left[1 - \overline{f^{\textrm{start}}}\right] ,\cdots\right).

These priors and likelihoods can be applied to several types of characters as described in :ref:`inferprefs_algorithm`. One possibility is that :ref:`chartype_DNA` A second possibility is :ref:`chartype_codon` A third possibility is :ref:`chartype_codon_to_aa`. A fourth possibility is :ref:`chartype_aa`.

The actual preferences are inferred by using `MCMC`_ to sample the product of the priors and likelihoods described above. This `MCMC` is done using the same computational approach and convergence criteria described in :ref:`MCMC_inference` But rather than reporting credible intervals, we report the posterior probability that :math:`\Delta\pi_{r,x} > 0`; this gives us an estimate of the probability that :math:`x` is favored at :math:`r` in selection condition :math:`s2` versus :math:`s1`.

.. include:: weblinks.txt
