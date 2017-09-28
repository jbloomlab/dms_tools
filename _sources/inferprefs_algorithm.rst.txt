.. _inferprefs_algorithm:
  
===================================================
Algorithm to infer site-specific preferences
===================================================
Here is a description of the algorithm that :ref:`dms_inferprefs` uses to infer the site-specific preferences.

.. contents::
   :depth: 3

Definition of the site-specific preferences.
--------------------------------------------

We assume each site :math:`r` in our gene has a preference
:math:`\pi_{r,x} \ge 0` for each character :math:`x`, where the
characters can be nucleotides, codons, or amino acids. The preferences
are defined as follows: if :math:`\mu_{r,x}` is the frequency of
:math:`x` at :math:`r` in the mutant library pre-selection and
:math:`f_{r,x}` is the frequency of :math:`x` in the library
post-selection, then

.. math::
   :label: pirx

   f_{r,x} = \frac{\mu_{r,x} \times \pi_{r,x}}{\sum_y \mu_{r,y} \times \pi_{r,y}}.

This equation implies that :math:`1 = \sum_x \pi_{r,x}`.

The experimental data from deep sequencing.
-------------------------------------------

We create mutant libraries such that each site :math:`r` is potentially
mutated from the wildtype identity  to any of the other possible
characters. We use deep sequencing to count the appearances of each
character :math:`x` at each site :math:`r` in this mutant library; since
this sequencing is performed prior to selection for gene function, we
refer to it as *pre*-selection. Let :math:`N_r^{\textrm{pre}}` be the total number of sequencing
reads covering :math:`r`, and let :math:`n_{r,x}^{\textrm{pre}}`  be the number that report :math:`x`
(note that
:math:`\mbox{$N_r^{\textrm{pre}}$}= \sum_x \mbox{$n_{r,x}^{\textrm{pre}}$}`).

We impose an experimental selection on the mutant library with some
biologically relevant pressure that favors some mutations and disfavor
others. We then use deep sequencing to count the characters in this
selected library; since this sequencing is performed after selection, we
refer to it as *post*-selection. Let :math:`N_r^{\textrm{post}}` be the total number of sequencing
reads covering :math:`r`, and let :math:`n_{r,x}^{\textrm{post}}` be the number that report :math:`x`
(note that
:math:`\mbox{$N_r^{\textrm{post}}$}= \sum_x \mbox{$n_{r,x}^{\textrm{post}}$}`).

We allow for the possibility that the deep sequencing is not completely
accurate; for instance, perhaps some of the reads that report a mutant
character reflect a sequencing error rather than a true mutation. The
rates of such errors can be quantified by control sequencing of the
unmutated gene, so that any counts at :math:`r` of
:math:`x \ne \mbox{$\operatorname{wt}\left(r\right)$}` reflect sequencing errors. It
is possible that the error rates for sequencing the pre-selection and
post-selection libraries are different, as for instance would arise in
sequencing an RNA virus for which the post-selection libraries must be
reverse-transcribed to DNA prior to sequencing. Let :math:`N_r^{\textrm{err,pre}}` be the total number
of sequencing reads covering :math:`r` in the pre-selection error
control, and let :math:`n_{r,x}^{\textrm{err,pre}}` be the number that report :math:`x` (note that
:math:`\mbox{$N_r^{\textrm{err,pre}}$}= \sum_x \mbox{$n_{r,x}^{\textrm{err,pre}}$}`).
Make similar definitions of :math:`N_r^{\textrm{err,post}}` and :math:`n_{r,x}^{\textrm{err,post}}` for the post-selection error control.

For notational compactness, we define vectors that contain the counts of
all characters at each site :math:`r` for each sample:
:math:`\mbox{$\mathbf{n_r^{\textbf{pre}}}$}= \left(\cdots, \mbox{$n_{r,x}^{\textrm{pre}}$}, \cdots\right)`,
:math:`\mbox{$\mathbf{n_r^{\textbf{post}}}$}= \left(\cdots, \mbox{$n_{r,x}^{\textrm{post}}$}, \cdots\right)`,
:math:`\mbox{$\mathbf{n_r^{\textbf{err,pre}}}$}= \left(\cdots, \mbox{$n_{r,x}^{\textrm{err,pre}}$}, \cdots\right)`,
and
:math:`\mbox{$\mathbf{n_r^{\textbf{err,post}}}$}= \left(\cdots, \mbox{$n_{r,x}^{\textrm{err,post}}$}, \cdots\right)`.

Assumption that mutations and errors are rare.
----------------------------------------------

The samples described above allow for the possibility of errors
as well as the actual mutations.
We assume that the the mutagenesis and error rates rate at each site are
sufficiently low that most characters are wildtype, so for instance
:math:`\mbox{$N_r^{\textrm{pre}}$}\sim \mbox{$n_{r,\operatorname{wt}\left(r\right)}^{\textrm{pre}}$}\gg \mbox{$n_{r,x}^{\textrm{pre}}$}`
for all :math:`x \ne \mbox{$\operatorname{wt}\left(r\right)$}`, and
:math:`\mbox{$N_r^{\textrm{err,pre}}$}\sim \mbox{$n_{r,\operatorname{wt}\left(r\right)}^{\textrm{err,pre}}$}\gg \mbox{$n_{r,x}^{\textrm{err,pre}}$}`
for all :math:`x \ne \mbox{$\operatorname{wt}\left(r\right)$}`. This allows us to
ignore as negligibly rare the possibility that the same site experiences
both a mutation and an error in a single molecule, which simplifies the
analysis below.

Actual frequencies and likelihoods of observing experimental data.
------------------------------------------------------------------

We have defined :math:`\mu_{r,x}` and :math:`f_{r,x}` as the frequencies
of :math:`x` at site :math:`r` in the pre-selection and post-selection
libraries, respectively. Also define :math:`\epsilon_{r,x}` to be the
frequency at which :math:`r` is identified as :math:`x` in sequencing
the error control for the pre-selection library, and let
:math:`\rho_{r,x}` be the frequency at which :math:`r` is identified as
:math:`x` in sequencing the error control for the post-selection
library. For notational compactness, we define vectors of these
frequencies for all characters at each site :math:`r`:
:math:`\mbox{$\boldsymbol{\mathbf{\mu_r}}$}= \left(\cdots, \mu_{r,x}, \cdots\right)`,
:math:`\mbox{$\boldsymbol{\mathbf{f_r}}$}= \left(\cdots, f_{r,x}, \cdots\right)`,
:math:`\mbox{$\boldsymbol{\mathbf{\epsilon_r}}$}= \left(\cdots, \epsilon_{r,x}, \cdots\right)`,
and
:math:`\mbox{$\boldsymbol{\mathbf{\rho_r}}$}= \left(\cdots, \rho_{r,x}, \cdots\right)`.
We also define
:math:`\mbox{$\boldsymbol{\mathbf{\pi_r}}$}= \left(\cdots, \pi_{r,x}, \cdots\right)`.
Note that Equation :eq:`pirx` implies that

.. math::
   :label: pir

   \mbox{$\boldsymbol{\mathbf{f_r}}$}= \frac{\mbox{$\boldsymbol{\mathbf{\mu_r}}$}\circ \mbox{$\boldsymbol{\mathbf{\pi_r}}$}}{\mbox{$\boldsymbol{\mathbf{\mu_r}}$}\cdot \mbox{$\boldsymbol{\mathbf{\pi_r}}$}}

where :math:`\circ` is the `Hadamard product`_.

Unless we exhaustively sequence every molecule in each library, the
observed counts will not precisely reflect their actual frequencies in
the libraries, since the deep sequencing only observes a sample of the
molecules. If the deep sequencing represents a random sampling of a
small fraction of a large library of molecules, then the likelihood of
observing some specific set of counts will be multinomially distributed
around the actual frequencies. So

.. math::
   :label: pr_nrpre

   \Pr\left(\mbox{$\mathbf{n_r^{\textbf{pre}}}$}\mid \mbox{$N_r^{\textrm{pre}}$}, \mbox{$\boldsymbol{\mathbf{\mu_r}}$}, \mbox{$\boldsymbol{\mathbf{\epsilon_r}}$}\right) = \operatorname{Multinomial}\left(\mbox{$\mathbf{n_r^{\textbf{pre}}}$}; \mbox{$N_r^{\textrm{pre}}$}, \mbox{$\boldsymbol{\mathbf{\mu_r}}$}+ \mbox{$\boldsymbol{\mathbf{\epsilon_r}}$}- \mbox{$\boldsymbol{\mathbf{\delta_r}}$}\right)

where :math:`\operatorname{Multinomial}` denotes the `Multinomial`_ distribution,
:math:`\mbox{$\boldsymbol{\mathbf{\delta_r}}$}= \left(\cdots, \delta_{x, \operatorname{wt}\left(r\right)}, \cdots\right)`
is a vector for which the element corresponding to :math:`\operatorname{wt}\left(r\right)` is one and all other
elements are zero (:math:`\delta_{xy}` is the `Kronecker delta`_), and we
have assumed that the probability that a site experiences both a
mutation and an error is negligibly small. Similarly,

.. math::
   :label: pr_nrerrpre

   \Pr\left(\mbox{$\mathbf{n_r^{\textbf{err,pre}}}$}\mid \mbox{$N_r^{\textrm{err,pre}}$}, \mbox{$\boldsymbol{\mathbf{\epsilon_r}}$}\right) = \operatorname{Multinomial}\left(\mbox{$\mathbf{n_r^{\textbf{err,pre}}}$}; \mbox{$N_r^{\textrm{err,pre}}$}, \mbox{$\boldsymbol{\mathbf{\epsilon_r}}$}\right)

.. math::
   :label: pr_nrerrpost

   \Pr\left(\mbox{$\mathbf{n_r^{\textbf{err,post}}}$}\mid \mbox{$N_r^{\textrm{err,post}}$}, \mbox{$\boldsymbol{\mathbf{\rho_r}}$}\right) = \operatorname{Multinomial}\left(\mbox{$\mathbf{n_r^{\textbf{err,post}}}$}; \mbox{$N_r^{\textrm{err,post}}$}, \mbox{$\boldsymbol{\mathbf{\rho_r}}$}\right)

.. math::
   :label: pr_nrpost

   \begin{eqnarray}
   \Pr\left(\mbox{$\mathbf{n_r^{\textbf{post}}}$}\mid \mbox{$N_r^{\textrm{post}}$}, \mbox{$\boldsymbol{\mathbf{\mu_r}}$}, \mbox{$\boldsymbol{\mathbf{\pi_r}}$}, \mbox{$\boldsymbol{\mathbf{\rho_r}}$}\right) &=& \operatorname{Multinomial}\left(\mbox{$\mathbf{n_r^{\textbf{post}}}$}; \mbox{$N_r^{\textrm{post}}$}, \mbox{$\boldsymbol{\mathbf{f_r}}$}+ \mbox{$\boldsymbol{\mathbf{\rho_r}}$}- \mbox{$\boldsymbol{\mathbf{\delta_r}}$}\right) \\
   &=& \operatorname{Multinomial}\left(\mbox{$\mathbf{n_r^{\textbf{post}}}$}; \mbox{$N_r^{\textrm{post}}$}, \frac{\mbox{$\boldsymbol{\mathbf{\mu_r}}$}\circ \mbox{$\boldsymbol{\mathbf{\pi_r}}$}}{\mbox{$\boldsymbol{\mathbf{\mu_r}}$}\cdot \mbox{$\boldsymbol{\mathbf{\pi_r}}$}} + \mbox{$\boldsymbol{\mathbf{\rho_r}}$}- \mbox{$\boldsymbol{\mathbf{\delta_r}}$}\right).
   \end{eqnarray}

Priors over the unknown parameters.
-----------------------------------

We specify `Dirichlet`_ priors over the four parameter vectors for each
site :math:`r`:

.. math::
   :label: pr_pir

   \Pr\left(\boldsymbol{\mathbf{\pi_r}}\right) = \operatorname{Dirichlet}\left(\boldsymbol{\mathbf{\pi_r}}; \alpha_{\pi} \times \mathbf{1}\right)

.. math::
   :label: pr_mur

   \Pr\left(\boldsymbol{\mathbf{\mu_r}}\right) = \operatorname{Dirichlet}\left(\boldsymbol{\mathbf{\mu_r}}; \alpha_{\mu} \times \mathcal{N}_x \times \boldsymbol{\mathbf{a_{r,\mu}}}\right)

.. math::
   :label: pr_epsilonr

   \Pr\left(\boldsymbol{\mathbf{\epsilon_r}}\right) = \operatorname{Dirichlet}\left(\boldsymbol{\mathbf{\epsilon_r}}; \alpha_{\epsilon} \times \mathcal{N}_x \times \boldsymbol{\mathbf{a_{r,\epsilon}}}\right) 

.. math::
   :label: pr_rhor

   \Pr\left(\boldsymbol{\mathbf{\rho_r}}\right) = \operatorname{Dirichlet}\left(\boldsymbol{\mathbf{\rho_r}}; \alpha_{\rho} \times \mathcal{N}_x \times \boldsymbol{\mathbf{a_{r,\rho}}}\right)

where :math:`\mathbf{1}` is a vector of ones, :math:`\mathcal{N}_x` is
the number of characters (i.e. 64 for codons, 20 for amino acids, 4 for
nucleotides), the :math:`\alpha` parameters (i.e. :math:`\alpha_{\pi}`,
:math:`\alpha_{\mu}`, :math:`\alpha_{\epsilon}`, and
:math:`\alpha_{\rho}`) are scalar concentration parameters with values
:math:`> 0` specified by the user, and the :math:`\mathbf{a_r}` vectors
have entries :math:`> 0` that sum to one.

We specify the prior vectors :math:`\mathbf{a_{r,\mu}}`, :math:`\mathbf{a_{r,\epsilon}}`, and :math:`\mathbf{a_{r,\rho}}` in terms of the average per-site mutation or error rates over the entire library. 

Our prior assumption is that the rate of sequencing errors depends on the number of nucleotides being changed -- for nucleotide characters all mutations have only one nucleotide changed, but for codon characters there can be one, two, or three nucleotides changed. 
Specifically, the average per-site error rate for mutations with :math:`m` nucleotide changes :math:`\overline{\epsilon_m}` and :math:`\overline{\rho_m}` in the pre-selection and post-selection controls, respectively, are defined as

.. math::
   :label: avgepsilonm

   \overline{\epsilon_m} = \frac{1}{L}\sum\limits_r \frac{1}{\mbox{$N_r^{\textrm{err,pre}}$}}\sum\limits_{x} \mbox{$n_{r,x}^{\textrm{err,pre}}$}\times \delta_{m,D_{x,\operatorname{wt}\left(r\right)}} 

.. math::
   :label: avgrhom

   \overline{\rho_m} = \frac{1}{L}\sum\limits_r \frac{1}{\mbox{$N_r^{\textrm{err,post}}$}}\sum\limits_{x} \mbox{$n_{r,x}^{\textrm{err,post}}$}\times \delta_{m,D_{x,\operatorname{wt}\left(r\right)}}

where :math:`L` is the number of sites with deep mutational scanning data, :math:`r` ranges over all of these sites, :math:`x` ranges over all characters (nucleotides or codons), and :math:`D_{x,\operatorname{wt}\left(r\right)}` is the number of nucleotide differences between :math:`x` and the wildtype character :math:`\operatorname{wt}\left(r\right)`. Note that :math:`1 = \sum\limits_m \overline{\epsilon_m} = \sum\limits_m \overline{\rho_m}`. 

Given these definitions, we define the prior vectors for the error rates as

.. math::
   :label: arepsilon

   \mbox{$\boldsymbol{\mathbf{a_{r,\epsilon}}}$} = \left(\cdots, \sum\limits_m \frac{\overline{\epsilon}_{m}}{\mathcal{C}_m} \times \delta_{m,D_{x,\operatorname{wt}\left(r\right)}},\cdots\right),

.. math::
   :label: arrho

   \mbox{$\boldsymbol{\mathbf{a_{r,\rho}}}$} = \left(\cdots, \sum\limits_m \frac{\overline{\rho}_{m}}{\mathcal{C}_m} \times \delta_{m,D_{x,\operatorname{wt}\left(r\right)}}, \cdots\right)

where :math:`\mathcal{C}_m` is the number of mutant characters with :math:`m` changes relative to the wildtype (so for nucleotides :math:`\mathcal{C}_0 = 1` and :math:`\mathcal{C}_1 = 3`, while for codons :math:`\mathcal{C}_0 = 1`, :math:`\mathcal{C}_1 = 9`, :math:`\mathcal{C}_2 = \mathcal{C}_3 = 27`) and :math:`\delta_{xy}` is again the `Kronecker delta`_.

Our prior assumption is that the mutagenesis is done so that all mutant characters are introduced at equal frequency -- note that this assumption is only true for codon characters if the mutagenesis is done at the codon level. In estimating the mutagenesis rate, we need to account for the facts that the observed counts of non-wildtype characters in the pre-selection mutant library will be inflated by the error rate represented by :math:`\overline{\epsilon_m}`. We therefore estimate the per-site mutagenesis rate as

.. math::
   :label: avgmu

   \overline{\mu} = \left(\frac{1}{L}\sum\limits_r \frac{1}{\mbox{$N_r^{\textrm{pre}}$}}\sum\limits_{x\ne \operatorname{wt}\left(r\right)} \mbox{$n_{r,x}^{\textrm{pre}}$}\right) - \sum\limits_{m \ge 1} \overline{\epsilon_m}.

Note that a value of :math:`\overline{\mu} \le 0` would suggest that mutations are no more prevalent (or actually less prevalent) in the pre-selection mutant library then the pre-selection error control. Such a situation would violate the assumptions of the experiment, and so the algorithm will halt if this is the case. The prior vector for the mutagenesis rate is then 

.. math::
   :label: armu

   \mbox{$\boldsymbol{\mathbf{a_{r,\mu}}}$} = \left(\cdots, \frac{\overline{\mu}}{\mathcal{N}_x - 1} + \delta_{x,\operatorname{wt}\left(r\right)} \times \left[1 - \overline{\mu}\right] ,\cdots\right).


.. _chartype_DNA:

Characters are nucleotides; preferences are for nucleotides.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
One scenario is that the deep sequencing counts are for nucleotide
characters :math:`x`, and that we want to determine the preference
:math:`\pi_{r,x}` for each nucleotide at each site :math:`r`. This is
the most natural approach if the mutant library is generated by a
nucleotide-level mutagenesis process, such as error-prone PCR. Note that
in this case the priors correspond to the assumption that all single-nucleotide
mutations and errors are equally likely. In this scenario, we use the prior
vectors specified by Equations :eq:`arepsilon`, :eq:`arrho`, and :eq:`armu`,
with the summations over :math:`x` covering the four nucleotide identities
(``A``, ``C``, ``G``, and ``T``) and the summations over :math:`m` covering
the two possible number of nucleotide changes to a nucleotide character (0 and 1).

.. _chartype_codon:

Characters are codons; preferences are for codons.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
A second
scenario is that the deep sequencing counts are for codon
characters :math:`x`, and that we want to determine the preference
:math:`\pi_{r,x}` for each codon at each site :math:`r`. This is a
natural approach if the mutant library is generated by introducing
all possible codon mutations at each site, as is done by techniques
like the one described by `Firnberg2012`_. The priors above based on the assumption
that all possible codon mutations are made at equal frequency, so
these priors are only appropriate if codons are mutagenized to
``NNN`` (they would need to be adjusted if codons are only mutated to
a subset of possibilities such as ``NNK``). In this scenario, we use the prior
vectors specified by Equations :eq:`arepsilon`, :eq:`arrho`, and :eq:`armu`,
with the summations over :math:`x` covering the 64 codon identities
(``AAA``, ``AAC``, ``AAG``, ``AAT``, ``ACA``, etc) and the summations over :math:`m` covering
the four possible number of nucleotide changes to a codon character (0, 1, 2, and 3).

.. _chartype_codon_to_aa:

Characters are codons; preferences are for amino acids.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
A third
possibility is that the deep sequencing counts are for codon
characters :math:`x`, but that we want to determine the preferences
:math:`\pi_{r,a}` for amino-acid characters :math:`a`. This would be
the appropriate approach if the mutant library is generated by
introducing all possible codon mutations at each site and we are
assuming that all selection acts at the protein level (so all codons
for the same amino-acid should have the same preference). In this
case, we define the vector-valued function :math:`\mathbf{A}` to map
from codons to amino acids, so that

.. math::
   :label: A

   \mathbf{A}\left(\mathbf{w}\right) = \left(\cdots, \sum\limits_x \delta_{a,\mathcal{A}\left(x\right)} \times w_x, \cdots\right)

where :math:`\mathbf{w}` is a 64-element vector giving the values for
each codon :math:`x`, :math:`\mathbf{A}\left(\mathbf{w}\right)` is a
20-element vector giving the values for each amino acid :math:`a` (or
a 21-element vector if stop codons are considered a possible
amino-acid identity), and :math:`\mathcal{A}\left(x\right)` is the
amino acid encoded by codon :math:`x`. The parameter vectors :math:`\boldsymbol{\mathbf{\pi_r}}`, :math:`\boldsymbol{\mathbf{\mu_r}}`, :math:`\boldsymbol{\mathbf{\epsilon_r}}`,
and :math:`\boldsymbol{\mathbf{\rho_r}}` are then all of length 20 (or 21) for the amino acid vectors.
The first of these vectors is still a symmetric `Dirichlet`_, and the
priors for the remaining three are
:math:`\mathbf{A}\left(\mbox{$\boldsymbol{\mathbf{a_{r,\mu}}}$}\right)`,
:math:`\mathbf{A}\left(\mbox{$\boldsymbol{\mathbf{a_{r,\epsilon}}}$}\right)`,
and
:math:`\mathbf{A}\left(\mbox{$\boldsymbol{\mathbf{a_{r,\rho}}}$}\right)`
for the :math:`\boldsymbol{\mathbf{a_{r,\mu}}}`, :math:`\boldsymbol{\mathbf{a_{r,\epsilon}}}`, and :math:`\boldsymbol{\mathbf{a_{r,\rho}}}` vectors defined for codons immediately above. The
likelihoods are computed by similarly transforming the count vectors
in Equations :eq:`pr_nrpre`, :eq:`pr_nrpost`,
:eq:`pr_nrerrpre`, and :eq:`pr_nrerrpost` to
:math:`\mathbf{A}\left(\mbox{$\mathbf{n_r^{\textbf{pre}}}$}\right)`,
:math:`\mathbf{A}\left(\mbox{$\mathbf{n_r^{\textbf{post}}}$}\right)`,
:math:`\mathbf{A}\left(\mbox{$\mathbf{n_r^{\textbf{err,pre}}}$}\right)`,
and
:math:`\mathbf{A}\left(\mbox{$\mathbf{n_r^{\textbf{err,post}}}$}\right)`.

.. _chartype_aa:

Characters are amino acids; preferences are for amino acids.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
A fourth scenario is that the deep sequencing counts are for amino-acid
characters :math:`x`, and that we want to determine the preference
:math:`\pi_{r,x}` for each amino acid at each site :math:`r`. This scenario
would only be used if there has already been some post-processing of the deep mutational
scanning data, since the actual sequence data will always report either
nucleotides or codons. If at all possible, you should prefer to analyze the data
using the approach :ref:`chartype_codon_to_aa`. However, in some cases you may receive third-party
data that has already been processed to amino acids. In this scenario,
we use the prior assumption that all amino acid mutations are introduced at equal frequency,
and that all amino-acid errors happen at equal frequency. This first assumption will
only be true for certain mutagenesis strategies. The second assumption is very unlikely
to be true -- but we can't do better without analyzing the data at the codon level.
So in this scenario, we use the prior
vectors specified by Equations :eq:`arepsilon`, :eq:`arrho`, and :eq:`armu`,
with the summations over :math:`x` covering the 20 amino acids (or 21 if stop codons
are allowed),
and the summations over :math:`m` covering
the two possible number of amino-acid changes to an amino-acid character (0 and 1).


.. _MCMC_inference:

Inferring the preferences by MCMC.
---------------------------------------------------------

Given the sequencing data and the parameters that specify the priors,
the posterior probability of any given set of the parameters :math:`\boldsymbol{\mathbf{\pi_r}}`, :math:`\boldsymbol{\mathbf{\mu_r}}`, :math:`\boldsymbol{\mathbf{\epsilon_r}}`,
and :math:`\boldsymbol{\mathbf{\rho_r}}`
is given by the product of Equations :eq:`pr_nrpre`,
:eq:`pr_nrpost`, :eq:`pr_nrerrpre`,
:eq:`pr_nrerrpost`, :eq:`pr_pir`, :eq:`pr_mur`,
:eq:`pr_epsilonr`, and :eq:`pr_rhor`. We use `MCMC`_
implemented by `PyStan`_ 
to sample from this posterior distribution. We
monitor for convergence of the sampling of the site-specific preferences
 by running four chains, and ensuring that the mean over all :math:`\pi_{r,x}`
values of the potential scale
reduction statistic :math:`\hat{R}` of `GelmanRubin1992`_ is :math:`\le 1.1` and that the
mean effective sample size is :math:`\ge 100`, or that :math:`\hat{R} \le 1.15` and the effective sample size is :math:`\ge 300`. We repeatedly increase the number of steps 
until convergence occurs. The inferred preferences are
summarized by the posterior mean and median-centered 95% credible
interval for each :math:`\pi_{r,x}`.

.. include:: weblinks.txt
