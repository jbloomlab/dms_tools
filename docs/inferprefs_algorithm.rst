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

   \Pr\left(\mbox{$\boldsymbol{\mathbf{\pi_r}}$}\right) &=& \operatorname{Dirichlet}\left(\alpha_{\pi} \times \mathbf{1}\right)

.. math::
   :label: pr_mur

   \Pr\left(\mbox{$\boldsymbol{\mathbf{\mu_r}}$}\right) &=& \operatorname{Dirichlet}\left(\alpha_{\mu} \times \mathcal{N}_x \times \mbox{$\boldsymbol{\mathbf{a_{r,\mu}}}$}\right)

.. math::
   :label: pr_epsilonr

   \Pr\left(\mbox{$\boldsymbol{\mathbf{\epsilon_r}}$}\right) &=& \operatorname{Dirichlet}\left(\alpha_{\epsilon} \times \mathcal{N}_x \times \mbox{$\boldsymbol{\mathbf{a_{r,\epsilon}}}$}\right) 

.. math::
   :label: pr_rhor

   \Pr\left(\mbox{$\boldsymbol{\mathbf{\rho_r}}$}\right) &=& \operatorname{Dirichlet}\left(\alpha_{\rho} \times \mathcal{N}_x \times \mbox{$\boldsymbol{\mathbf{a_{r,\rho}}}$}\right)

where :math:`\mathbf{1}` is a vector of ones, :math:`\mathcal{N}_x` is
the number of characters (i.e. 64 for codons, 20 for amino acids, 4 for
nucleotides), the :math:`\alpha` parameters (i.e. :math:`\alpha_{\pi}`,
:math:`\alpha_{\mu}`, :math:`\alpha_{\epsilon}`, and
:math:`\alpha_{\rho}`) are scalar concentration parameters with values
:math:`> 0` specified by the user, and the :math:`\mathbf{a_r}` vectors
(i.e. , , ) are vectors with entries :math:`> 0` that sum to one.

The method used to specify the prior vectors :math:`\mathbf{a_{r,\mu}}`, :math:`\mathbf{a_{r,\epsilon}}`, and :math:`\mathbf{a_{r,\rho}}` depends on the
type of data being analyzed. We consider three possibilities:

.. _chartype_DNA:

Characters are nucleotides; preferences are for nucleotides.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
One possibility is that the deep sequencing counts are for nucleotide
characters :math:`x`, and that we want to determine the preference
:math:`\pi_{r,x}` for each nucleotide at each site :math:`r`. This is
the most natural approach if the mutant library is generated by a
nucleotide-level mutagenesis process, such as error-prone PCR. We use
priors based on the expectation is that all possible
single-nucleotide mutations are made at equal frequency in the mutant
library, and that all single-nucleotide errors are equally likely. So
we define

.. math::
   :label: avgmu

   \overline{\mu} = \sum\limits_r \frac{1}{\mbox{$N_r^{\textrm{pre}}$}}\sum\limits_{x\ne \operatorname{wt}\left(r\right)} \mbox{$n_{r,x}^{\textrm{pre}}$}

as the average mutation rate in the pre-selection mutant library,

.. math::
   :label: avgepsilon

   \overline{\epsilon} = \sum\limits_r \frac{1}{\mbox{$N_r^{\textrm{err,pre}}$}}\sum\limits_{x\ne \operatorname{wt}\left(r\right)} \mbox{$n_{r,x}^{\textrm{err,pre}}$}

as the average error rate in the pre-selection error control, and

.. math::
   :label: avgrho

   \overline{\rho} = \sum\limits_r \frac{1}{\mbox{$N_r^{\textrm{err,post}}$}}\sum\limits_{x\ne \operatorname{wt}\left(r\right)} \mbox{$n_{r,x}^{\textrm{err,post}}$}

as the average error rate in the post-selection error control. We
then define

.. math::
   :label: ntarmu

   \mbox{$\boldsymbol{\mathbf{a_{r,\mu}}}$}&=& \left(\cdots, \frac{\overline{\mu}}{\mathcal{N}_x - 1} + \delta_{x,\operatorname{wt}\left(r\right)} \times \left[1 - \overline{\mu}\right] ,\cdots\right),

.. math::
   :label: ntarepsilon

   \mbox{$\boldsymbol{\mathbf{a_{r,\epsilon}}}$}&=& \left(\cdots, \frac{\overline{\epsilon}}{\mathcal{N}_x - 1} + \delta_{x,\operatorname{wt}\left(r\right)} \times \left[1 - \overline{\epsilon}\right] ,\cdots\right),

.. math::
   :label: ntarrho

   \mbox{$\boldsymbol{\mathbf{a_{r,\rho}}}$}&=& \left(\cdots, \frac{\overline{\rho}}{\mathcal{N}_x - 1} + \delta_{x,\operatorname{wt}\left(r\right)} \times \left[1 - \overline{\rho}\right] ,\cdots\right).

.. _chartype_codon:

Characters are codons; preferences are for codons.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
A second
possibility is that the deep sequencing counts are for codon
characters :math:`x`, and that we want to determine the preference
:math:`\pi_{r,x}` for each codon at each site :math:`r`. This is a
natural approach if the mutant library is generated by introducing
all possible codon mutations at each site, as is done by techniques
like the one described by `Firnberg2012`_. We use priors based on the expectation
that all possible codon mutations are made at equal frequency. Note
that these priors are only appropriate if codons are mutagenized to
``NNN`` (they would need to be adjusted if codons are only mutated to
a subset of possibilities such as ``NNK``). Assuming that all codon
mutations are created, we define :math:`\boldsymbol{\mathbf{a_{r,\mu}}}` by Equations :eq:`avgmu` and
:eq:`ntarmu` above, with the summations over :math:`x` now
encompassing all possible codons. We assume that sequencing errors
occur at the nucleotide level, so the prior expectation for the error
rate depends on the number :math:`N_{x,\operatorname{wt}\left(r\right)}` of
nucleotide changes between codon :math:`x` and the wildtype codon  at
site :math:`r`. Specifically, let

.. math::
   :label: avgepsilonm

   \overline{\epsilon_m} &=& \sum\limits_r \frac{1}{\mbox{$N_r^{\textrm{err,pre}}$}}\sum\limits_{x} \mbox{$n_{r,x}^{\textrm{err,pre}}$}\times \delta_{m,N_{x,\operatorname{wt}\left(r\right)}} 

.. math::
   :label: avgrhom

   \overline{\rho_m} &=& \sum\limits_r \frac{1}{\mbox{$N_r^{\textrm{err,post}}$}}\sum\limits_{x} \mbox{$n_{r,x}^{\textrm{err,post}}$}\times \delta_{m,N_{x,\operatorname{wt}\left(r\right)}}

where :math:`m` can be 0, 1, 2, or 3 (the possible number of
nucleotide changes to a codon); note that
:math:`1 = \overline{\epsilon_0} + \overline{\epsilon_1} + \overline{\epsilon_2} + \overline{\epsilon_3} = \overline{\rho_0} + \overline{\rho_1} + \overline{\rho_2} + \overline{\rho_3}`.
We then define

.. math::
   :label: codonarepsilon

   \mbox{$\boldsymbol{\mathbf{a_{r,\epsilon}}}$}&=& \left(\cdots, \sum\limits_{m=0}^3 \frac{\overline{\epsilon}_{m}}{\mathcal{C}_m} \times \delta_{m,N_{x,\operatorname{wt}\left(r\right)}},\cdots\right),

.. math::
   :label: codonarrho

   \mbox{$\boldsymbol{\mathbf{a_{r,\rho}}}$}&=& \left(\cdots, \sum\limits_{m=0}^3 \frac{\overline{\rho}_{m}}{\mathcal{C}_m} \times \delta_{m,N_{x,\operatorname{wt}\left(r\right)}}, \cdots\right)

where :math:`\mathcal{C}_m` is the number of codons with :math:`m`
changes relative to some wildtype codon (so
:math:`\mathcal{C}_0 = 1`, :math:`\mathcal{C}_1 = 9`,
:math:`\mathcal{C}_2 = \mathcal{C}_3 = 27`).

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

Inferring the posterior probability over the preferences.
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
 by running three chains, and ensuring that the potential scale
reduction statistic :math:`\hat{R}` of `GelmanRubin1992`_ is :math:`\le 1.05` and that the
effective sample size is :math:`\ge 200` (the number of steps is
increased until convergence occurs). The inferred preferences are
summarized by the posterior mean and median-centered 95% credible
interval for each :math:`\pi_{r,x}`.

.. include:: weblinks.txt
