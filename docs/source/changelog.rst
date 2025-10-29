.. _changelog:

Release History
===============

Version 3.1.0 (October 29, 2025)
-------------------------------------

**New Features**

* Added an ``--error_bar`` option to the plot command to display the standard error of the mean (SEM) on various plots (:math:`D`, :math:`D_{rand}`, :math:`f`, :math:`\tau`).

**Changes & Fixes**

* Fixed an issue where the CLI was not available when ``VacHopPy`` was installed via ``pip``.
* Updated default values for ``cos_margin`` and ``force_margin`` parameters.
* Removed the automatic opening of the HTML file in the ``trajectory`` command to prevent errors on certain systems.


====

Version 3.0.0 (October 18, 2025)
---------------------------------

Initial public release of ``VacHopPy``.

* Core analysis engine for vacancy trajectory tracing.
* Ensemble analysis for multi-temperature simulations.
* Tools for calculating diffusivity, correlation factors, and attempt frequencies.
