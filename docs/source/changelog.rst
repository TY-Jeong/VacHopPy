.. _changelog:

Release History
===============

Version 3.1.1 (January 2, 2026)
-------------------------------

**Documentation**

* Updated the citation instructions in the documentation.

Version 3.1.0 (October 29, 2025)
--------------------------------

**New Features**

* Added an ``--error_bar`` option to the plot command to display the standard error of the mean (SEM) on various plots (:math:`D`, :math:`D_{rand}`, :math:`f`, :math:`\\tau`).

**Fixes**

* Fixed an issue where the CLI was not available when ``VacHopPy`` was installed via ``pip``.
* Fixed a segmentation fault with NumPy (``>=2.3``) and SciPy (``>=1.16``).

**Changes**

* Updated default values for ``cos_margin`` and ``force_margin``.
* Disabled automatic opening of the HTML output in the ``trajectory`` command to avoid errors on some systems.

Version 3.0.0 (October 18, 2025)
--------------------------------

Initial public release of ``VacHopPy``.

* Core analysis engine for vacancy trajectory tracing.
* Ensemble analysis for multi-temperature simulations.
* Tools for calculating diffusivity, correlation factors, and attempt frequencies.
