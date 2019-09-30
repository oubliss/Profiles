.. Profile_Set documentation master file, created by
   sphinx-quickstart on Tue Feb 26 09:13:51 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Profiles: powered by CASS
=======================================

.. toctree::
   :maxdepth: 2
   :caption: Contents:
   :titlesonly:

   Profile
   Raw_Profile
   Thermo_Profile
   utils
   Profile_Set
   Wind_Profile



Key Equations
==================

1) Position data to wind

.. math::

   \begin{aligned}
   \mathrm{Speed} &= A \sqrt{\tan{[\arccos{(\cos{(\mathrm{PITCH})} - \cos{(\mathrm{ROLL})})}]}} + B \\ 
   \\
   \mathrm{Direction} &= \arctan{\Bigg[\frac{-\cos{(\mathrm{YAW})}\sin{(\mathrm{ROLL})}}{\sin{(\mathrm{YAW})}\sin{(\mathrm{ROLL})} + \sin{(\mathrm{PITCH})}\cos{(\mathrm{YAW})}\cos{(\mathrm{ROLL})}} \Bigg]}
   \end{aligned}


2) Thermistor resistance to temperature

.. math::

   \begin{aligned}
   T &= \frac{1}{A + B (log{R}) + C (log{R})^3}
   \end{aligned}
