# SoS-ZCP

 [![DOI](https://zenodo.org/badge/374856817.svg)](https://zenodo.org/badge/latestdoi/374856817) [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4968507.svg)](https://doi.org/10.5281/zenodo.4968507)


SoS-ZCP is an algorithm that uses the exact ray tracing to evaluate the transverse coma aberration of an optical system. Our algorithm applies the Newton-Raphson method  to calculate the value of tilt that minimizes the amount of coma produced by a secondary mirror that is lateral shifted . This code was developed to calculate the zero coma point for every one of the classical telescopes of the OAN-SPM at Sierra San Pedro Mártin (B.C., México), but is suitable to be used in non-classical telescopes, like those of the TAOS-II project which is in development at the OAN-SPM. The case of the TAOS-II is important because the classical solution for the calculation of the zero coma point does not take correcting aspherical surfaces into account. This program has been tested with Python 2.7 and the libraries numpy 1.8.0 rc1 and matplotlib 1.3.1.

## Features
The main goal of this algorithm is to calculate the zero coma point not only for classical telescopes, but also for non-classical telescopes. You can do the following:
* Define the number of surfaces according to your requirements.
* As we work with a numerical derivative to find the point of intersection and the normal vector, you can use any surface that you need. You only have to define the function.
* Save the points of the intersection and their director cosines to use this information later.
* Calculate the zero coma point and the angle that compensates the aberration for a given translation of the secondary mirror.

You will find three examples of the algorithm that illustrates three different problems. The first example is the calculation of the zero coma point for a Non-Classical Telescope. The second example, we do the same but now for a Classical Telescope. The last example shows the compensation of 9 different fields.

## Name

We named the algorithm as SoS-ZCP, SoS for Serpents of Stars because we use vectors to represent the rays of lights and we can see them like serpents. Citlacoatl means “Serpents of Stars” and the M. A. Patricia Ávila (Avila, 2014) tells us the history of a great warrior who had to fulfill his mission of giving vitality and building up the great Mexica Empire.

## Citation

## Documentation
* Avila, P. L. (2014). Serpiente de estrellas. Createspace Independent Publishing Platform.
