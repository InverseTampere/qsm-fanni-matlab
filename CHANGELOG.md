# Changelog

## 2017-11-22 Version 1.1.0

- Crushed cylinder-triangle intersection detection bugs.
- Added option to detect edge hits.
- Updated *test.m* file to include leaf over-sampling, to reach the goal leaf area.
- Updated *test.m* resulting export files and rendered image.
- Updated the constructor of `QSMBCylindrical` to accept struct-based model data, compatible with [TreeQSM v2.30](https://github.com/InverseTampere/TreeQSM/releases).
- Updated `LeafModelTriangle` class with additional methods:
	- `compute_geometry(...)` is a method for computing leaf vertices and faces based on the leaf geometry and either transformation parameters given as input, or parameters of the accepted leaves.
	- `export_geometry(...)` is a generalization of the old `export_obj(...)` method, allowing new export formats to be added to the same method in the future.
	- Other methods have been updated to accommodate the new methods, as well as, the demonstration file *test.m*.
- Wrote documentation in *README.md*.

## 2017-05-30 Version 1.0.0-alpha

- Initial pre-release version.
- Well commented in source, but lacking proper documentation.