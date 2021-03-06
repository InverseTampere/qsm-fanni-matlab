# Changelog

## 2019-07-12 Version 1.3.0

- Separated QSM-Blocks as a separate repository as it is now shared with various projects.
	- QSMB, QSMBCylindrical and CubeVoxelization depencies must be downloaded manually and added to the path.
- Added warning for when leaf candidate area is smaller than the leaf target area in the qsm_fanni function.
- Moved lisence boilerplate after the main function help text in qsm_fanni.m
- Set float format from *f* to *g* to exclude trailing zeros when exporting.
- Simplified *test.m* as QSMBCylindrical can now be initialized with the 'example' for creating the simple QSM tree.
- Added a `trim_slack()` call at the end of the qsm_fanni function to clear empty rows in object properties.
- Fixed *leaf_parent* property orientation from row to column vector.

## 2018-04-17

- Fixed leaf scale bug that caused scaling to be done twice when computing exact leaf geometry.
- Updated README with a link to a demonstration video.

## 2018-02-13 Version 1.2.0

- Added second leaf export format: *Extended OBJ*.
- Added annotations to *test.m*.

## 2017-12-08 Version 1.1.1

- Fixed bug where the `draw_cylinders()` method of the `QSMBCylindrical` class did not have access to the `rotation_matrix()` function.
- Extended the documentation to include a detailed description of the test case defined in *src/test.m* file.
- Added syntax highlighting in *README.md*.
- Added accepted leaf geometry computation to the example in *src/test.m* file.

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