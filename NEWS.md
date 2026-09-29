# methylTFRAnnotationBuilder 0.99.2

* Removed the legacy GC helpers `calculate_gcdist()`, `get_gcdist()`,
  `compute_gc()`, `convert_to_bins()` and `convert_to_matrix()`. They
  were superseded by `processMotifs2Matrix()` and `computeGCgenome()`,
  and nothing in the package used them.
* README no longer describes a CpG-restricted genome GC table
  (`cpgSites()` / `gc_sites`); annotations use the full genome-wide
  table.
* `createMethylTFRPackageScaffold()` sets the organism biocView from the
  assembly instead of always writing `Homo_sapiens`.
* Documentation cleaned up; `DESCRIPTION` lists all optional packages
  used by `getGenomeObject()` and `prepareMotifmatchr()`, and the
  funders with the `fnd` role.
* Added `CITATION.cff`.
