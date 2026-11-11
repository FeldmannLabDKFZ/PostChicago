# PostChicago 0.99.6

## Vignette improvements
* added additional description and explanation regarding k to vignette and plotInteractions()
* included a detailed comparison with Chicago output
* included more view points

## Minor improvements
makeQCplots() now option to display numbers and individual interaction hetamaps limited to maximum 50000 interactions

# PostChicago 0.99.5

* added citation to the new preprint describing the package

# PostChicago 0.99.4

## Bug fixes
* fixed duplicate rows in getMatrix

## Improvements
* sped up getMatrix()

# PostChicago 0.99.3

## Bug fixes
* Fixed oneGeneOnePeak handling of gaps between fragments
* Fixed the handling of name extraction for the files loaded in loadGrl()

## Minor improvements
* Added the option of also creating sample specific interaction bed files via makeBedGraphs 
* Changed matrix normalization to a reference sample in getMatrix to normalize to reference sample summit
* changed plotSigIntsStats to plotStats to simplify naming

# PostChicago 0.99.2

* Initial version submitted to Bioconductor