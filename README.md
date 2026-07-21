![](https://raw.githubusercontent.com/Enserink-lab/screenwerk-2.0/refs/heads/develop/doc/logo.png)

# screenwerk: a modular pipeline for drug sensitivity screens

------------------------------------------------------------------------

An R package with a modular library for the design and analysis of high-throughput drug combination screens.

code name: screenwerk dev: 2.2.2-1 author: roberthanes

![](https://raw.githubusercontent.com/Enserink-lab/screenwerk-2.0/refs/heads/main/doc/screenwerk-2.2.png)

The R-package consists of individual modules that can be used to (a) set-up a screen and generate a dispensing file, (b) to read raw measurements and consolidate different data sets, (c) perform a quality control, (d) analyze, (e) report and visualize the results of a drug (combination) screen.

### INSTALLATION

**screenwerk** can be installed directly from github using devtools.

```{r screenwerk, eval = FALSE}
install.packages('devtools')
library(devtools)

# Install R-package from github
devtools::install_github('Enserink-lab/screenwerk-2.0', build = TRUE, build_opts = c("--no-resave-data", "--no-build-vignettes"))
# Load 'screenwerk' package
library(screenwerk)
```

### Example files and data

A set of example files can be found in inst/extdata/ and used as a template for building individual data sets to set up and run a custom drug sensitivity screen.

The folder inst/extdata/**library** contains examples of input files that have been used to set up a explanatory drug sensitivity screen. The folder inst/extdata/**raw** contains a set of raw measurement files based on the dispensing layout generated during the set up of the explanatory drug sensitivity screen.

In addition, a set of `.rds` files (R objects) can be found in inst/extdata/**rds** illustrating the data structure and data output of individual modules.

### **Example results and plots**

The output from individual modules based on the example files above, can be found in the folder doc/results and used as a reference on what output files are being generated. A set of csv files can be found in doc/results/**Dynamic Drug Activity Range**, doc/results/**EC50s** and doc/results/**synergyscores**. A set of plots can be found in the folder oc/results/**graphs**.

The output from the quality control module can be found in the folder doc/results/**quality control**.

The graphical output illustrating the dispensing layout can be found in the folder doc/**export**/**plates**.

### Workflow

Each module consists of an individual function that can be used to accomplish a certain tasks. Those modules are divided into (a) pre-experimental and (b) post-experimental modules. The pre-experimental modules are used to initiate the set up of a drug screen, design and generate dispensing files with instructions for a dispensing machine/robot. The post-experimental modules are used once a drug screen has been run in the lab. Those modules are used to handle different data sets, such as the consolidation of the dispensing or experimental layout with the raw measurements from the drug screen, normalize the data and prepare the data set for downstream analysis. Once the data has been processed, it can be used first by analytical modules for the analysis of the drug screen, followed by the use of reporting modules for visual and statistical reporting of the results from the drug sensitivity screen.

### *(a) pre-experimental: setting-up a screen and generating a dispensing file*

In order to be able to generate a dispensing file for a drug combination screen, we need a list of combinations, among a few other files. All files are provided as .csv files. The list of combinations can be easily generated using one of the functions in this package.

First we need a list of drugs with all the doses to be used in the screen, which will serve as a list of doses to be combined with each other at selected doses. This list can be provided in a wide format, as it would be done using Microsoft Excel, or any other spreadsheet tool.

![](https://raw.githubusercontent.com/Enserink-lab/screenwerk-2.0/refs/heads/develop/doc/figures/doses-wide.png)

Now we need to import this file into R using the code below:

```{r, eval = FALSE}
listofDoses <- read.csv(file=file.path("inst/extdata/library/listofdoses.csv"), check.names=FALSE, header=TRUE, stringsAsFactors=FALSE, na.strings="", sep=",", dec=".", skip=0)
```

and run the following function, which converts the list from a human readable wide-format to a long-format:

```{r, eval = FALSE}
listofDoses <- generateListofDoses(listofDoses, .doseIdentifier = "dose", .dropCol = TRUE)
```

It is possible to specify the identifier for the doses with the argument *.doseIdentifier*.

In the example above with the imported list of doses, the doses are provided in the columns labelled "1st Dose", "2nd Dose", ... "6th Dose". With the argument *.doseIdentifier = "dose"*, all columns containing the case-insensitive word "dose" will be identified as columns containing the different doses for each of the drugs. In addition, the argument *.dropCol = TRUE* drops columns that are not needed, but might be present in the the file.

This will generate a data set that looks like the one below:

![](https://raw.githubusercontent.com/Enserink-lab/screenwerk-2.0/refs/heads/develop/doc/figures/doses-long.png)

**Combining drugs: creating a list of combinations**

Before we can generate the drug combinations from the generated list of doses, we need to import a list of drugs as shown below:

```{r, eval = FALSE}
listofDrugs <- read.csv(file=file.path("inst/extdata/library/listofdrugs.csv"), check.names=FALSE, header=TRUE, stringsAsFactors=FALSE, na.strings="", sep=",", dec=".", skip=0)
```

![](https://raw.githubusercontent.com/Enserink-lab/screenwerk-2.0/refs/heads/develop/doc/figures/drugs.png)

Now that we have a list of drugs and a list of doses, we can use both to generate a list of combinations:

```{r, eval = FALSE}
listofCombinations <- combineDrugs(
  listofDrugs, listofDoses,
  .combineDoses=c(2:5),
  .noReplicates = 3, .drugRepAttrib = "single",
  .pairBy = "group",
  .inclusive = NULL, .exclusive = NULL
)
```

In the function above, we can provide the list of drugs together with the list of doses we just generated. With *.combineDoses* we can specify, which of the doses we want to combine and therby define the n-by-n combination matrix. In the example above with *.combineDoses=c(2:5)* we select doses 2 to 5 to be combined out of a total of 6 doses, leaving out the lowest and highest dose. This particular selection of doses will create a 4-by-4 combination matrix.

![](https://raw.githubusercontent.com/Enserink-lab/screenwerk-2.0/refs/heads/develop/doc/figures/n-by-n-design.png)

In order to run a full combination matrix, we would have to select all doses with *.combineDoses=c(1:6)*, leading to a 6-by-6 combination matrix.

The same function above can be used to set up a (a) single drug or (b) drug combination screen through the argument *.combineDoses*.

**Setting up a single drug screen**

In order to set up a single drug screen, the same function can be used with *.combineDoses = FALSE*. This will only create a list of single drug treatments.

```{r, eval = FALSE}
listofCombinations <- combineDrugs(
  listofDrugs, listofDoses,
  .combineDoses=FALSE,
  .noReplicates = 3, .drugRepAttrib = "all"
)
```

If the *.combineDoses* argument is missing, a drug combination screen will be set up with all the available doses combined.

![](https://raw.githubusercontent.com/Enserink-lab/screenwerk-2.0/refs/heads/develop/doc/figures/combinations-single.png)

**Setting up a drug combination screen**

In order to set up a drug combination screen, the specific doses to be combined need to be specified through *.combineDoses.*

```{r, eval = FALSE}
listofCombinations <- combineDrugs(
  listofDrugs, listofDoses,
  .combineDoses=c(2:5),
  .noReplicates = 3, .drugRepAttrib = "single"
  )
```

If the doses to be combined are not specified, by leaving out the argument from the function call, then the entire dose range will be combined.

```{r, eval = FALSE}
listofCombinations <- combineDrugs(
  listofDrugs, listofDoses,
  .noReplicates = 3, .drugRepAttrib = "single"
  )
```

![](https://raw.githubusercontent.com/Enserink-lab/screenwerk-2.0/refs/heads/develop/doc/figures/combinations-double.png)

**Creating replicates for single drug treatments or both, single drug and drug combination treatments**

With *.noReplicates* we can specify how many replicates we want, and with *.drugRepAttrib* we can specify, if only single drug treatments or both, single drug and drug combination treatments should be replicated. In th example above, we chose with *.noReplicates = 3* and *.drugRepAttrib = "single"* to generate triplicates for all single drug treatments only.

Alternatively, this argument can be set to *.drugRepAttrib = "all"*, in which case all drug treatments will be replicated at the given number.

```{r, eval = FALSE}
listofCombinations <- combineDrugs(
  listofDrugs, listofDoses,
  .combineDoses=c(2:5),
  .noReplicates = 3, .drugRepAttrib = "all"
)
```

**Combining drugs based on solvent or individual groups**

In addition, the pairing of drugs can be individually selected with the argument *.pairBy*, which allows to combine drugs either by solvent *.pairBy = "solvent"* or by custom groups.

```{r, eval = FALSE}
listofCombinations <- combineDrugs(
  listofDrugs, listofDoses,
  .combineDoses=c(2:5),
  .noReplicates = 3, .drugRepAttrib = "single",
  .pairBy = "solvent",
)
```

In the case of grouping by solvent, drugs dissolved in the same solvent (belonging to the same solvent group) will be combined with each other, but not with drugs belonging to a different solvent group. This grouping ensures that drugs are not combined (dispensed into the same well) with incompatible solvents that migtht have an adverse impact on the activity and effect of the other drug.

In any other case, drugs can be paired by individual groups specified by the user. This allows selective and custom combinations of drugs. This can be accomplished by providing a column in the list of drugs with a custom column $$'GROUP'$$ containing the individual groups, as shown in an example below.

![](https://raw.githubusercontent.com/Enserink-lab/screenwerk-2.0/refs/heads/develop/doc/figures/drug-groups.png)

In the case above, drugs belonging to the same group (column containing the same group) will be combined with each other, but not with any other drugs that belong to a different group. In other words, drugs that have matching groups will be combined.

```         
⟁ Note: In the figue shown above, drug with the number '1' belonging to the group 'A' will be comined with every other drug that has or contains the same group as drug '1'. In that example, drug '1' will be combined with drug number '20', '21', '22', '30' and '54', but not with all the other drugs having different groups. Furthermore, drugs with identical groups will also be combined with each other. In that case, drugs '20', '21','22, '30' and '54' will not only be combined with all the other drugs, but also with each other, since those drugs not only contain all other groups, but are also identical in group to each other, respectively. 
```

```{r, eval = FALSE}
listofCombinations <- combineDrugs(
  listofDrugs, listofDoses,
  .combineDoses=c(2:5),
  .noReplicates = 3, .drugRepAttrib = "single",
  .pairBy = "group",
)
```

**Including and excluding individual drugs from being combined**

The arguments *.inclusive* and *.exclusive* allow optionally only selected drugs to be either combined or excluded from being combined. In the first case, only drugs that are listed in the argument *.inclusive* will be combined, all other drugs will be ignored.

In contrast, all drugs that are listed in the argument *.exclusive* will be excluded from being combined. If the same drugs should be listed both as *.inclusive* and *.exclusive*, the later argument will be ignored and the inclusion of the drug prioritized.

**Reducing the design of a drug combination matrix**

Now that we have generated a list of combinations based either on a full or limited n-by-n matrix design, we can further reduce the number of drug combination treatments by filtering only drug treatments that fall into a given matrix pattern. As of now, only the x-design and all diagonal designs are supported. The x-design is still kept as the default, as it is the only pattern with the highest accuracy despite a reduction in the number of drug treatments for a given drug pair. This is due to the fact, that the x-design is the only design that covers the entire dose range, compared to most other designs.

The full matrix design can be reduced by the following function:

```{r, eval = FALSE}
listofCombinations <- reduceDesign(listofCombinations, .design="x")
```

![](https://raw.githubusercontent.com/Enserink-lab/screenwerk-2.0/refs/heads/develop/doc/figures/reduced-design.png)

In order to use a diagonal pattern, the the **'Diagonal'** matrix design can be selected in two directions, from top-left to bottom-right with either one of the arguments **'Diagonal:Top-Left'** or **'Diagonal:Bottom-Right'** or from top-right to bottom-left with either one of the arguments **'Diagonal:Top-Right'** or **'Diagonal:Bottom-Left'**

```{r, eval = FALSE}
listofCombinations <- reduceDesign(listofCombinations, .design="Diagonal:Top-Left")
# or
listofCombinations <- reduceDesign(listofCombinations, .design="Diagonal:Bottom-Right")
```

```{r, eval = FALSE}
listofCombinations <- reduceDesign(listofCombinations, .design="Diagonal:Top-Right")
# or
listofCombinations <- reduceDesign(listofCombinations, .design="Diagonal:Bottom-Left")
```

The reduced design is optional and only applied when the space on the plate is limited, or if the experimental design has limitation in terms of resources, i.e. drugs or samples. Otherwise, it is always recommended to run a full treatment design.

```         
⟁ Note: depending on the size of the screen and the number of combinations, this task can take some time to successfully finish. The status of reducing the drug treatments will be shown on screen throughout the process.
```

Now we are ready to use the list of combinations we just generated, to build a dispensing data set. However, in order to do that, we need to provide a few additional files.

**Excluding well from the plate map layout**

First, let's generate a list of wells we want to exclude from being dispensed into. This can be achieved using the function below:

```{r, eval = FALSE}
listofExWells <- excludeWells(1536, outer.wells=TRUE)
```

This will generate the following list of wells:

```{r, eval = FALSE}
listofExWells
```

``` r
> listofExWells
[1] "A1"   "A2"   "A3"   "A4"   "A5"   "A6"   "A7"   "A8"   "A9"   "A10"  "A11"  "A12"  "A13"  "A14"  "A15"  "A16"  "A17"  "A18"  "A19"  "A20"  "A21"  "A22"  "A23"  "A24"  "A25"
[26] "A26"  "A27"  "A28"  "A29"  "A30"  "A31"  "A32"  "A33"  "A34"  "A35"  "A36"  "A37"  "A38"  "A39"  "A40"  "A41"  "A42"  "A43"  "A44"  "A45"  "A46"  "A47"  "A48"  "B1"   "B48"
[51] "C1"   "C48"  "D1"   "D48"  "E1"   "E48"  "F1"   "F48"  "G1"   "G48"  "H1"   "H48"  "I1"   "I48"  "J1"   "J48"  "K1"   "K48"  "L1"   "L48"  "M1"   "M48"  "N1"   "N48"  "O1"  
[76] "O48"  "P1"   "P48"  "Q1"   "Q48"  "R1"   "R48"  "S1"   "S48"  "T1"   "T48"  "U1"   "U48"  "V1"   "V48"  "W1"   "W48"  "X1"   "X48"  "Y1"   "Y48"  "Z1"   "Z48"  "AA1"  "AA48"
[101] "AB1"  "AB48" "AC1"  "AC48" "AD1"  "AD48" "AE1"  "AE48" "AF1"  "AF2"  "AF3"  "AF4"  "AF5"  "AF6"  "AF7"  "AF8"  "AF9"  "AF10" "AF11" "AF12" "AF13" "AF14" "AF15" "AF16" "AF17"
[126] "AF18" "AF19" "AF20" "AF21" "AF22" "AF23" "AF24" "AF25" "AF26" "AF27" "AF28" "AF29" "AF30" "AF31" "AF32" "AF33" "AF34" "AF35" "AF36" "AF37" "AF38" "AF39" "AF40" "AF41" "AF42"
[151] "AF43" "AF44" "AF45" "AF46" "AF47" "AF48"
```

This function is complementary and returns a list of wells to be excluded from a 6, 12, 24, 48, 96, 384 or 1536-well microplate in an experimental setup of a high-throughput single or drug combination sensitivity screen (DSS).

In the first argument, we specify the plate format we plan to use, usually denoted by the number of wells. With the argument *wells* we can list all the wells that we want to exclude from the experimental set-up. This can be achieved by providing the column or row name, which are designated letters for rows (A, B, C, ..) and designated numbers for columns (1, 2, 3, ..). Alternatively, in addition to entire rows and columns, individual wells (A1, B2, C3, ..) can be marked for exclusion.

In order to simplify the exclusion of wells falling into specific patterns, such as outer wells of a microplate, the argument *outer.wells* is TRUE, can be used. In that case all outer wells of given plate type will be marked for exclusion.

```         
⟁ Note: allowing outer.wells to be TRUE for a 6-well plate will consequently mark all wells for exclusion on that plate!
```

Any list of wells to be excluded from an experimental set up can also be imported by a simple csv file, instead of being generated. this can be useful, if always the same wells are being excluded from an experimental set up of a drug sensitivity screen.

A list of excluded wells provided as a csv file could look like this:

```         
A1
A2
A3
A4
A5
A6
A7
A8
A9
A10
A11
A12
A13
A14
A15
A16
A17
A18
A19
A20
A21
A22
A23
A24
```

Now that this is done, we need to provide the following files, in addition to the ones already generated:

- a list of drugs

- a list of volumes

- a list of controls

- a list of stock concentrations

- one or more source plates, optional

The list of drugs can look as shown below:

![](https://raw.githubusercontent.com/Enserink-lab/screenwerk-2.0/refs/heads/develop/doc/figures/drugs.png)

The list of volumes, can be provided same as the list of doses in a wide-format:

![](https://raw.githubusercontent.com/Enserink-lab/screenwerk-2.0/refs/heads/develop/doc/figures/volumes.png)

The list of controls:

![](https://raw.githubusercontent.com/Enserink-lab/screenwerk-2.0/refs/heads/develop/doc/figures/controls.png)

The list of stock concentrations:

![](https://raw.githubusercontent.com/Enserink-lab/screenwerk-2.0/refs/heads/develop/doc/figures/stock-concentrations.png)

Once we have all files in place, we can import them:

```{r, eval = FALSE}
listofDrugs <- read.csv(file=file.path("inst/extdata/library/listofdrugs.csv"), check.names=FALSE, header=TRUE, stringsAsFactors=FALSE, na.strings="", sep=",", dec=".", skip=0)
listofVolumes <- read.csv(file=file.path("inst/extdata/library/listofvolumes.csv"), check.names=FALSE, header=TRUE, stringsAsFactors=FALSE, na.strings="", sep=",", dec=".", skip=0)
listofCtrls <- read.csv(file=file.path("inst/extdata/library/listofctrls.csv"), check.names=FALSE, header=TRUE, stringsAsFactors=FALSE, na.strings="", sep=",", dec=".", skip=0)
listofStockConcentrations <- read.csv2(file=file.path("inst/extdata/library/listofstockconcentrations.csv"), check.names=FALSE, header=TRUE, stringsAsFactors=FALSE, na.strings="", sep=",", dec=".", skip=0)
```

**Importing the plate map of the drug library from a source plate**

Additionally, we need to import the plate map of the 'source' plate(s), which holds the drugs at their stock concentrations and serves as a drug library from which the drugs will be dispensed. The plate map of the source plate contains the drugs and their concentrations along with the coordinates and names of all wells of the plate from which drugs will be dispensed.

The plate map of the source plate will be used to create instructions for the pipetting robot from what well on the source plate the drugs have to be taken in order to be dispensed onto the 'destination' plate, which is the final plate on which the drug treatments are being run.

```{r, eval = FALSE}
sourcePlate <- read.csv2(file=file.path("inst/extdata/library/sourcePlate.csv"), check.names=FALSE, header=TRUE, stringsAsFactors=FALSE, na.strings="", sep=",", dec=".", skip=0)
```

```{r, eval = FALSE}
sourcePlate
```

![](https://raw.githubusercontent.com/Enserink-lab/screenwerk-2.0/refs/heads/develop/doc/figures/sourceplate-drugs.png)

![](https://raw.githubusercontent.com/Enserink-lab/screenwerk-2.0/refs/heads/develop/doc/figures/sourceplate-ctrls.png)

The use of a source plate is optional and only required, if the drug library is created and maintained by the user. In that case the custom drug library is either being used directly for dispensing of the drugs in-house or shipped to a third-party facility that will use the custom drug library for dispensing of the drugs based on the user generated instructions.

In any other case, the set up and generation of dispensing files will also work without the use of a source plate, however in that case the instructions from where the drugs are being dispensed will be missing. This might be the case, if the dispensing of the drugs is carried out and handled by a third-party or a dispensing facility that usually takes care of the handling of the drug library and the pipetting of the drugs.

```         
⟁ Note: In theany case, where no source plate is being used and imported, the columns, such as 'Source.Plate.Barcode' and 'Source.Well', will remain empty in the dipsensing file.
```

Alternatively, it is possible to import the plate map of the source plate in the format of a plate, rather than a list, either through the proprietary *.PlateMap* file format or through one or more *.csv* files.

**Importing the plate map from .PlateMap files**

The custom plate map can also be imported from *.PlateMap* files with the following function:

```{r, eval = FALSE}
sourcePlate <- importPlateMap("inst/extdata/library/C008.PlateMap", .sourcePlateConv = TRUE)
```

With the argument .*`sourcePlateConv`* the plate map can be formatted to the requirements of the source plate format, in which empty wells as well as components that are not compounds, are being removed. In any other case the entire plate map is returned as a data set.

In the example above, the plate map of a single source plate is being imported. If the drug library consists of multiple plates, it is possible to specify the folder from which multiple plate maps will be read.

```{r, eval = FALSE}
sourcePlate <- importPlateMap("inst/extdata/library/", .fileFormat = ".PlateMap", .sourcePlateConv = TRUE)
```

With the argument *`.fileFormat`* the format of the file from which the plate map will be read can be specified.

![](https://raw.githubusercontent.com/Enserink-lab/screenwerk-2.0/refs/heads/develop/doc/figures/M006.png)

![](https://raw.githubusercontent.com/Enserink-lab/screenwerk-2.0/refs/heads/develop/doc/figures/C008.png)

**Importing the plate map from .csv files**

The plate map can also be imported manually as a .csv file. This can be also done by specifying a single file

```{r, eval = FALSE}
sourcePlate <- importPlateMap("inst/extdata/library/drugs.csv", .fileFormat = ".csv", .sourcePlateConv = TRUE)
```

or by specifying the folder location from which multiple csv files will be imported.

```{r, eval = FALSE}
sourcePlate <- importPlateMap("inst/extdata/library/", .fileFormat = ".csv", .sourcePlateConv = TRUE)
```

```         
⟁ Note: If compounds from multiple plates are being imported, each .csv file containing the plate map must be labeled denoting the type of content of the plate map and a numeric suffix, such as drugs-1.csv, drugs-2.csv, concentrations-1.csv and concentrations-2.csv.
```

Please note that this function allows the import of plate maps with multiple components, such as *compounds*, *cells* and/or *other*. However, doing so, requires a specific file nomenclature, if the plate maps are being imported from individual .csv files. Those csv files follow a format represented by the plate layout. The first row denotes the column labels (numeric) for a given plate, while the first column denotes the row labels (alphabetic) for a given plate. The actual values for compounds (drugs) and the corresponding concentrations, and/or cells are then provided for each well following the layout of the plate type. The plate map for compounds, compound concentrations and cells needs to come from a separate file each. If multiple compounds are found in a given well, a single file for each compound needs to be provided.

The .csv file with the plate map of the drugs could look like

![](https://raw.githubusercontent.com/Enserink-lab/screenwerk-2.0/refs/heads/develop/doc/figures/sourceplate-drugs-csv.png)

while the corresponding drug concentrations

![](https://raw.githubusercontent.com/Enserink-lab/screenwerk-2.0/refs/heads/develop/doc/figures/sourceplate-drugconcentrations-csv.png)

**Group of components across multiple source plates**

The *.plateMapGrp* allows to group individual components where necessary. A case scenario would be, in which multiple cell lines need to be imported that share the same drug treatment, but are found across different plates, or in which multiple cell lines are found across a single plate. The benefit of grouping is that it does not require to import the same compounds multiple times for each cell line, but rather will replicate components that are not grouped for each grouped component.

```{r, eval = FALSE}
# Grouping of cell lines across single and multiple plates
# Note: First two cell lines are grouped on the same plate, while the last two cell lines on another plate
importPlateMap("inst/extdata/library/", .fileFormat = ".csv", 
               .plateMapGrp = list(type = "Cells", name = c("A375", "WM1366", "WM45.1", "FEMXV"), group = c("A", "A", "B", "B")))
```

For example: (a) If the same drugs and drug concentrations are to be used for multiple cell lines, the cells can be grouped, while the drugs left ungrouped. This will replicate the same drug treatment for each of the cell lines. (b) If multiple cell lines are found across a single plate. in case where multiple cell lines are spread across a number of plates, all cell lines found on an individual plate can be assigned to the same group. The drug treatment will then be replicated across that group of cell lines. Please note that only one component can be grouped at a time.

The grouping can be extended for the generation of plate barcodes. Instead of assigning individual grouping labels, a unique barcode can be assigned to each group, denoting the plate map to a single plate.

```{r, eval = FALSE}
# Extending labels with plate barcodes, same set two plates
importPlateMap("inst/extdata/library/", .fileFormat = ".csv", 
               .plateMapGrp = list(type = "Cells", name = c("A375", "WM1366", "WM45.1", "FEMXV"), group = c("0921A1", "0921A1", "0921A2", "0921A2")))

```

```{r, eval = FALSE}
# Extending labels with plate barcodes, two different sets of plates
importPlateMap("path/to/folder/", .fileFormat = ".csv", 
               .plateMapGrp = list(type = "Cells", name = c("A375", "WM1366", "WM45.1", "FEMXV"), group = c("0921A1", "0921A1", "0921B1", "0921B1")))
```

This does not apply to plate maps imported from .PlateMap files. Due to the internal structure of the file format, multiple components can be imported simultaneously through a single file. However, the name of the .PlateMap file will denote the plate id, or plate barcode.

**Generating a dispensing data set**

Once we have imported all files required, we can define the last arguments

```{r, eval = FALSE}
.ctrlReplicates = 8
.addUntreated = list(name = "Untreated", replicates = 8)
.finalWellVolume = list(volume = 5, unit = "μl")
.plateFormat = 1536
.destinationPlateID = "PID3-0726"
```

and run the following function:

```{r, eval = FALSE}
dispensingData <- generateDispensingData(listofCombinations, listofDrugs, listofDoses, listofVolumes, listofCtrls, listofStockConcentrations, sourcePlate,
listofExWells, .ctrlReplicates = 8, .addUntreated = list(name = "Untreated", replicates = 8), .finalWellVolume = list(volume = 5, unit = "μl"), .plateFormat = 1536, .destinationPlateID = "PID3-0726", .randomizeDispensing = TRUE, .backfilling = TRUE, .probeDispensing = FALSE)
```

This function will generate a dispensing data set with a map that will serve as a blueprint for the experimental drug screen. It distributes the drug treatments onto a given plate design.

Since the use of the `sourcePlate` is optional, it can be set to *FALSE*, which will still generate the dispensing data, however without any information from where the drugs should be dispensed from.

The option `.ctrlReplicates` allows to define how many replicates for each of the controls provided in `listofCtrls` should be generated. The higher the number of replicates, the better the approximation of the true response for the negative and positive controls will be.

In addition, it is possible to add an untreated control to the screen design, in which no drug or any other chemical compound will be added. This can be done with the argument `.addUntreated` by providing a list containing an element with the `name` for the untreated controls (i.e. 'Untreated') and the number of `replicates`. If no untreated control is needed, this argument can be left out or explicitly set to *FALSE*.

In order for the function to be able to map the stock concetrations to the final drug concentrations for each drug treatment the final volume across all wells found on the destination (treatment) plate. This volume can be provided with the argument `.finalWellVolume` as a list containing the volume and the unit.

The `.plateFormat` represents the type of plate format, which is given by the total number of wells for the used plates on which the screen will be conducted. The argument accepts a predefined numeric value designating the plate format based on either 6, 12, 24, 48, 96, 384 or 1536 wells.

The `.destinationPlateID` servers as an unique plate identifier, which can also be used as a plate barcode.

The function allows for the randomization of the drug treatments across all plates. This can be achieved by setting the argument `.randomizeDispensing` to *TRUE*. If *FALSE*, drug treatments will be dispensed in sequential order. The default is set to *TRUE*.

With the argument `.backfilling` it is possible to obtain the same solvent concentration across all the wells, which is usually the highest concentration of combined solvents used in a treatment where two or more drugs are combined. This argument is useful ensuring that the same conditions are present in each well and that the differences in response to individual drug treatments is not due to different solvent concentrations. By setting this argument to *TRUE*, each well in which only one unit of solvent is being dispensed, such as for single drug treatments, another unit of the solvent will be backfilled to the well. This way, all th wells will have the same base condition by having the same solvent concentration.

**Probe dispensing**

Another useful argument is `.probeDispensing` with the possibility to probe the dispensing without actually generating any dispensing data. This resembles a "dry-run", in which the function is being executed as if it was being run. However, instead of generating any dispensing data, a concise summary is displayed for a potential screen set-up, by providing the same data as one would do for an actual set-up. The summary will provide information about the total number of drug treatments, the unique number of drug combinations and single drug treatments, the number of excluded wells as well as the number of controls per plate, and the total number of plates for a given dispensing set-up. Particularly being able to estimate the number of plates for a given experimental design can be essential in the planning of a drug screen, without having to generate data for each possible screen design.

```{r, eval = FALSE}
generateDispensingData(listofCombinations, listofDrugs, listofDoses, listofVolumes, listofCtrls, listofStockConcentrations, sourcePlate, listofExWells, .ctrlReplicates = 8, .addUntreated = list(name = "Untreated", replicates = 8), .finalWellVolume = list(volume = 5, unit = "μl"), .plateFormat = 1536, .destinationPlateID = "PID3-0726", .randomizeDispensing = TRUE, .backfilling = TRUE, .probeDispensing = TRUE)
```

Alternatively the function can be run with a partial set of files for probing by either specifying individual files

```{r, eval = FALSE}
generateDispensingData(listofCombinations = listofCombinations, listofDrugs = listofDrugs, listofCtrls = listofCtrls, listofExWells = listofExWells, .ctrlReplicates = 8, .addUntreated = list(name = "Untreated", replicates = 8), .plateFormat = 1536, .probeDispensing = TRUE)
```

or by assigning to non-relevant files a NULL object.

```{r, eval = FALSE}
generateDispensingData(listofCombinations, listofDrugs, NULL, NULL, listofCtrls, NULL, NULL, listofExWells, .ctrlReplicates = 8, .addUntreated = list(name = "Untreated", replicates = 8), .finalWellVolume = list(volume = 5, unit = "μl"), .plateFormat = 1536, .destinationPlateID = "PID3-0726", .randomizeDispensing = TRUE, .backfilling = TRUE, .probeDispensing = TRUE)
```

**Obtaining a summary of the dispensing data**

The summary can also be obtained after a drug screen has been set-up by screenwerk through the function `summary(object, ...)` on the dispensing data.

```{r, eval = FALSE}
summary(dispensingData)
```

This will generate a summary of the dispensing as the probing would do.

> ```         
> Summary for dispensing ID: PID3-0726
>
> Number of drug treatments: 3 456
> Number of unique drug combinations: 2 430
> Number of unique single drug treatments: 342
> Number of excluded wells per plate: 156
> Number of controls per plate: 40
> Number of total plates: 3
> Number of drugs: 57,  Number of doses: 6/5
> ```

**Extracting the dispensing data**

The dispensing data is usually saved by screenwerk in an object, along with the original data, all the data lists and the summary. In order the extract the output from the dispensing object, the function `print(object, ...)` can be used.

```{r, eval = FALSE}
print(dispensingData)
```

This output can then be assigned to a variable or exported as a .csv file.

```{r, eval = FALSE}
finalDispensingData <- print(dispensingData)

write.csv2(finalDispensingData, file = file.path(exportDirectory, paste("Dispensing Data", "_", .destinationPlateID, "_V1", ".csv", sep = "")), row.names = FALSE, quote = FALSE)
```

**Export the dispensing data to dispensing** **files**

In order to export the dispensing data to files containing only the instructions for the dispensing robot, the function `save(object, .saveto, ..., .sets, .labels, .split = FALSE, .by, .format)` can be used.

```{r, eval = FALSE}
save(dispensingData, .saveto=file.path("../myDispensing/files"), .sets = 2, .labels = "alphabetic", .split = FALSE, .format = "Echo")

```

This function allows to generate multiple identical sets from a single dispensing layout, which is often a requirement whenever a drug screen is being performed on multiple samples, e.g. cell lines. In certain instances, an experiment requires to be repeated, either to obtain experimental replicates, or due to a failed experiment. In either case, a single dispensing layout can be replicated based to an unlimited number of sets as needed.

The number of sets can be specified with the argument `.sets`.

The labels can be assigned to each dispensing set either through a predefined selection of labels by using either a set of alphabetic or numeric values. This can be accomplished by setting `.labels` = "alphabetic", or `.labels` = "numeric", respectively. Otherwise, the labels have to be explicitly stated for each individual set.

Individual labels can be assigned as shown below.

```{r, eval = FALSE}
save(dispensingData, .saveto = "../myDispensing/files", .sets = 3, .labels = c("REP1", "REP2", "CTRL"), .split = FALSE, .format = "Echo")

```

In that case, the number of `.labels` have to match the number of `.sets`.

Please note, that if the number of sets is not provided, or explicitly set as `FALSE`, any provided labels will be ignored. If the argument is set as TRUE, instead of a numeric value, the number of sets will depend on the number of labels provided. In that case, for each label a set will be created. If the number of sets is large than the number of labels provided, the function will result in an error, while if the number of sets is smaller than the number of labels provided, only the first *n* labels will be used.

**Splitting the dispensing data to individual dispensing files**

Some dispensing robots require individual files as input for each source plate they are pipetting from. For that case, the function can split the dispensing data into individual files for each source plate with the argument `.split` set to `TRUE`. The default is `FALSE`, in which all sets are saved to one single file.

```{r, eval = FALSE}
save(dispensingData, .saveto=file.path("../myDispensing/files"), .sets = 2, .labels = "alphabetic", .split = TRUE, .format = "Echo")
```

Please note, that if splitting of the dispensing data is requested by setting the argument `.split` to `TRUE` and not specifying `.by` what group the data should be split, the *default* splitting parameter will be used, which is splitting by source plate. Otherwise the data can be split either `.by` = `"sourceplate`", by `"set"` or by `"solvent"` in which the drug has been dissolved.

```{r, eval = FALSE}
save(dispensingData, .saveto=file.path("../myDispensing/files"), .sets = 2, .labels = "alphabetic", .split = TRUE, by = "sourceplate", .format = "Echo")
```

The dispensing data can only be split by source plate, if one or more source plates have been used during the set-up of the drug screen and the generation of the dispensing data, see section *Generating a dispensing data set*. If no source plates have been used, the function will return a warning and export the data as a single file instead.

Alternatively, the data can be split by the number of sets provided or by the solvents used for each drug.

```{r, eval = FALSE}
save(dispensingData, .saveto=file.path("../myDispensing/files"), .sets = 2, .labels = "alphabetic", .split = TRUE, by = "set", .format = "Echo")
```

In case where no sets have been provided, but for whatever reasons a splitting by sets has been requested, the function will return a warning, ignore the splitting argument and export the data as a single file instead.

```{r, eval = FALSE}
save(dispensingData, .saveto=file.path("../myDispensing/files"), .sets = 2, .labels = "alphabetic", .split = TRUE, by = "solvent", .format = "Echo")
```

If the data is split by solvent, then for each solvent that was used to dissolve the drugs into, a separate file will be created containing only the dispensing instructions for those drugs that have been dissolved in a given solvent. This splitting method is particularly useful, if drugs dissolved in certain solvents need to be kept separated or on separate source plates due to different properties, such as evaporation potential.

**Saving dispensing files in dispensing machine dependent formats**

Certain dispensing robots / machines have specific requirements how the dispensing instruction should look like. This specific formatting can be achieved with the argument `.format`, which currently supports the Echo Acoustic Liquid Handlers from Labcyte Inc., Beckman Coulter Life Sciences. The formatting specific for this robot can be obtained by specifying the instrument name (E5XX-1366) or model (Echo 550).

```{r, eval = FALSE}
save(dispensingData, .saveto=file.path("../myDispensing/files"), .sets = 2, .labels = "alphabetic", .split = FALSE, .format = "Echo")
```

```{r, eval = FALSE}
save(dispensingData, .saveto=file.path("../myDispensing/files"), .sets = 2, .labels = "alphabetic", .split = FALSE, .format = "Echo 550")
```

```{r, eval = FALSE}
save(dispensingData, .saveto=file.path("../myDispensing/files"), .sets = 2, .labels = "alphabetic", .split = FALSE, .format = "E5XX-1366")
```

Alternatively, the complete dispensing data can be saved, which allows machines not listed or currently supported in this package to pick individual instructions for dispensing. This is accomplished by setting the argument `.format` = `"full"` or `.format` = `"complete"`.

```         
⟁ Note: The dispensing files will be save as a comma-separated values (csv) file, in compliance with the standard outlined in RFC 4180 (2005). For more information, see https://www.ietf.org/rfc/rfc4180.txt
```

**Plotting the plate map / layout from the dispensing data**

The plate map and experimental layout can be plotted from the generated dispensing data. This can be particularly helpful in visualizing the distribution of the drug treatments across the plates.

For that the generic plot function can be used on the dispensing data (as long as the data is an objects of class *dispensingData*).

```{r, eval = FALSE}
plot(dispensingData, .saveto=file.path("../myDispensing/plots"))
```

This should produce an image as shown below:

![](https://raw.githubusercontent.com/Enserink-lab/screenwerk-2.0/refs/heads/develop/doc/export/plates/PID3-0726-Plate1.png)

![](https://raw.githubusercontent.com/Enserink-lab/screenwerk-2.0/refs/heads/develop/doc/export/plates/PID3-0726-Plate2.png)

![](https://raw.githubusercontent.com/Enserink-lab/screenwerk-2.0/refs/heads/develop/doc/export/plates/PID3-0726-Plate3.png)

![](https://raw.githubusercontent.com/Enserink-lab/screenwerk-2.0/refs/heads/develop/doc/export/plates/PID3-0826-Plate3.png)

**Saving the dispensing data to a file**

It is recommended that the dispensing data is saved or archived as an R object to a file for reference or as a backup.

```{r, eval = FALSE}
saveRDS(dispensingData, file = file.path("../myDispensing/data/dispensingData.rds"))

```

That way the dispensing data cab always be loaded back into R and used for generating further sets at a later point in time, if needed.

### ***(b) post-experimental: data processing and data analysis***

#### **Reading raw data and consolidating data sets**

Once the drug screen has been run and all the screen plates have been read, we can use the function below to read all raw data files, in order to be able to build a working data set by consolidating the dispensing data with the raw measurements.

The raw measurements can be read with

```{r, eval = FALSE}
rfs <- readRAWData(.readfrom = dataDirectory, .fileformat = c(".csv", ".txt"), .format = "EnVision")
```

The function `readRAWData` is used to read raw measurement files from selected plate readers. It supports the most common text-based file types (.txt, .csv) in various export-formats. Those formats are machine specific and follow a proprietary layout.

The function is capable of detecting the export-format and identifying the used text delimiter, such as comma-, semicolon- or tab-separated, as well as the plate format and find the raw data of interest. For more information, see package documentation.

The format of the read-out files can be specified with the argument `.format`. This argument is machine specific, since each manufacturer and plate reader will use a different, sometimes proprietary, format.

As of now, files from the following machines are supported: **EnVision** and **VICTOR X** Multimode Plate Reader from PerkinElmer, Inc.

The specific format can then be selected with either `.format` = `"EnVision"` or with `.format` = `"VICTOR"`.

Alternatively, the data can be provided in a raw data format, independent of the export-format and the plate reader used. In that case the raw measurements are provided as a plate layout.

```{r, eval = FALSE}
rfs <- readRAWData(.readfrom = dataDirectory, .fileformat = c(".csv", ".txt"), .format = "raw")
```

With the raw measurements read, we can now build the final reference data set, used for downstream analysis of the drug sensitivity screen. This can be achieved by consolidating the raw measurements with the dispensing data.

Before we can do that, we need to import a barcode reference list with the names of the samples used in the drug screen and by associating them to the corresponding plate id and set.

```{r, eval = FALSE}
.barcodeReference <- read.csv(file=file.path(libDirectory, "platebarcode.csv"), check.names=FALSE, header=TRUE, stringsAsFactors=FALSE, colClasses=c("PlateID"="character"), comment.char = "#", blank.lines.skip	= TRUE, na.strings="", sep=",", dec=".", skip=0)
```

The imported barcode reference could look like shown below.

![](https://raw.githubusercontent.com/Enserink-lab/screenwerk-2.0/refs/heads/develop/doc/figures/barcode-reference.png)

> ```         
>     PlateID Set Number     Sample
> 1 PID3-0726   A      3 OVCAR-8 R1
> 2 PID3-0726   B      3 OVCAR-8 R2
> 3 PID3-0726   C      3 OVCAR-8 R3
> ```

**Consolidating the raw measurements with the dispensing data**

Now we can consolidate the raw measurements with the dispensing data, using the function below:

```{r, eval = FALSE}
clData <- consolidateData(dispensingData = dispensingData, rawMeasurements = rfs, .barcodeReference)
```

**Exporting the consolidated data**

Once the data has been consolidated and a reference data set with the raw measurements build, the data set can be exported at this stage for the analysis with third-party tools.

This can be achieved by using again *save* function, as shown below:

```{r, eval = FALSE}
save(clData, .saveto = file.path("../myExportFiles/files", "consolidatedData.csv"), .fileformat = ".csv", .sep = ";")
```

With the argument `.format` the data can be exported for specific analytical applications in the required format. For now, the list of supported formats include: BREEZE. The file type can be set with the argument `.fileformat`, which will set the file extension accordingly. In addition, the field separator (deliminator) can be set with the argument `.sep`.

```{r, eval = FALSE}
save(clData, .saveto = file.path("../myExportFiles/files", "consolidatedData.csv"), .fileformat = ".csv", .format = "breeze")
```

#### Running the quality control (QC) analysis

Now that we have the raw measurements annotated, we can use them to run the quality control (QC).

The function `qc` is used to assess the quality of a drug sensitivity screen by looking at the variance and signal distribution between individual controls.

```{r, eval = FALSE}
qcdata <- qc(consolidatedData, .saveto = file.path(resultDirectory), .ctrls=list(positive = "BzCl", negative = "all"), .qcMethod = "all")
```

The QC module consists of different quality control methods, which all offer a different approach to the assessment of the data and consequently offer a different perspective on the quality of the data.

**variance** : assessing the variance between individual controls both, across all plates, as well as by individual plate **emptywells** : assessing the signal of empty and untreated wells, this will also include any excluded wells **firstcolumn** : assessing the signal of wells in the first column of each plate **zprime** : assessing the z'-factor based on the distribution between the positive and negative control.

The z'-factor is a metric providing a measure of quality for high-throughput screens. The calculations are based on the z'-factor described in the original paper by Zhang, 1999, J Biomol Screen (see references).

It is possible to run multiple QC methods at once, or **all** by simply specifying `.qcMethod = "all"`.

With the parameter *.ctrls*, it is possible to specify a list with a set of positive and negative controls for (QC) analysis. The **positive** control needs to be provided with a single positive control, while the **negative** control can be a set with any number of negative controls. Alternative the **negative** control can be set with **all** leading to the inclusion of all negative controls present in the drug sensitivity screen.

**Assessing the variance between individual controls**

```{r, eval = FALSE}
# Assessing the variance between individual controls
qc(consolidatedData, .saveto = "path/to/folder/", .ctrls = list(positive="BzCl", negative=c("DMSO", "H2O", "TFA", "Untreated")), .qcMethod = "variance")
```

**Assessing the signal of empty wells across plates**

```{r, eval = FALSE}
# Assessing the signal of empty wells across plates
qc(consolidatedData, .saveto = "path/to/folder/", .qcMethod = "emptywells")
```

**Assessing the z'-factor for a given set of controls**

```{r, eval = FALSE}
# Assessing the z'-factor for individual plates by looking at the distribution between the positive and negative controls
qc(consolidatedData, .saveto = "path/to/folder/", .ctrls = list(positive="BzCl", negative="DMSO"), .qcMethod = "zprime")
```

**Running multiple QC assessments**

```{r, eval = FALSE}
# Running multiple QC assessments
qc(consolidatedData, .saveto = "path/to/folder/", .ctrls = list(positive="BzCl", negative="DMSO"), .qcMethod = c("variance", "emptywells", "firstcolumn"))
```

**Running all QC assessments with selected controls**

```{r, eval = FALSE}
# Running all QC assessments with selected controls
# Note: If zprime is being used, only the first two controls will be used as the positive and negative control. Set the order of controls in that case accordingly.
qc(consolidatedData, .saveto = "path/to/folder/", .ctrls = list(positive="BzCl", negative="DMSO"), .qcMethod = "all")
```

**Running all QC assessments with selected (negative) controls**

```{r, eval = FALSE}
# Running all QC assessments with multiple (negative) controls
# Note: If zprime is being used, only the first two controls will be used as the positive and negative control. Set the order of controls in that case accordingly.
qc(consolidatedData, .saveto = "path/to/folder/", .ctrls = list(positive="BzCl", negative=c("DMSO", "H2O", "DMSO+H2O", "TFA", "Untreated")), .qcMethod = "all")

```

**Running all QC assessments with all (negative) controls**

```{r, eval = FALSE}
# Running all QC assessments with all (negative) controls
# Note: If zprime is being used, only the first two controls will be used as the positive and negative control. Set the order of controls in that case accordingly.
qc(consolidatedData, .saveto = "path/to/folder/", .ctrls = list(positive="BzCl", negative="all"), .qcMethod = "all")
```

Any plots generated during the QC analysis will be saved to the location provided with `.saveto`.

**Plotting a plate heatmap of the signal distribution across experimental / control plates**

In some instances the experimental set-up has additional QC measures implemented, such as the use of control plates. Those plates are generally used to assess the quality of the dispensing, such as the dispensing of cells, across multiple plates or timepoints.

In order to assess the signal on those plates we can plot the raw measurements and compare the signal across the plate, as well as between individual plates.

For that we can simply use the same function we did before for reading the raw measurements and use the generic `plot` function to plot a heatmap of the plates.

```{r, eval = FALSE}
# Reading raw data of control plates
ctrlPlates <- readRAWData(.readfrom, .fileformat = c(".csv", ".txt"), .format = "EnVision")
```

and now we can plot the data

```{r, eval = FALSE}
# Plotting heatmap of control plates
plot(ctrlPlates, .saveto = "path/to/folder/")
```

#### Data processing: normalization of the data

In order to be able to analyze the screen, it is necessary to normalize the data. Furthermore, the data needs to be processed for various tasks in the downstream analysis of the drug sensitivity screen. This module carries out three essential steps in data processing: normalizing, splitting and assembling (re-formatting) the data.

Raw measurements (CPS, counts/s) are normalized to median positive and median negative controls using the equation: (drug – positive control) / (negative control – positive control).

All this can be achieved with `processData`, which normalizes the raw measurements to the positive and negative control, and subsequently splits the data into individual data sets, one for the controls, single drug treatments and combination treatments, in case of a drug combination screen. Furthermore, a data table and a matrix is assembled for each drug pair, representing the dose response between two drugs at each dose.

```{r, eval = FALSE}
processedData <- processData(consolidatedData, .ctrls=list(positive="BzCl"))
```

The controls are essential for the normalization of the raw measurements to the positive and negative control. The positive control has to be specified, while the negative control will be retrieved from the solvent in which the drug has been dissolved. This requires an individual set of controls for each solvent and combination of solvents.

```         
⟁ Note: Since drug treatments, such as single drug treatments, are backfilled with one unit of the control corresponding to the solvent in which the drug has been dissolved, each well has the highest concentration of solvent used across all the treatments in the drug screen.
```

#### Single drug dose response modeling (DRM)

Now that we have normalized and formatted the data, we can use the data to run the first analytic module of screenwerk, which is the assessment of the single drug treatments.

The function *`runDRM`* estimates the dose-response of the single drug treatments. It performs curve fitting using a *four-parameter log-logistic function (LL.4)*, and estimates the EC10, EC50 and EC90. The function also allows to plot the single drug curves by viability and inhibition for each single drug.

```{r, eval = FALSE}
doseRespModel <- runDRM(processedData, .saveto = "path/to/folder/", .plot = TRUE)
```

The function is used to run the first of a series of analyses, specifically, it looks at the dose response of each single drug. It performs curve fitting based on a four-parameter log-logistic model and estimates the EC50s, along with the EC10 and EC90. The curve fitting is performed with lower limits set to 0 and the upper limit to 1. In cases where the curve fitting fails, such as due to lack of response resulting in a flat line, or due to outliers, or any other irregularities, the fitting is performed without restrictions to either one of the limits.

All estimates of the EC50s is exported as a table in a csv file format. Plots are exported as composite of all drugs for each sample as individual png image files.

Files are saved either to the specified location or the default working environment, with the corresponding folder structure: 'results/EC50s' and 'results/graphs/single drug response'

The dose-response analysis is performed with the function **drm** from the R-package *drc* as described in the original paper by Ritz, 2015, PLOS ONE (see references).

**Additional plotting of dose response from single drug treatments**

This module is a complementary component of screenwerk providing a set of custom plots for the visualization of drug-dose responses. `customPlotting` is a function that generates a set of plots that provide an overview of the dose-response between drugs and samples. Furthermore, it plots the dose-response matrix for drug combinations along with the dose-response curves for each individual drug pair.

```{r, eval = FALSE}
customPlotting(processedData, doseRespModel, .saveto = "path/to/folder/")
```

This module is meant as a library that contains a set of custom plot that have been requested over time by different users. Additional requests for custom plots can be added to this function over time. This function is meant to be used for the generation of additional plots to the DRM, in order to provide an overview of the responses of all treatments between individual drugs and samples. This is of particularly benefit for large drug screens in which a large number of drugs and samples have been screened.

Files are saved either to the specified location or the default working environment, with the corresponding folder structure: 'results/graphs/custom plots'

#### Estimating the dynamic range of drug activity based on the single drug response

This is another complementary module of screenwerk assessing the dynamic drug-activity range for individual drug-dose responses. `dynamicRange` is a function that estimates the dynamic drug-activity range (DDAR) across a number of doses for each drug response and generates a set of plots.

```{r, eval = FALSE}
ECDR <- dynamicRange(doseRespModel, .saveto = "path/to/folder/")
```

The function can use one of two data sets, either an object of class '`processedData`' or an object of class '`drm`'.

```{r, eval = FALSE}
ECDR <- dynamicRange(processedData, .saveto = "path/to/folder/")
```

The function `dynamicRange` is used to provide an assessment of the drug activity range across the selected range of doses. It is an essential indicator, especially for drug combination screens, in which only a selected range of doses are combined to assess synergies. If doses are combined at which drugs enfold their full inhibitory potential, it won't leave enough room for potential combinatory effects.

The dynamic range is considered the dose-response range between the ED10 and ED90. This function will indicate the expected and the observed drug activity range for each drug and sample. It will generate plots based on the fitted dose-response models as well as unfitted curves.

Files are saved either to the specified location or the default working environment, with the corresponding folder structure: 'results/Dynamic Drug Activity Range' and 'results/graphs/dynamic range'

#### Plotting kinetic dose response curves

This is yet another complementary module of screenwerk providing a set of kinetic plots for the visualization of drug-dose responses.

```{r, eval = FALSE}
kineticPlots(doseRespModel, .saveto = "path/to/folder/")
```

The function can use one of two data sets, either an object of class '`processedData`' or an object of class '`drm`'.

```{r, eval = FALSE}
kineticPlots(processedData, .saveto = "path/to/folder/")
```

The function `kineticPlots` is used to provide an overview of the **raw** and **unfitted** responses of all treatments between individual drugs and samples. This is of particularly benefit for large drug screens in which a large number of drugs and samples have been screened.

Files are saved either to the specified location or the default working environment, with the corresponding folder structure: 'results/graphs/single drug response'

#### Plotting waterfall plots

Another complementary module of screenwerk providing a set of waterfall plots for the visualization of single drug treatment responses.

```{r, eval = FALSE}
waterfallPlots(doseRespModel, .saveto = "path/to/folder/")
```

The function can use one of two data sets, either an object of class '`processedData`' or an object of class '`drm`'.

```{r, eval = FALSE}
waterfallPlots(processedData, .saveto = "path/to/folder/")
```

The function `waterfallPlots` is used to provide an overview of the responses of all treatments between individual drugs and samples as waterfall plots based on the area under the curve (AUC).

Files are saved either to the specified location or the default working environment, with the corresponding folder structure: 'results/graphs/waterfall plot'

#### Plotting heatmap of single drug responses

A complementary module of screenwerk providing a heatmap plot for the visualization of single drug treatment responses.

```{r, eval = FALSE}
heatmapPlots(doseRespModel, .saveto = "path/to/folder/")
```

The function can use one of two data sets, either an object of class '`processedData`' or an object of class '`drm`'.

```{r, eval = FALSE}
heatmapPlots(processedData, .export = c("dendrogram", "rownames"), .saveto = "path/to/folder/")
```

The function `heatmapPlots` is used to provide an overview of the responses for all drugs and samples as a heatmap based on the area under the curve (AUC).

With the argument `.export` a set of predefined heatmap annotation (such as dendrogram and row names) can be plotted separately in addition to the heatmap.

Files are saved either to the specified location or the default working environment, with the corresponding folder structure: 'results/graphs/heatmaps'

#### Drug interactions assessment from drug combination treatments

After the single drug treatments have been analyzed and assessed, the anlysis can be continued with the drug combination treatments. This can be done with the function `bayesynergy`.

```{r, eval = FALSE}
bayesdata <- bayesynergy(processedData, .saveoutput = TRUE, .plot = TRUE, .saveto = "path/to/folder/")
```

The function `bayesynergy` is used to assess the interaction between two drugs across their dose ranges based on a bayesian semi-parametric model. The analysis is based on the package [`bayesynergy`](http://127.0.0.1:15262/help/library/bayesynergy/help/bayesynergy) as described in the original paper by Rønneberg, 2021, Brief Bioinform (see references).

With the argument `.saveoutput` = `TRUE` the function offers to save the output from `bayesynergy` to a file. The default is `FALSE`, by which the output is not saved. Furthermore, a set of plots can be generated with the argument `.plot` = `TRUE` and saved to the folder location specified with `.saveto`. The default is `FALSE`, by which no plots are generated.

```{r, eval = FALSE}
bayesdata <- bayesynergy(processedData)
```

The function returns an class S3 object of type list with the volumetric surfaces (VUS), the EC50s, summary statistics and additional quality parameters, along with a measure of synergy (bayesfactor).

Files are saved either to the specified location or the default working environment, with the corresponding folder structure: 'results/data/bayesynergy' and 'results/graphs/bayesynergy'

#### Scoring of drug interactions

Once the drug interactions have been estimated, the data can be used to score drug interactions and rank them based on their individual synergy scores.

```{r, eval = FALSE}
synMetrics <- synergyScoring(bayesdata, .saveoutput = TRUE, .plot = TRUE, .saveto = "path/to/folder/")
```

The function `synergyScoring` is used to score and rank drug interactions between drugs of a combination treatment. If plotting is set to `.plot` = `TRUE`, the function will generate three types of plots: (1) one in which all drug synergy and antagonism scores are plotted for each sample, (2) another one in which only the synergy score is plotted by drug fore each sample and (3) one in which only synergy and antagonism scores are plotted for each individual drug. Furthermore, the data output can be save to a file with `.saveoutput` = `TRUE`. The default is `FALSE` for both arguments, by which no plots will be generated and no data will be saved to a file.

```{r, eval = FALSE}
synMetrics <- synergyScoring(bayesdata)
```

Files are saved either to the specified location or the default working environment, with the corresponding folder structure: 'results/synergyscores/bayesynergy' and 'results/graphs/synergyscores/bayesynergy'

**References**

Ritz, C., Baty, F., Streibig, J. C., Gerhard, D. (2015) Dose-Response Analysis Using R. PLOS ONE 10 (12), e0146021. DOI: 10.1371/journal.pone.0146021

Rønneberg, L., Cremaschi, A., Hanes, R., Enserink, J. M., Zucknick, M. (2021) bayesynergy: flexible Bayesian modelling of synergistic interaction effects in in vitro drug combination experiments. Brief Bioinform 22 (6), bbab251, DOI: 10.1093/bib/bbab251

Zhang, J. H., Chung, T. D., Oldenburg, K. R. (1999): A Simple Statistical Parameter for Use in Evaluation and Validation of High Throughput Screening Assays. J Biomol Screening 4 (2), S. 67–73. DOI: 10.1177/108705719900400206
