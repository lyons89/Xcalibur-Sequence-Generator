# Thermo Xcalibur Sequence Generator

A Shiny app for building and downloading sample-sequence tables for Thermo mass spectrometers. It creates sample and pool names, inserts blank runs, assigns vial positions, and selects instrument methods based on the chosen instrument and acquisition method. The downloaded CSV is formatted for import into an Xcalibur workflow, with a `Bracket Type=4` line before the table.

> \*\*Important:\*\* This app was originally designed for use in the UNC Metabolomics and Proteomics Core Facility. 

> \*\*Important:\*\* Review each generated sequence and confirm that the methods, vial positions, injection volume, and run order match your instrument and laboratory SOP before using it. This app generates a sequence file; it does not connect to or control a mass spectrometer. 


## Features

* Supports the methods `Astral-DIA`, `Astral-DDA`, and `OE480-DDA`.
* Builds file names from the date, instrument, column, user-provided identifier, method, and sample number.
* Supports technical injections (replicates) for each sample.
* Supports three randomization modes:

  * **None** — retain the generated sample order.
  * **All** — shuffle all sample injections together.
  * **Block** — shuffle samples within conditions, then create successive blocks containing one sample per condition when available. Blank runs are inserted at block boundaries.
* Adds pool runs and blanks according to the requested pool count.
* Assigns vial positions using plate layouts defined in the script.
* Displays the generated sequence in a table and downloads it as a CSV file.
* Allows an injection volume and a comment to be applied to all rows.

## Requirements

The app uses R and these R packages:

* `tidyverse`
* `shiny`
* `DT`
* `markdown`

Install missing packages from an R console:

```r
install.packages(c("tidyverse", "shiny", "DT", "markdown"))
```

The instrument-method paths in `app.R` are Windows paths under `C:\\\\Xcalibur\\\\methods\\\\...`. They are written into the output file as method assignments; make sure these paths and method names correspond to the methods installed on the target instrument computer.

## Run the app

Save `app.R` in a working directory, open that directory as the R working directory, and run:

```r
shiny::runApp()
```

Alternatively, open `app.R` in RStudio and select **Run App**.

## How to use it

1. Enter a **Unique Identifier** to distinguish the run or project.
2. Select an acquisition **Method**.
3. Enter the number of **Samples**, **technical injections**, and **Pools**. A technical-injection value of `0` means no replicate suffixes; a pool count of `0` means no pool runs.
4. Enter the **Starting Vial Position** using the format expected for the selected instrument. Examples shown by the app are `R:B1` for Astral and `RB1` for OE480.
5. Set the injection volume and optional comment.
6. Choose a randomization mode. For **Block**, paste one condition per sample, in the same order as the samples were entered. Separate conditions with new lines.
7. Click **Generate Sequence** to display the table.
8. Optionally enter a file-name label, then click **Download Sequence** to save the CSV.

### Example operation

For example, to make an Astral-DIA sequence for four samples, one technical injection per sample, two pools, and block randomization:

* Method: `Astral-DIA`
* Unique Identifier: `Pilot01`
* Number of Samples: `4`
* Technical injections: `1`
* Number of Pools: `2`
* Starting Vial Position: `R:B1`
* Randomization: `Block`
* Conditions, one per line:

```text
Control
Treatment
Control
Treatment
```

After clicking **Generate Sequence**, inspect the table and download it. The exact order of samples in a randomized sequence will vary between generations.

## How it works

### Names and methods

The app uses the current date (`YYYYMMDD`), instrument, column, identifier, method suffix, and sample number to construct sample file names. When technical injections are requested, it adds replicate suffixes such as `\_R1` and `\_R2`. Pool names receive `\_Pool\_R1`, `\_Pool\_R2`, and so on. The selected acquisition method determines the sample method path, blank method path, and column label embedded in the sequence.

### Randomization and conditions

* **None** leaves samples in their generated order.
* **All** randomly permutes the sample injections.
* **Block** requires one condition for each original sample. Conditions are repeated for technical injections. The app shuffles samples within each condition, repeatedly takes one available sample from each condition to form a block, and shuffles the order within each block. If a condition has fewer samples, later blocks omit it once its samples are exhausted.

Conditions are ignored for **None** and **All**. Randomization is not seeded in the app, so a new generation may produce a different order.

### Blanks and pools

The sequence begins with two blanks. The app then places sample injections and pools according to the pool count:

* **0 pools:** samples, followed by two blanks.
* **1 pool:** one pool near the beginning, then samples, followed by two blanks.
* **2 pools:** one pool near the beginning and one near the end, with blanks around the pool placements.
* **3 or more pools:** the first pool is near the beginning, middle pools are spaced through the sample list, and the final pool is near the end.

For block randomization, a blank is also inserted between adjacent blocks. The code uses a fixed blank position for blank injections.

### Vial positions and CSV output

The script defines plate layouts separately for Astral and OE480 and rotates the ordered list of sample positions so assignment starts at the requested vial position. Sample vial assignment is based on the sample number embedded in its generated name; technical replicates of the same sample therefore share a vial position. Pool injections are assigned a position after the sample positions when the pool-position assignment is enabled by the code.

The download is a comma-separated file with a first line containing `Bracket Type=4`, followed by columns for file name, path, instrument method, position, injection volume, and comment. The download name is based on the current date and the optional **File Name** label.

## Implementation notes

* The app runs locally in R/Shiny and generates a file; it does not validate that the destination instrument has the configured methods or that a sequence is operationally safe.
* The current UI declares a minimum of `2` for the technical-injection and pool numeric inputs, although the sequence-generation code supports `0` for both. If your Shiny/browser setup prevents entering values below the declared minimum, change those two `numericInput()` minimums in `app.R` to `0`.
* The starting-position controls are implemented as `textInput()` elements. The method-change observer currently calls `updateSelectInput()` for those controls; if the default position does not update when switching instruments, replace those calls with `updateTextInput()`.
* The code assigns a fixed location to all blanks. Confirm this matches your laboratory's vial-loading convention.

